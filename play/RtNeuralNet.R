# Covariate candidates: Climate, Mask use, Mobility (Apple/Google)

# source('data/covariates.R')
# source('data/2.signals.R')
# source('data/variants.R')

b117.ts.partialGen <- function(region) {
  out <- tibble(date=c(), b117prop=c(), state=c())
  fit_tibble <- fit_growth(region)
  for (state_ab in data.variants.b117.region.states[[region]])
    out <- rbind(out, mutate(fit_tibble, state = state_ab))
  out
}
b117.ts.gen <- function() {
  out <- tibble()
  for (region in 1:10) {
    out <- rbind(out, b117.ts.partialGen(region))
  }
  out %>% rename(b117prop = prop)
}

covariates <- inner_join(cc.masked, data.mobility.apple, by = c('state', 'date')) %>%
  mutate(month = month(date)) %>%
  inner_join(data.weather.temp, by = c('state', 'month')) %>%
  inner_join(indicators.can, by = c('state', 'date')) %>%
  inner_join(b117.ts.gen(), by = c('state', 'date')) %>%
  select(state, date, masked, driving, walking, transit, temp, cases, b117prop) %>%
  group_by(state) %>% arrange(date) %>%
  mutate(severity = c(0, diff(cases)) / (data.population.states[state] / 100000)) %>%
  mutate(severity = rollmean(severity, k=7, na.pad=TRUE),
         driving = rollmean(driving, k=12, na.pad=TRUE),
         walking = rollmean(walking, k=12, na.pad=TRUE),
         transit = rollmean(transit, k=12, na.pad=TRUE)) %>%
  drop_na()

library(keras)
library(tidyverse)
library(caret)

# Model characteristics
window_len <- 25
num_epochs <- 10
atomic_horizon <- 1

set.seed(2424)

library(gridExtra)
source('deliver/PlotTheme.R')
.plot.covariates <- function(stateab) {
  joint.table <- covariates %>% filter(state == stateab, !is.na(masked)) %>%
    select(date, driving, masked, temp, transit, walking, b117prop, severity)
  
  # bigplot <- ggplot(joint.table, aes(date)) + geom_line(aes(y=value)) +
  #   facet_wrap(~name, scales = 'free_y') +
  #   ylab('') +
  #   poster.theme
  plot.mob <- ggplot(joint.table %>% pivot_longer(c('driving', 'walking', 'transit'), names_to = 'Transportation')) +
    geom_line(aes(date, value, color = Transportation)) +
    geom_hline(yintercept = 0, color = 'gray40') +
    ylab('Baseline % Difference') +
    ggtitle('Mobility Trends') + poster.theme

  plot.temp <- ggplot(joint.table %>% filter(day(date) == 15)) +
    geom_point(aes(date, temp)) +
    geom_smooth(aes(date, temp)) +
    ylab('Temperature (°F)') +
    ggtitle('Monthly Average Temperature') + poster.theme

  plot.masks <- ggplot(joint.table) +
    geom_line(aes(date, masked)) +
    ylab('Proportion of peers seen masked') +
    ylim(0.6, 1) +
    ggtitle('Mask use rating') + poster.theme

  plot.b117 <- ggplot(joint.table) +
    geom_line(aes(date, b117prop)) +
    ylab('Proportion of B.1.1.7') +
    ggtitle('B.1.1.7 Variant Propagation') + poster.theme


  grid.arrange(plot.mob, plot.temp, plot.masks, plot.b117, nrow=2, ncol=2)
}

# Inputs: Temperature, all 3 mobility metrics,
# population density (constant)
# mask usage percentage, past Rt
#  -> For each RNN iteration: 1 + 3 + 1 + 1 + 1: 7 parameters
# Outputs: Mobility metrics, mask usage, Rt
#  -> For each RNN iteration: 3 + 1 + 1: 5 outputs
size.input <- 9
size.output <- 6

popdensity.zscore.gen <- function(states) {
  all_pops <- data.population.density.states %>% keep(~ .x < 2000) %>% as.numeric
  the_mean <- mean(all_pops)
  the_sd <- sd(all_pops)
  to_zscore <- ~ (.x - the_mean) / the_sd
  states %>% map(~ data.population.density.states[.x]) %>%
    map(to_zscore) %>%
    reduce(c)
}

fittings <- readRDS('NewFits.rds')
fittedRts <- fittings %>% imap_dfr(~ mutate(.x[[3]]$history, state = .y)) %>%
  filter(date >= '2020-12-01') %>%
  group_by(state) %>% arrange(date) %>%
  mutate(Rt = rollmean(Rt, k=7, na.pad=TRUE)) %>%
  drop_na


popdensity.zscores <- popdensity.zscore.gen(unique(fittedRts$state))
  
fittedCovariates <- covariates %>%
  filter(state %in% names(fittings), date >= '2020-12-15') %>%
  inner_join(fittedRts %>% select(date, state, Rt), by = c('date', 'state')) %>%
  mutate(popdensity = sapply(state, function(x) popdensity.zscores[x]),
         temp = (temp - mean(temp)) / sd(temp),
         severity = (severity - mean(severity)) / sd(severity),
         driving = driving / 2 + 0.5,
         walking = walking / 2 + 0.5,
         transit = transit / 2 + 0.5,
         masked = (masked - min(masked)) / ((max(masked) - min(masked))))

# model <- keras_model_sequential() %>% 
#   layer_dense(input_shape = dim(X_train)[2:3], units = window_len) %>%
#   layer_simple_rnn(units = 10) %>%
#   layer_dense(units = 1, activation = 'relu')



# [samples, time_steps, features]
# [~600, 15, 9] to [~600, 5]
num_samples <- ((max(fittedCovariates$date) - min(fittedCovariates$date)) / ddays() - window_len + 2 - atomic_horizon) * length(unique(fittedCovariates$state))

weights <- c(2, 1, 1, 1, 2, 8)
weights <- weights / sum(weights)

model <- keras_model_sequential() %>%
  layer_gru(units = window_len * size.input * 2,
             input_shape = c(window_len, size.input),
            dropout = 0.1) %>%
  layer_dense(units = 20, activation = 'relu') %>%
  layer_dense(units = size.output * atomic_horizon, activation = 'linear') %>%
  layer_activation_leaky_relu() %>%
  compile(loss = 'mae', optimizer = 'adam', metrics = c('accuracy'),
          loss_weights = rep(weights, atomic_horizon))

X_all <- array(0, c(num_samples, window_len, size.input))
y_all <- array(0, c(num_samples, size.output * atomic_horizon))
.samples_filled <- 0
for (st8 in unique(fittedCovariates$state)) {
  all_ts <- fittedCovariates %>%
    filter(state == st8) %>%
    ungroup %>%
    arrange(date) %>%
    select(masked:temp, b117prop:popdensity) %>%
    as.data.frame %>% array
  numrows <- dim(all_ts)[1] - window_len - atomic_horizon + 1
  for (row.start in 1:numrows) {
    X_all[.samples_filled + 1,,] <- all_ts[row.start:(row.start+window_len-1),] %>% data.matrix
    y_all[.samples_filled + 1,] <- all_ts[(row.start + window_len):(row.start + window_len + atomic_horizon - 1),c(1:4, 7, 8)] %>%
      data.matrix %>% t %>% as.vector
    .samples_filled <- .samples_filled + 1
  }
}

# Randomly ablate 10% for the test set
cross_train <- function(proportion_out = 0.1, epochs = 40) {
  test_indexes <- sample(1:num_samples, round(num_samples * proportion_out))
  X_test <- X_all[test_indexes,,]
  y_test <- y_all[test_indexes,]
  X_train <- X_all[(1:num_samples)[-test_indexes],,]
  if (any(is.na(X_train))) {
    print(str_glue("NA found in {(1:num_samples)[-test_indexes]}"))
  }
  y_train <- y_all[(1:num_samples)[-test_indexes],]
  fitting <- fit(model, X_train, y_train, validation_data = list(X_test, y_test), epochs = epochs)
  list(fit = fitting, X_test = X_test, y_test = y_test, validation_data = list(X_test, y_test))
}

demonstrate_rt <- function(st8, start_date, selector = function(x) select(x, masked, severity, driving, walking, transit, Rt)) {
  ts_array <- fittedCovariates %>%
    filter(state == st8) %>%
    ungroup %>% filter(date >= start_date)
  output <- tibble()
  start_row <- 1
  while (start_date + days(start_row-1 + window_len) <= max(ts_array$date)) {
    X <- array(0, c(1, window_len, size.input))
    X[1,,] <- (ts_array %>%
                 select(masked:temp, b117prop:popdensity) %>%
                 as.data.frame %>% array)[start_row:(start_row + window_len - 1),] %>% data.matrix
    y <- predict(model, X) %>% as.numeric %>% matrix(nrow = size.output) %>% t
    colnames(y) <- c('masked', 'driving', 'walking', 'transit', 'severity', 'Rt')
    ydates <- seq(start_date + days(start_row-1 + window_len),
                  start_date + days(start_row-1 + window_len) + days(atomic_horizon - 1),
                  by = 'd')
    updater <- as_tibble(y) %>%
      selector %>%
      mutate(date = ydates)
    ts_array <- rows_upsert(ts_array, updater, by = 'date')
    colnames(y) <- c('P.masked', 'P.driving', 'P.walking', 'P.transit', 'P.severity', 'P.Rt')
    output <- output %>% rbind(as_tibble(y) %>%
                                 mutate(date = ydates))
    start_row <- start_row + atomic_horizon
  }
  output %>% inner_join(fittedCovariates %>% filter(state == st8))
}

# 0.1567
# 0.2232