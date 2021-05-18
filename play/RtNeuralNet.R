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
  mutate(severity = rollmean(severity, k=7, na.pad=TRUE))

library(keras)
library(tidyverse)
library(caret)

# Model characteristics
window_len <- 15
batch_size <- 20
num_epochs <- 10

set.seed(9913)

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
# population density (constant), number of new cases per 100,000 (pandemic awareness / response, which could be associated with mask use)

fittedRt <- list(
  CA = try.ca.final[[3]]$history %>% mutate(state = 'CA'),
  MA = try.ma.final[[3]]$history %>% mutate(state = 'MA'),
  MI = try.mi.final[[3]]$history %>% mutate(state = 'MI'),
  WA = try.wa.final[[3]]$history %>% mutate(state = 'WA'),
  FL = try.fl.final[[3]]$history %>% mutate(state = 'FL')
)
fittedRts <- reduce(fittedRt, rbind)

fittedCovariates <- covariates %>%
  filter(state %in% names(fittedRt), date >= '2021-01-01') %>%
  inner_join(fittedRts %>% select(date, state, Rt), by = c('date', 'state'))

covMatrix <- matrix(nrow = nrow(fittedCovariates) - window_len + 1, ncol = window_len + 1)
for (i in 1:nrow(fittedCovariates)) {
  fittedCovariates
}

model <- keras_model_sequential() %>% 
  layer_dense(input_shape = dim(X_train)[2:3], units = window_len) %>%
  layer_simple_rnn(units = 10) %>%
  layer_dense(units = 1, activation = 'relu')
