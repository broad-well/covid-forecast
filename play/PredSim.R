source('play/CompartmentSim.R')
source('play/VaccToImmunity.R')
source('play/RtNeuralNet.R')

prepare_pred_tibble <- function(st8, last_real) {
  ts_array <- fittedCovariates %>%
    filter(state == st8) %>%
    ungroup %>%
    arrange(date) %>%
    filter(date <= last_real) %>%
    tail(window_len) %>%
    select(date, masked:temp, b117prop:popdensity, state)
}

psdata.b117prop <- b117.ts.gen()

predict_next_day <- function(pred_tibble) {
  X <- array(0, dim = c(atomic_horizon, window_len, size.input))
  X[1:atomic_horizon,,] <- pred_tibble %>%
    arrange(date) %>%
    select(-date, -state) %>%
    tail(window_len) %>%
    data.matrix %>% array
  
  y <- predict(model, X)
  
  colnames(y) <- c('masked', 'driving', 'walking', 'transit', 'severity', 'Rt')
  start_date <- max(pred_tibble$date) + days(1)
  ydates <- seq(start_date,
                start_date + days(atomic_horizon - 1),
                by = 'd')
  new_rows <- as_tibble(y) %>%
    select(masked, severity, driving, walking, transit, Rt) %>%
    mutate(date = ydates,
           state = pred_tibble$state[1],
           popdensity = pred_tibble$popdensity[1],
           month = month(date)) %>%
    inner_join(psdata.b117prop, by = c('state', 'date')) %>%
    inner_join(data.weather.temp, by = c('state', 'month')) %>%
    mutate(temp = (temp - mean(covariates$temp)) / sd(covariates$temp)) %>%
    select(-month)
  
  rbind(pred_tibble, new_rows)
}

plot_predict_Rt <- function(st8, last_real, horizon) {
  prep <- prepare_pred_tibble(st8, last_real)
  for (i in 1:horizon) {
    prep <- predict_next_day(prep)
  }
  ggplot(mutate(prep, predicted = date > last_real)) +
    geom_line(aes(date, Rt, color=predicted))
}

predict_Rt_ts <- function(st8, last_real, horizon) {
  prep <- prepare_pred_tibble(st8, last_real)
  for (i in 1:horizon) {
    prep <- predict_next_day(prep)
  }
  prep_pred <- filter(prep, date > last_real)
  xts(prep_pred$Rt, prep_pred$date)
}

# Next step: Acquire Vt prediction from ARIMA

predict_Vt_ts <- function(st8, last_real, horizon) {
  first_vac_date <- (indicators.can %>% filter(state == st8, vacs > 0) %>% arrange(date))$date[1]
  vac_df <- indicators.can %>%
    filter(state == st8, date >= first_vac_date) %>%
    arrange(date) %>%
    fill(vacs) %>%
    filter(date <= last_real)
  
  vac_ts <- xts(vac_df$vacs, vac_df$date)
  vac_arima <- arima(vac_ts, order = c(4, 2, 2))
  # cic.discrete(vac_ts, min(vac_df$date), max(vac_df$date))
  vac_forecast <- pmin(forecast(vac_arima, h=horizon + 40)$mean, data.population.states[st8])
  dv <- c(vac_df$vacs[1], diff(c(vac_df$vacs, vac_forecast)))
  vac_fc_dates <- seq(last_real + days(1), last_real + days(horizon + 40), by = 'days')
  
  vac_full_ts <- xts(dv, c(vac_df$date, vac_fc_dates))
  cic_ts <- cic.discrete(vac_full_ts, min(index(vac_full_ts)), max(index(vac_full_ts)))
  
  pred_start <- as.Date(last_real + days(1))
  pred_end <- as.Date(last_real + days(horizon))
  diff(cic_ts)[seq(pred_start, pred_end, by = 'days')]
}

# Put it all together and pass it to run.simulation

clean_prime <- function(prime) {
  prime$params$ipc <- prime$params$ipc[!is.na(prime$params$ipc)]
  prime
}

simulate_from_Rt_Vt <- function(st8, horizon = 30) {
  prime <- fittings[[st8]][3][[1]]$prime %>% clean_prime
  last_real <- as.Date(prime$compartments$date)
  pred_end <- last_real + days(horizon)
  
  # Override compartments with truth
  real_row <- filter(indicators.can, state == st8, date == last_real)
  prime$compartments$Q <- real_row$cases.smooth
  prime$compartments$D <- real_row$deaths.smooth
  prime$compartments$H <- real_row$hosp.smooth
  
  # Assume ipc stays constant
  prime$params$ipc <- c(prime$params$ipc, rep(last(prime$params$ipc),
                                              max(0, horizon - length(prime$params$ipc))))
  
  Rt <- predict_Rt_ts(st8, last_real, horizon) %>% as.numeric
  Vt <- predict_Vt_ts(st8, last_real, horizon) %>% as.numeric
  
  sim <- run.simulation(prime$params, prime$queues, prime$compartments,
                 Rt, Vt, last_real + days(1))
  
  list(Rt=Rt, Vt=Vt, sim=sim, pred_end=pred_end, state=st8)
}

