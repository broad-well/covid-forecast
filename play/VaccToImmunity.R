source('data/vaccinations.R')
library(tidycensus)
library(xts)
library(lubridate)
library(forecast)

# Individual efficacy curve (t=0 is fully vaccinated)
iec <- function(t) {
  coefs <- coef(eff.curve)
  coefs['feff'] / (1 + exp(coefs['r'] * (t + 25) + coefs['n']))
}

ts.test <- c(60, 5, 10, 30, 20)
# 60iec(t) + 5iec(t-1) + 10iec(t-2) + 30iec(t-3)
# Then t=0 corresponds to first date

# Cumulative efficacy curve
# ts: Time series of daily vaccinations
cec <- function(ts) function(t) {
  sum(ts * sapply(t - (0:(length(ts) - 1)), iec))
}

# More optimized (?) cumulative immunity curve
# Measured on 03-31: estimated 15.3x speedup (CI 14.08x-16.51x) compared to cic.continuous
cic.discrete <- function(ts, date.min, date.max) {
  output <- numeric((date.max - date.min) / ddays() + 1)
  offset <- (min(index(ts)) - date.min) / ddays()
  for (i in 1:length(ts)) {
    daily <- as.numeric(ts[i])
    output <- output + daily * iec((1:length(output)) - i - offset)
  }
  output.dates <- seq(date.min, date.max, by = 'days')
  xts(output, output.dates)
}

cic.continuous <- function(ts, date.min, date.max) {
  prefix <- (date.min - min(index(ts))) / ddays()
  domain <- prefix:(prefix + (date.max - date.min) / ddays())
  xts(sapply(domain, cec(ts)), seq(date.min, date.max, by = 'days'))
}

state.eff.curve <- function (state, prefix = -28,
                             gen.ts = cic.discrete, date.min = NA, date.max = NA) {
  fulvac.df <- data.vac.owid %>%
    filter(location == state) %>%
    select(date, people_fully_vaccinated) %>%
    fill(people_fully_vaccinated) %>%
    mutate(daily.vac = c(NA, diff(people_fully_vaccinated))) %>%
    drop_na()
  
  fulvac.ts <- xts(fulvac.df$daily.vac, fulvac.df$date)
  
  as.Date(ifelse(is.na(date.min), min(fulvac.df$date) - days(25), date.min)) -> date.min
  as.Date(ifelse(is.na(date.max), max(fulvac.df$date), date.max)) -> date.max
  
  cic.ts <- gen.ts(fulvac.ts, date.min, date.max)
  cec.df <- tibble(date = index(cic.ts),
                   ci = as.numeric(cic.ts)) %>% left_join(fulvac.df, by = 'date')
  cec.df
}

plot.eff.states <- function (states) {
  key <- read_file('auth/census.key')
  pop.est.df <- get_estimates('state', 'population', state = states, key = key) %>%
    filter(variable == 'POP')
  pop.est <- setNames(pop.est.df$value, pop.est.df$NAME)
  
  net.table <- tibble()
  for (state in states) {
    df <- state.eff.curve(state)
    df$state <- state
    net.table <- net.table %>% rbind(df)
  }
  ggplot(net.table, aes(x=date)) +
    geom_line(aes(y=people_fully_vaccinated / pop.est[state]), size=2) +
    geom_line(aes(y=ci / pop.est[state]), color = 'gray40', size=2) +
    facet_wrap(vars(state))
}

plot.owid.pred <- function(states, order = c(7, 2, 6)) {
  vac.df <- data.vac.owid %>%
    filter(location %in% states) %>%
    select(date, location, people_fully_vaccinated)

  make.state.xts <- function(.state) {
    vac.df.state <- vac.df %>% filter(location == .state) %>%
      mutate(daily = rollmean(c(0, diff(people_fully_vaccinated)), k = 7, na.pad = TRUE)) %>%
      drop_na(daily)
    xts(vac.df.state$daily, vac.df.state$date)
  }
  
  plot.state.pred <- function(.xts, name) {
    pred <- arima(.xts, order = c(7, 2, 6))
    print(pred)
    autoplot(forecast(pred, h = 21))
  }
  
  # sapply(states, function (s) (make.state.xts(s), s))
  plots <- list()
  
  for (.state in states) {
    plots[[.state]] <- make.state.xts(.state)
  }
  plots
}
