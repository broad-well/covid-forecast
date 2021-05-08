source('data/covariates.R')
source('data/cases.R')

library(lubridate)

fit_temp_cases <- function (statename, degree = 1) {
  monthly_cases <- data.cases.raw %>%
    filter(state == statename, date >= '2020-04-01') %>%
    group_by(month(date)) %>%
    summarize(cases = (max(cases) - min(cases)) / as.numeric(max(date) - min(date), unit = 'days')) %>%
    rename(month = `month(date)`)
  monthly_temp <- data.weather.temp %>%
    filter(Location == statename)
  cases.temp <- inner_join(monthly_cases, monthly_temp, by = 'month')
  list(cases.temp, lm(log2(cases) ~ poly(temp, degree), data = cases.temp))
}

plot.fit <- function (data, fit, date_metric='month') {
  pdata <- tibble(pred = fitted(fit), month = 1:12) %>% inner_join(data, by=date_metric)
  ggplot(pdata, aes(x=month)) + geom_line(aes(y=cases), color='red') +
    geom_line(aes(y=2^pred), color='blue')
}
