library(tidyverse)
library(zoo)

data.cases.raw <- read_csv('https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-states.csv')
data.cases.smoothed <- data.cases.raw %>%
  group_by(state) %>%
  summarize(cases.smooth = rollmean(cases, k = 7, na.pad = TRUE),
            deaths.smooth = rollmean(deaths, k = 7, na.pad = TRUE),
            date = date)

# DC[t] / DC[t-1]
# exp(log(DC[t]) - log(DC[t-1]))
data.cases.rt <- data.cases.smoothed %>%
  group_by(state) %>%
  summarize(daily = c(0, diff(cases.smooth)), date = date, rt = lead(daily) / daily)

data.cases.cdc <- read_csv('https://data.cdc.gov/api/views/9mfq-cb36/rows.csv?accessType=DOWNLOAD') %>%
  mutate(date = mdy(submission_date))
