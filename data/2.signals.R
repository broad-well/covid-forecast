source('cases.R')
source('hosp.deaths.R')

indicators <- data.cases.cdc %>%
  mutate(date = mdy(submission_date)) %>%
  inner_join(data.hosp.raw) %>%
  group_by(state) %>%
  arrange(date) %>%
  summarize(dc = rollmean(new_case, k=7, na.pad=TRUE),
            dd = rollmean(new_death, k=7, na.pad=TRUE),
            hosp = rollmean(hospitalized, k=7, na.pad=TRUE),
            dh = rollmean(admission, k=7, na.pad=TRUE),
            date = date) %>%
  filter(date >= '2020-07-01')

indicators.can <- read_csv(paste(
  'https://api.covidactnow.org/v2/states.timeseries.csv?apiKey=',
  read_file('auth/covidactnow.key'), sep = '')) %>%
  select(state, date, fips, actuals.cases:actuals.negativeTests,
         actuals.hospitalBeds.currentUsageCovid,
         actuals.newCases,
         actuals.vaccinationsInitiated,
         actuals.vaccinationsCompleted,
         metrics.testPositivityRatio,
         actuals.newDeaths) %>%
  rename(cases = actuals.cases, deaths = actuals.deaths,
         test.pos = actuals.positiveTests, test.neg = actuals.negativeTests,
         hosp = actuals.hospitalBeds.currentUsageCovid,
         dc = actuals.newCases,
         vacs = actuals.vaccinationsCompleted,
         test.pos.ratio = metrics.testPositivityRatio,
         dd = actuals.newDeaths) %>%
  group_by(state) %>%
  arrange(date) %>%
  mutate(cases.smooth = rollmean(cases, k=7, na.pad=TRUE),
         hosp.smooth = rollmean(hosp, k=7, na.pad=TRUE),
         deaths.smooth = rollmean(deaths, k=7, na.pad=TRUE))
  filter(date >= '2020-06-01')
