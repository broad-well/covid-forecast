library(tidyverse)

data.weather.temp.raw <- read_csv('https://www.ncdc.noaa.gov/cag/statewide/mapping/110-tavg.csv', skip=3)
data.weather.temp <- data.weather.temp.raw %>%
  mutate(month = Date %% 100) %>%
  group_by(month, Location) %>%
  summarize(temp = mean(`1901-2000 Mean`))

data.ages <- read_csv('https://www2.census.gov/programs-surveys/popest/tables/2010-2019/state/asrh/sc-est2019-agesex-civ.csv') %>%
  select(STATE, NAME, POPEST2019_CIV, AGE) %>%
  filter(STATE > 0, AGE < 200) %>%
  mutate(age.group = cut(AGE, breaks = c(-1, 17, 49, 64, 200))) %>%
  group_by(NAME, age.group) %>%
  summarize(size = sum(POPEST2019_CIV)) %>%
  mutate(state = if (NAME == 'District of Columbia') 'DC' else state.abb[match(NAME, state.name)]) %>%
  pivot_wider(names_from=age.group, values_from=size)

cc.masked <- covidcast_signal('fb-survey', 'smoothed_wothers_masked', '2020-11-24', '2021-05-10', 'state')