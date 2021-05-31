library(tidyverse)

data.weather.temp.raw <- read_csv('https://www.ncdc.noaa.gov/cag/statewide/mapping/110-tavg.csv', skip=3)
data.weather.temp <- data.weather.temp.raw %>%
  mutate(month = Date %% 100) %>%
  group_by(month, Location) %>%
  summarize(temp = mean(`1901-2000 Mean`)) %>%
  mutate(Location = state.abb[match(Location, state.name)]) %>%
  rename(state = Location)

data.ages <- read_csv('https://www2.census.gov/programs-surveys/popest/tables/2010-2019/state/asrh/sc-est2019-agesex-civ.csv') %>%
  select(STATE, NAME, POPEST2019_CIV, AGE) %>%
  filter(STATE > 0, AGE < 200) %>%
  mutate(age.group = cut(AGE, breaks = c(-1, 17, 49, 64, 200))) %>%
  group_by(NAME, age.group) %>%
  summarize(size = sum(POPEST2019_CIV)) %>%
  mutate(state = if (NAME == 'District of Columbia') 'DC' else state.abb[match(NAME, state.name)]) %>%
  pivot_wider(names_from=age.group, values_from=size)

cc.masked <- covidcast_signal('fb-survey', 'smoothed_wothers_masked', '2020-11-24', '2021-05-10', 'state') %>%
  as_tibble %>% rename(state = geo_value, date = time_value) %>% arrange(date) %>%
  group_by(state) %>% mutate(masked = rollmean(value, k = 7, na.pad = TRUE) / 100, state = toupper(state))

data.mobility.apple <- read_csv('data/raw/apple-mobility.csv') %>%
  filter(region %in% state.name) %>%
  group_by(region) %>%
  pivot_longer(`2020-06-01`:`2021-05-14`, names_to = 'date', values_to = 'mobility') %>%
  select(region, date, transportation_type, mobility) %>%
  pivot_wider(names_from = 'transportation_type', values_from = 'mobility') %>%
  rename(state = region) %>%
  mutate(driving = na.approx(driving),
         walking = na.approx(walking),
         transit = if (length(unique(transit)) == 1) NA else na.approx(transit)) %>%
  mutate(driving = rollmean(driving, k=10, na.pad=TRUE) / 100 - 1,
         walking = rollmean(walking, k=10, na.pad=TRUE) / 100 - 1,
         transit = rollmean(transit, k=10, na.pad=TRUE) / 100 - 1,
         state = state.abb[match(state, state.name)],
         date = ymd(date))


library(covidcast)
data.population.density.raw <- read_csv('https://www2.census.gov/programs-surveys/decennial/2020/data/apportionment/apportionment.csv')
data.population.density <- data.population.density.raw %>%
  filter(`Geography Type` == 'State', Year == 2020) %>% 
  select(Name, `Resident Population Density`) %>%
  mutate(state = name_to_abbr(Name) %>% unname) %>%
  rename(density = `Resident Population Density`)
data.population.count <- data.population.density.raw %>%
  filter(`Geography Type` == 'State', Year == 2020) %>%
  select(Name, `Resident Population`) %>%
  rename(state = Name, pop = `Resident Population`) %>%
  mutate(state = state.abb[match(state, state.name)])
data.population.states <- data.population.count$pop
names(data.population.states) <- data.population.count$state
data.population.density.states <- data.population.density$density
names(data.population.density.states) <- data.population.density$state
