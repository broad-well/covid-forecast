library(tidyverse)

data.vac.raw <- read_csv('https://github.com/govex/COVID-19/raw/master/data_tables/vaccine_data/us_data/time_series/vaccine_data_us_timeline.csv')
data.vac.owid <- read_csv('https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/vaccinations/us_state_vaccinations.csv')

# Effectiveness over days from https://www.nejm.org/doi/full/10.1056/NEJMoa2101765
data.vac.eff.pfizer = tibble(
  point = c(0, .46, .60, .92, .92),
  lower = c(0, .40, .53, .88, .88),
  upper = c(0, .51, .66, .95, .95),
  week = c(0, 2.5, 3.5, 5.5, 7)
)

# Alternative from https://www.fda.gov/media/144246/download
eff.pfizer <- tibble(
  day = c(0, 7, 14, 21, 28, 35, 42, 49, 56, 63, 70, 77),
  pfizer.risk = c(21314, 21230, 21054, 20481, 19314, 18377, 17702, 17186, 15464, 14038, 12169, 9591),
  pfizer.inf = c(0, 21, 37, 39, 41, 42, 42, 43, 44, 47, 48, 48),
  placebo.risk = c(21258, 21170, 20970, 20366, 19209, 18218, 17578, 17025, 15290, 13876, 11994, 9471),
  placebo.inf = c(0, 25, 55, 73, 97, 123, 143, 166, 192, 212, 235, 249)) %>%
  mutate(eff = 1 - (c(0, diff(pfizer.inf)) / pfizer.risk) /
           (c(0, diff(placebo.inf)) / placebo.risk)) %>%
  filter(day < 60)

# eff.pfizer.pred <- glm(eff ~ day, data = eff.pfizer, family = 'binomial')

# Moderna data from https://www.fda.gov/media/144434/download (limited, Table 15)
eff.moderna <- tibble(
  day = c(0, 14, 21, 28 + 14),
  eff = c(0, 0.508, 0.921, 0.941)
)

# what if...
eff.curve <- nls(eff ~ feff / (1 + exp(r*day + n)),
                 data = rbind(eff.pfizer %>% select(day, eff) %>% filter(day < 60), eff.moderna),
                 start = list(r = -0.117, n = 1.3, feff = 0.92),
                 upper = c(Inf, Inf, 0.92),
                 algorithm = 'port',
                 trace = TRUE,
                 control = nls.control(maxiter = 40000, printEval = TRUE))

data.vac.ages.raw <- read_csv('data/raw/vac.dem.trends.csv') %>%
  rename(vac = `People who are fully vaccinated`, age.group = `Age Group`)
