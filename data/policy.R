library(tidyverse)
library(lubridate)

data.policy.ox.raw <- read_csv('https://github.com/OxCGRT/covid-policy-tracker/raw/master/data/OxCGRT_latest.csv',
                               col_types = cols(RegionName = 'c', RegionCode = 'c')) %>%
  select(CountryCode, RegionName, RegionCode, Jurisdiction, Date, `H7_Vaccination policy`, H7_Flag) %>%
  filter(CountryCode == 'USA' & Jurisdiction == 'STATE_TOTAL') %>%
  mutate(Date = ymd(Date))

