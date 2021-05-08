library(tidyverse)
library(RSocrata)
library(lubridate)

data.surv.raw <- read.socrata('https://data.cdc.gov/resource/vbim-akqf.json?cdc_case_earliest_dt=2021-03-16')