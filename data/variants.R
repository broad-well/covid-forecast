library(tidyverse)

data.variants.raw <- read_csv('https://raw.githubusercontent.com/broad-well/covid-variants-us-states/main/us-states.csv')
data.variants.raw.b117 <- data.variants.raw %>%
  select(date, state, b117) %>%
  group_by(state) %>%
  mutate(cases=b117, daily=c(0, diff(b117) / as.numeric(diff(date), units = "days")))

# library(openxlsx)

# download.file('https://www.cdc.gov/coronavirus/2019-ncov/cases-updates/Data-from-variant-barg-graph-210311.csv', destfile='variantProportions.xlsx')
# data.variants.proportions <- readxl::read_xlsx('variantProportions.xlsx')
# data.variants.proportions.b117 <- data.variants.proportions %>%
#   filter(Lineage == 'B.1.1.7') %>%
#   mutate(date=as.Date(`Day of Collection Biweekly  Ending Date`), prop=`% of Total`) %>%
#   select(date, prop)

data.variants.b117 <- data.frame(date = c(16, 30, 44, 58, 72, 86), prop = c(0.021, 0.025, 0.044, 0.115, 0.274, 0.447))

prop.nls <- nls(prop ~ 1 / (1 + exp(k * date + X)), data = data.variants.b117, start = list(k=-0.2, X=10))

# ggplot(prop.test) + geom_point(aes(date, prop)) + geom_ribbon(aes(date, y=pred, ymax=c90, ymin=c10), alpha = .25) + geom_line(aes(date, pred), alpha = 0.5) + theme_bw(base_family='FreeSans') + labs(title = 'Proportion of B.1.1.7 in the U.S. over 2021')
