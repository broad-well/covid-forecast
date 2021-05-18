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

data.variants.b117 <- tibble(
  date = c(30, 44, 58, 72, 86, 100, 114),
  prop = c(0.012, 0.044, 0.115, 0.274, 0.446, 0.595, 0.66))

prop.nls <- nls(prop ~ mx / (1 + exp(k * date + X)), data = data.variants.b117, start = list(k=-0.2, X=10, mx=0.8))

data.variants.b117 %>% mutate(date = ymd('2020-12-31') + days(date)) -> data.variants.b117
# ggplot(prop.test) + geom_point(aes(date, prop)) + geom_ribbon(aes(date, y=pred, ymax=c90, ymin=c10), alpha = .25) + geom_line(aes(date, pred), alpha = 0.5) + theme_bw(base_family='FreeSans') + labs(title = 'Proportion of B.1.1.7 in the U.S. over 2021')

data.variants.b117.dates <- c(30, 44, 58, 72, 86, 100, 114)
data.variants.b117.regions <- list(
  c(0.004, 0.034, 0.115, 0.217, 0.370, 0.443, 0.509),
  c(0.029, 0.099, 0.146, 0.285, 0.371, 0.482, 0.545),
  c(0.009, 0.021, 0.098, 0.250, 0.462, 0.621, 0.665),
  c(0.018, 0.066, 0.150, 0.365, 0.545, 0.658, 0.705),
  c(0.006, 0.031, 0.081, 0.265, 0.499, 0.676, 0.731),
  c(0.013, 0.042, 0.143, 0.330, 0.505, 0.693, 0.741),
  c(0.000, 0.016, 0.077, 0.136, 0.447, 0.638, 0.724),
  c(0.000, 0.029, 0.047, 0.187, 0.354, 0.536, 0.568),
  c(0.003, 0.018, 0.048, 0.131, 0.251, 0.434, 0.577),
  c(0.019, 0.006, 0.079, 0.104, 0.258, 0.366, 0.523))
data.variants.b117.region.states <- list(
  c('CT', 'ME', 'MA', 'NH', 'RI', 'VT'),
  c('NJ', 'NY'),
  c('DE', 'DC', 'MD', 'PA', 'VA', 'WV'),
  c('AL', 'FL', 'GA', 'KY', 'MS', 'NC', 'SC', 'TN'),
  c('IL', 'IN', 'MI', 'MN', 'OH', 'WI'),
  c('AR', 'LA', 'NM', 'OK', 'TX'),
  c('IA', 'KS', 'MO', 'NE'),
  c('CO', 'MT', 'ND', 'SD', 'UT', 'WY'),
  c('AZ', 'CA', 'HI', 'NV'),
  c('AK', 'ID', 'OR', 'WA')
)

fit_growth <- function(region) {
  tib <- tibble(date = data.variants.b117.dates, prop = data.variants.b117.regions[[region]])
  prop.nls <- nls(prop ~ mx / (1 + exp(k * date + X)), data = tib, start = list(k=-0.2, X=15, mx=0.7))
  delta_days <- (-90):180
  predictor <- function(d) coef(prop.nls)['mx'] / (1 + exp(coef(prop.nls)['k'] * d + coef(prop.nls)['X']))
  predicted <- sapply(delta_days, predictor)
  tibble(date = ymd('2020-12-31') + days(delta_days), prop = predicted)
}

plot_growth <- function(region) {
  ggplot(fit_growth(region)) +
    geom_line(aes(date, prop)) +
    geom_point(data = tibble(date = ymd('2020-12-31') + days(data.variants.b117.dates),
                             prop = data.variants.b117.regions[[region]]), aes(date, prop)) +
    theme_bw(base_family = 'Linux Libertine O', base_size = 15) +
    ylab('Proportion B.1.1.7') +
    ggtitle(do.call(paste, c(as.list(data.variants.b117.region.states[[region]]), sep = ', ')))
}

library(gridExtra)
plot_growths <- function() {
  do.call(grid.arrange, c(map(1:9, plot_growth), nrow=3, ncol=3))
}
