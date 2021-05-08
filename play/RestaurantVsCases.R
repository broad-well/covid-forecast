library(covidcast)
library(datasets)
library(stringr)

sig.rest <- covidcast_signal(data_source = 'safegraph', signal = 'restaurants_visit_prop',
                             start_day = '2020-05-01', end_day = '2021-03-20', geo_type = 'state')

fit_rest_cases <- function (stateq) {
  postal <- state.abb[match(stateq, state.name)]
  postal.lower <- str_to_lower(postal)
  rest.df <- sig.rest %>% filter(geo_value == postal.lower) %>%
    rename(restaurant=value, date=time_value)
  cases.df <- data.cases.rt %>% filter(`state` == stateq)
  combined.df <- inner_join(cases.df, rest.df, by = 'date') %>%
    mutate(daily.smooth = rollmean(daily, k=7, na.pad = TRUE),
           restaurant.smooth = rollmean(restaurant, k=7, na.pad = TRUE)) %>%
    drop_na(restaurant.smooth, daily.smooth) %>%
    dplyr::select(date, restaurant, restaurant.smooth, daily, daily.smooth, rt)
  print(combined.df)
  model <- lm(rt ~ poly(restaurant.smooth, 3) + log(restaurant.smooth), data = combined.df)
  list(data = combined.df, model = model)
}
