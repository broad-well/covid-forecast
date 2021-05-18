show_cases <- function(st) {
  indicators.can %>%
    filter(state %in% st, date > '2021-01-01') %>%
    group_by(state) %>% arrange(date) %>%
    ggplot(aes(date, c(NA, diff(cases)))) + geom_line() + facet_wrap(~state, scales = 'free_y')
}
