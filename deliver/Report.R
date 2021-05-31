library(gridExtra)
source('data/covariates.R')

report_state <- function(sim_output) {
  population <- data.population.states[[sim_output$state]]
  history <- sim_output$sim$history %>%
    mutate(Daily.Deaths = c(NA, diff(D)), Daily.Cases = c(NA, diff(Q)), Immune = R / population) %>%
    rename(Hospitalized = H) %>%
    select(date, Daily.Cases, Daily.Deaths, Hospitalized, Immune) %>%
    pivot_longer(c(Daily.Cases, Daily.Deaths, Hospitalized, Immune))
  
  ggplot(history) +
    geom_line(aes(date, value)) +
    facet_wrap(~ name, scales = 'free_y')
}

library(cowplot)
show_all_flowprobs <- function() {
  dists <- list(
    `S→E`=initial_params$SE,
    `E→Q`=initial_params$EQ,
    `E→H`=initial_params$EH,
    `E→R`=initial_params$ER,
    `H→R`=initial_params$HR,
    `Q→D`=initial_params$QD,
    `H→D`=initial_params$HD)
  order <- c('S→E', 'E→Q', 'E→H', 'E→R', 'H→R', 'H→D', 'Q→D')
  all_df <- dists %>%
    imap(~ tibble(day = 1:length(.x), p = .x, name = .y)) %>%
    reduce(rbind)
  
  cairo_pdf('deliver/flowprobs.pdf', width = 7, height = 4)
  plot(ggplot(all_df) + geom_line(aes(day, p)) +
    facet_wrap(~ factor(name, levels = order), scales = 'free') +
    labs(title = 'Probability Distributions of Inter-Compartmental Connections') +
    scale_y_continuous(labels = NULL) + xlab('Day') + ylab('') +
    global.theme)
  dev.off()
}

export_diagnoses <- function() {
  for (st8 in names(fittings.new)) {
    cairo_pdf(str_glue('deliver/fitting-{st8}.pdf'), width = 10, height = 2.5)
    plot(plot.diagnosis(fittings.new[[st8]]))
    dev.off()
  }
}

plot_rnn_ablation <- function(st8, start_date = ymd('2021-02-01')) {
  res.allOverridden <- demonstrate_rt(st8, start_date,
                                      selector = function(x) select(x, Rt))
  res.mobilityLeftOut <- demonstrate_rt(st8, start_date,
                                        selector = function(x) select(x, driving, walking, transit, Rt))
  res.mobilitySeverityLeftOut <- demonstrate_rt(
    st8, start_date, selector = function(x) select(x, driving, walking, transit, severity, Rt))
  
  res.allLeftOut <- demonstrate_rt(st8, start_date)
  
  res.all <- rbind(
    mutate(res.allOverridden, ablation = 'None'),
    mutate(res.mobilityLeftOut, ablation = 'Mobility'),
    mutate(res.mobilitySeverityLeftOut, ablation = 'Mobility, Severity'),
    mutate(res.allLeftOut, ablation = 'All')) %>%
    pivot_longer(c(P.masked:P.Rt, masked:Rt)) %>%
    mutate(predicted = startsWith(name, 'P.'),
           name = recode(name, P.masked = 'masked',
                         P.driving = 'driving',
                         P.walking = 'walking',
                         P.transit = 'transit',
                         P.severity = 'severity',
                         P.Rt = 'Rt')) %>%
    filter(!(name %in% c('b117prop', 'cases', 'temp')))
  
  the.plot <- ggplot(res.all) +
    geom_line(aes(date, value, color=predicted)) +
    facet_grid(rows = vars(ablation), cols = vars(name), scales = 'free_y')
  
  list(data = res.all, plot = the.plot)
}
