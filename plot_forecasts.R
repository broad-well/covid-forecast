library(covidHubUtils)
library(dplyr)

all_models <- c(
  "COVIDhub-ensemble",
  "Caltech-CS156",
  "CovidAnalytics-DELPHI",
  "GT-DeepCOVID",
  "IHME-CurveFit",
  "LANL-GrowthRate",
  "UCLA-SuEIR",
  "UMass-MechBayes")
#  "OneQuietNight-ML")

locations <- c("06", "25", "46", "12")
checkpoints <- seq.Date(as.Date("2020-10-31"), as.Date("2021-01-23"), by="28 days")

load_checkpoint <- function(checkpoint) {
  return(load_latest_forecasts(models = all_models,
                               last_forecast_date = checkpoint,
                               forecast_date_window_size = 7,
                               source = "zoltar",
                               locations = locations,
                               types = c("quantile", "point"),
                               targets = paste(1:4, "wk ahead inc case")))
}

fdat <- bind_rows(lapply(checkpoints, load_checkpoint))
fdat <- load_latest_forecasts(models = all_models,
                      # forecast_dates = seq.Date(as.Date("2020-11-28"), as.Date("2021-01-23"), by="28 days"),
                       last_forecast_date = "2020-12-14",
                      forecast_date_window_size = 7, 
                      source = "zoltar", #
                       locations = locations,
                       types = c("quantile", "point"),
                       targets = paste(1:4, "wk ahead inc case"))


truth <- load_truth(truth_source = "NYTimes", target_variable = "inc case",
                    locations = locations) %>% filter(target_end_date > "2020-07-01", target_end_date < "2021-02-01")

library(ggplot2)
plot_forecast(fdat, locations = "46",
              target_variable = "inc case",
              truth_data = truth,
              truth_source = "NYTimes",
              intervals = c(.5, .95),
              facet = location~model, fill_by_model = TRUE, plot = TRUE) + ylim(0, 10000)


test_fdat <- load_forecasts(models = c( "Caltech-CS156"),
                            forecast_dates = seq.Date(as.Date("2020-11-01"), as.Date("2021-01-10"), by="14 days"),
                            locations = locations,
                            types = c("quantile", "point"),
                            targets = paste(1:3, "wk ahead inc case"))
