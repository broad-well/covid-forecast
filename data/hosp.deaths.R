library(RSocrata)

data.hosp.raw <- read.socrata('https://healthdata.gov/resource/g62h-syeh.json?$select=state,date,inpatient_beds_used_covid,previous_day_admission_adult_covid_confirmed,previous_day_admission_adult_covid_suspected,previous_day_admission_pediatric_covid_confirmed,previous_day_admission_pediatric_covid_suspected') %>%
  mutate(hospitalized = as.numeric(inpatient_beds_used_covid),
         admission = as.numeric(previous_day_admission_adult_covid_suspected) +
                     as.numeric(previous_day_admission_adult_covid_confirmed) +
                     as.numeric(previous_day_admission_pediatric_covid_suspected) +
                     as.numeric(previous_day_admission_pediatric_covid_confirmed),
         date = ymd(date))
