rm(list = ls())

## Install and load necessary packages -----------------------------------------

req_pkgs <- c("cmdstanr", "dplyr", "tidybayes")
new_pkgs <- req_pkgs[!(req_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) {
  install.packages(new_pkgs, dependences = TRUE)
}
invisible(lapply(req_pkgs, library, character.only = TRUE))
# load helper functions
source("code/helper-funcs.R")

## Prepare data (cynomolgus data) ----------------------------------------------

# Notes on the data:
# -100 --> "-"
#  100 --> "+"

df <- read.csv("data/clean-data.csv")
df <- df %>% filter(animal_species == "CM")
df$individual <- as.factor(df$individual)
df$animal_species <- as.factor(df$animal_species)
df$study_name <- as.factor(df$study_name)
df$viral_isolate <- as.factor(df$viral_isolate)
df$data_type <- as.factor(df$data_type)

# convert log10 viral loads to log viral loads
data <- df %>%
  mutate(log_viral_load = case_when(log10_viral_load != -100 &
                                      log10_viral_load != 100 ~ log10_to_log(log10_viral_load),
                                    TRUE ~ log10_viral_load)) %>%
  mutate(log_LOD = log10_to_log(log10_LOD))

# index monkeys and studies
data <- compose_data(data)

# create vector of individual inoculation doses
# TCID_dose[i] is the dose of individual i
dose_df <- data.frame(individual_id = data$individual,
                      dose = data$TCID_dose)
TCID_dose <- arrange(unique(dose_df), individual_id)$dose

data_list <- list(
  n_monkeys = data$n_individual,
  n_datapoints = data$n,
  n_studies = data$n_study_name,
  individual_ids = data$individual,
  study_ids = data$study_name,
  TCID_dose = TCID_dose,
  day_post_infection = data$day_post_infection,
  log_viral_load = data$log_viral_load,
  logLOD = data$log_LOD
)

## Load and fit model ----------------------------------------------------------

model <- cmdstan_model("code/stan-model-main.stan")

fit <- model$sample(
  data = data_list,
  seed = 1234,
  chains = 4,
  parallel_chains = 4,
  refresh = 500,
  iter_sampling = 4000,
  adapt_delta = 0.99
)

## Save fit --------------------------------------------------------------------

fit$save_object(file = "outputs/fits/fit-cynomolgus-data.RDS")
