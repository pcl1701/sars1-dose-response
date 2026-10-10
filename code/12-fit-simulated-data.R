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

## Original doses, uneven distribution -----------------------------------------

# Notes on the data:
# -100 --> "-"
#  100 --> "+"

true_log10ID50s <- c(2, 3, 4, 5, 6)

for (j in 1:length(true_log10ID50s)) {

  df_all <- read.csv(paste0(
    "outputs/simulations/data/simulated-data-original-doses-uneven-distribution-fixed-log10ID50-",
    true_log10ID50s[j],
    ".csv"
  ))
  
  n_datasets <- 10
  
  for (i in 1:n_datasets) {
    
    print(paste("FITTING ORIGINAL DOSES, ORIGINAL DISTRIBUTION", i, "FIXED ID50", true_log10ID50s[j]))
    
    # filter by dataset id
    data_id <- i
    df <- df_all %>% filter(dataset_id == data_id)
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
    
    # load and fit model
    model <- cmdstan_model("code/stan-model-main.stan")
    fit <- model$sample(
      data = data_list,
      seed = 1234,
      chains = 4,
      parallel_chains = 4,
      refresh = 0,
      iter_sampling = 4000,
      adapt_delta = 0.99
    )
    
    # save fit
    data_id <- i
    fit$save_object(
      file = paste0(
        "outputs/simulations/fits/fit-simulated-data-original-doses-uneven-distribution-fixed-log10ID50-",
        true_log10ID50s[j],
        "-dataset-",
        data_id,
        ".RDS"
      )
    )
    
  }
}

## Original doses, balanced distribution ---------------------------------------

# Notes on the data:
# -100 --> "-"
#  100 --> "+"

true_log10ID50s <- c(2, 3, 4, 5, 6)

for (j in 1:length(true_log10ID50s)) {
  
  df_all <- read.csv(paste0(
    "outputs/simulations/data/simulated-data-original-doses-balanced-distribution-fixed-log10ID50-",
    true_log10ID50s[j],
    ".csv"
  ))
  
  n_datasets <- 10
  
  for (i in 1:n_datasets) {

    print(paste("FITTING ORIGINAL DOSES, BALANCED DISTRIBUTION", i, "FIXED ID50", true_log10ID50s[j]))
        
    # filter by dataset id
    data_id <- i
    df <- df_all %>% filter(dataset_id == data_id)
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
    
    # load and fit model
    model <- cmdstan_model("code/stan-model-main.stan")
    fit <- model$sample(
      data = data_list,
      seed = 1234,
      chains = 4,
      parallel_chains = 4,
      refresh = 0,
      iter_sampling = 4000,
      adapt_delta = 0.99
    )
    
    # save fit
    data_id <- i
    fit$save_object(
      file = paste0(
        "outputs/simulations/fits/fit-simulated-data-original-doses-balanced-distribution-fixed-log10ID50-",
        true_log10ID50s[j],
        "-dataset-",
        data_id,
        ".RDS"
      )
    )
    
  }
}

## New doses, balanced distribution ---------------------------------------

# Notes on the data:
# -100 --> "-"
#  100 --> "+"

true_log10ID50s <- c(2, 3, 4, 5, 6)

for (j in 1:length(true_log10ID50s)) {
  
  df_all <- read.csv(paste0(
    "outputs/simulations/data/simulated-data-new-doses-balanced-distribution-fixed-log10ID50-",
    true_log10ID50s[j],
    ".csv"
  ))
  
  n_datasets <- 10
  
  for (i in 1:n_datasets) {
    
    print(paste("FITTING NEW DOSES, BALANCED DISTRIBUTION", i, "FIXED ID50", true_log10ID50s[j]))
    
    # filter by dataset id
    data_id <- i
    df <- df_all %>% filter(dataset_id == data_id)
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
    
    # load and fit model
    model <- cmdstan_model("code/stan-model-main.stan")
    fit <- model$sample(
      data = data_list,
      seed = 1234,
      chains = 4,
      parallel_chains = 4,
      refresh = 0,
      iter_sampling = 4000,
      adapt_delta = 0.99
    )
    
    # save fit
    data_id <- i
    fit$save_object(
      file = paste0(
        "outputs/simulations/fits/fit-simulated-data-new-doses-balanced-distribution-fixed-log10ID50-",
        true_log10ID50s[j],
        "-dataset-",
        data_id,
        ".RDS"
      )
    )
    
  }
}

## Original doses, uneven distribution, LOD 1 ----------------------------------

# Notes on the data:
# -100 --> "-"
#  100 --> "+"

true_log10ID50s <- c(2, 3, 4, 5, 6)

for (j in 1:length(true_log10ID50s)) {
  
  df_all <- read.csv(paste0(
    "outputs/simulations/data/simulated-data-original-doses-uneven-distribution-fixed-log10ID50-",
    true_log10ID50s[j],
    ".csv"
  ))
  
  # use new LOD
  new_LOD <- 1
  df_all$log10_LOD <- new_LOD
  df_all$log10_viral_load <- df_all$log10_viral_load_og
  df_all$log10_viral_load[df_all$log10_viral_load < new_LOD] <- -100
  
  n_datasets <- 10
  
  for (i in 1:n_datasets) {
    
    print(paste("FITTING ORIGINAL DOSES, ORIGINAL DISTRIBUTION, LOD 1, ", i, "FIXED ID50", true_log10ID50s[j]))
    
    # filter by dataset id
    data_id <- i
    df <- df_all %>% filter(dataset_id == data_id)
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
    
    # load and fit model
    model <- cmdstan_model("code/stan-model-main.stan")
    fit <- model$sample(
      data = data_list,
      seed = 1234,
      chains = 4,
      parallel_chains = 4,
      refresh = 0,
      iter_sampling = 4000,
      adapt_delta = 0.99
    )
    
    # save fit
    data_id <- i
    fit$save_object(
      file = paste0(
        "outputs/simulations/fits/fit-simulated-data-original-doses-uneven-distribution-LOD1-fixed-log10ID50-",
        true_log10ID50s[j],
        "-dataset-",
        data_id,
        ".RDS"
      )
    )
    
  }
}

## Original doses, balanced distribution, LOD 1 --------------------------------

# Notes on the data:
# -100 --> "-"
#  100 --> "+"

true_log10ID50s <- c(2, 3, 4, 5, 6)

for (j in 1:length(true_log10ID50s)) {
  
  df_all <- read.csv(paste0(
    "outputs/simulations/data/simulated-data-original-doses-balanced-distribution-fixed-log10ID50-",
    true_log10ID50s[j],
    ".csv"
  ))

  # use new LOD
  new_LOD <- 1
  df_all$log10_LOD <- new_LOD
  df_all$log10_viral_load <- df_all$log10_viral_load_og
  df_all$log10_viral_load[df_all$log10_viral_load < new_LOD] <- -100
  
  n_datasets <- 10
  
  for (i in 1:n_datasets) {
    
    print(paste("FITTING ORIGINAL DOSES, BALANCED DISTRIBUTION, LOD 1, ", i, "FIXED ID50", true_log10ID50s[j]))
    
    # filter by dataset id
    data_id <- i
    df <- df_all %>% filter(dataset_id == data_id)
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
    
    # load and fit model
    model <- cmdstan_model("code/stan-model-main.stan")
    fit <- model$sample(
      data = data_list,
      seed = 1234,
      chains = 4,
      parallel_chains = 4,
      refresh = 0,
      iter_sampling = 4000,
      adapt_delta = 0.99
    )
    
    # save fit
    data_id <- i
    fit$save_object(
      file = paste0(
        "outputs/simulations/fits/fit-simulated-data-original-doses-balanced-distribution-LOD1-fixed-log10ID50-",
        true_log10ID50s[j],
        "-dataset-",
        data_id,
        ".RDS"
      )
    )
    
  }
}

## New doses, balanced distribution, LOD 1 -------------------------------------

# Notes on the data:
# -100 --> "-"
#  100 --> "+"

true_log10ID50s <- c(2, 3, 4, 5, 6)

for (j in 1:length(true_log10ID50s)) {
  
  df_all <- read.csv(paste0(
    "outputs/simulations/data/simulated-data-new-doses-balanced-distribution-fixed-log10ID50-",
    true_log10ID50s[j],
    ".csv"
  ))
  
  # use new LOD
  new_LOD <- 1
  df_all$log10_LOD <- new_LOD
  df_all$log10_viral_load <- df_all$log10_viral_load_og
  df_all$log10_viral_load[df_all$log10_viral_load < new_LOD] <- -100
  
  n_datasets <- 10
  
  for (i in 1:n_datasets) {
    
    print(paste("FITTING NEW DOSES, BALANCED DISTRIBUTION, LOD 1, ", i, "FIXED ID50", true_log10ID50s[j]))
    
    # filter by dataset id
    data_id <- i
    df <- df_all %>% filter(dataset_id == data_id)
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
    
    # load and fit model
    model <- cmdstan_model("code/stan-model-main.stan")
    fit <- model$sample(
      data = data_list,
      seed = 1234,
      chains = 4,
      parallel_chains = 4,
      refresh = 0,
      iter_sampling = 4000,
      adapt_delta = 0.99
    )
    
    # save fit
    data_id <- i
    fit$save_object(
      file = paste0(
        "outputs/simulations/fits/fit-simulated-data-new-doses-balanced-distribution-LOD1-fixed-log10ID50-",
        true_log10ID50s[j],
        "-dataset-",
        data_id,
        ".RDS"
      )
    )
    
  }
}

