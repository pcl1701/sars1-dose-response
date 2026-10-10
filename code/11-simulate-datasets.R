rm(list = ls())

## Install and load necessary packages -----------------------------------------

req_pkgs <- c("dplyr", "tidybayes")
new_pkgs <- req_pkgs[!(req_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) {
  install.packages(new_pkgs, dependences = TRUE)
}
invisible(lapply(req_pkgs, library, character.only = TRUE))
# load helper functions
source("code/helper-funcs.R")

## Load data and model fit -----------------------------------------------------

data <- read.csv("data/clean-data.csv")
data$individual <- as.factor(data$individual)
data$animal_species <- as.factor(data$animal_species)
data$study_name <- as.factor(data$study_name)
data$viral_isolate <- as.factor(data$viral_isolate)
data$data_type <- as.factor(data$data_type)

fit <- readRDS("outputs/fits/fit-all-data.RDS")
fit <- fit %>%
  recover_types(data)

## extract parameters ----------------------------------------------------------

set.seed(1234)

draws <- fit %>%
  spread_draws(hit_prob,
               alpha_t, log_tp5[study_name],
               alpha_N, log_N5[study_name],
               g[study_name], d[study_name],
               sigma_s[study_name], p_false)

## simulate original doses, uneven distribution --------------------------------

set.seed(1234)

true_log10ID50s <- c(2, 3, 4, 5, 6)
true_ID50s <- 10^true_log10ID50s
true_hit_probs <- log(2) / true_ID50s

for (j in 1:length(true_hit_probs)) {

  study_names <- c("BL", "JL", "NN", "YC", "LL16", "LL19")
  n_datasets <- 10

  sim_dataset <- setNames(
    data.frame(matrix(ncol = 11, nrow = 0)),
    c(
      "dataset_id",
      "individual",
      "animal_species",
      "study_name",
      "viral_isolate",
      "TCID_dose",
      "day_post_infection",
      "data_type",
      "log10_viral_load_og",
      "log10_viral_load",
      "log10_LOD"
    )
  )

  for (i in 1:n_datasets) {

    data_id <- i

    # sample parameter values
    n_draws <- length(unique(draws$.draw))
    sampled_draw <- sample(1:n_draws, size = 1)
    draws_filter <- draws %>% filter(.draw == sampled_draw)
    hit_prob <- true_hit_probs[j]
    alpha_t <- draws_filter$alpha_t[1]
    alpha_N <- draws_filter$alpha_N[1]
    p_false <- draws_filter$p_false[1]

    for (name in study_names) {

      # get study-specific parameters
      draws_study_filter <- draws_filter %>% filter(study_name == name)
      g <- draws_study_filter$g
      d <- draws_study_filter$d
      log_tp5 <- draws_study_filter$log_tp5
      log_N5 <- draws_study_filter$log_N5
      sigma_s <- draws_study_filter$sigma_s

      # filter data by study
      data_study_filter <- data %>% filter(study_name == name)

      # simulate data for each individual
      indiv_ids <- unique(data_study_filter$individual)
      for (indiv in indiv_ids) {

        # filter by individual
        data_indiv_filter <- data_study_filter %>% filter(individual == indiv)
        dose <- data_indiv_filter$TCID_dose[1]
        LOD <- data_indiv_filter$log10_LOD[1]

        # simulate infection, chance that monkey is not infected
        inf_prob <- DR_prob(hit_prob, dose)

        # simulate viral loads
        time_points <- data_indiv_filter$day_post_infection
        if (runif(1) > inf_prob) {
          # for failed infections
          log10_viral_loads_og <- rbinom(length(time_points), 1, prob=p_false)
          log10_viral_loads_og[log10_viral_loads_og == 0] <- -100 # true negative
          log10_viral_loads_og[log10_viral_loads_og == 1] <- 100 # false positive
          log10_viral_loads <- log10_viral_loads_og
        } else {
          # for successful infections
          # compute predicted trajectory
          t_p <- peak_time(hit_prob, dose, alpha_t, log_tp5)
          logN_0 <- log_init_RNA(hit_prob, dose, alpha_N, log_N5)
          # simulate data from predicted trajectory
          mean_viral_loads <- unlist(lapply(time_points,
                                            exp_growth_decay,
                                            g_i = g,
                                            d_i = d,
                                            t_p = t_p,
                                            logN_0 = logN_0))
          log_viral_loads_og <- rnorm(length(time_points),
                                   mean = mean_viral_loads,
                                   sd = sigma_s)
          log10_viral_loads_og <- log_to_log10(log_viral_loads_og)
          log10_viral_loads <- log10_viral_loads_og
          log10_viral_loads[log10_viral_loads < LOD] <- -100
        }

        # build individual simulate data
        indiv_sim_dataset <- data.frame(dataset_id = data_id,
                                        individual = indiv,
                                        animal_species = data_indiv_filter$animal_species[1],
                                        study_name = name,
                                        viral_isolate = data_indiv_filter$viral_isolate[1],
                                        TCID_dose = dose,
                                        day_post_infection = time_points,
                                        data_type = data_indiv_filter$data_type,
                                        log10_viral_load_og = log10_viral_loads_og,
                                        log10_viral_load = log10_viral_loads,
                                        log10_LOD = LOD)

        # combine with all simulate data
        sim_dataset <- rbind(sim_dataset, indiv_sim_dataset)
      }

    }

  }

  write.csv(
    sim_dataset,
    paste0(
      "outputs/simulations/data/simulated-data-original-doses-uneven-distribution-fixed-log10ID50-",
      true_log10ID50s[j],
      ".csv"
    ),
    row.names = FALSE
  )

}

## simulate original doses, balanced distribution ------------------------------

set.seed(1234)

true_log10ID50s <- c(2, 3, 4, 5, 6)
true_ID50s <- 10^true_log10ID50s
true_hit_probs <- log(2) / true_ID50s

study_names <- c("BL", "JL", "NN", "YC", "LL16", "LL19")
indiv_names <- unique(data$individual)
species_names <- unique(data$animal_species)
isolate_names <- unique(data$viral_isolate)

data_dose_scramble <- data

# 39 total individuals, across 4 doses
# balanced design:
# 10 at 10^3
# 10 at 10^5
# 10 at 10^6
# 9 at 10^6.9 (8571429)
indiv_doses <- c(rep(10^3, 10), rep(10^5, 10), rep(10^6, 10), rep(8571429, 9))

for (j in 1:length(true_hit_probs)) {

  study_names <- c("BL", "JL", "NN", "YC", "LL16", "LL19")
  n_datasets <- 10

  sim_dataset <- setNames(
    data.frame(matrix(ncol = 11, nrow = 0)),
    c(
      "dataset_id",
      "individual",
      "animal_species",
      "study_name",
      "viral_isolate",
      "TCID_dose",
      "day_post_infection",
      "data_type",
      "log10_viral_load_og",
      "log10_viral_load",
      "log10_LOD"
    )
  )

  for (i in 1:n_datasets) {

    # scramble the doses
    dose_assignment <- data.frame(
      individual = unique(data_dose_scramble$individual),
      TCID_dose = sample(indiv_doses)
    )
    data_dose_scramble <- data_dose_scramble %>%
      select(-TCID_dose) %>%
      left_join(dose_assignment, by = "individual")

    data_id <- i

    # sample parameter values
    n_draws <- length(unique(draws$.draw))
    sampled_draw <- sample(1:n_draws, size = 1)
    draws_filter <- draws %>% filter(.draw == sampled_draw)
    hit_prob <- true_hit_probs[j]
    alpha_t <- draws_filter$alpha_t[1]
    alpha_N <- draws_filter$alpha_N[1]
    p_false <- draws_filter$p_false[1]

    for (name in study_names) {

      # get study-specific parameters
      draws_study_filter <- draws_filter %>% filter(study_name == name)
      g <- draws_study_filter$g
      d <- draws_study_filter$d
      log_tp5 <- draws_study_filter$log_tp5
      log_N5 <- draws_study_filter$log_N5
      sigma_s <- draws_study_filter$sigma_s

      # filter data by study
      data_study_filter <- data_dose_scramble %>% filter(study_name == name)

      # simulate data for each individual
      indiv_ids <- unique(data_study_filter$individual)
      for (indiv in indiv_ids) {

        # filter by individual
        data_indiv_filter <- data_study_filter %>% filter(individual == indiv)
        dose <- data_indiv_filter$TCID_dose[1]
        LOD <- data_indiv_filter$log10_LOD[1]

        # simulate infection, chance that monkey is not infected
        inf_prob <- DR_prob(hit_prob, dose)

        # simulate viral loads
        time_points <- data_indiv_filter$day_post_infection
        if (runif(1) > inf_prob) {
          # for failed infections
          log10_viral_loads_og <- rbinom(length(time_points), 1, prob=p_false)
          log10_viral_loads_og[log10_viral_loads_og == 0] <- -100 # true negative
          log10_viral_loads_og[log10_viral_loads_og == 1] <- 100 # false positive
          log10_viral_loads <- log10_viral_loads_og
        } else {
          # for successful infections
          # compute predicted trajectory
          t_p <- peak_time(hit_prob, dose, alpha_t, log_tp5)
          logN_0 <- log_init_RNA(hit_prob, dose, alpha_N, log_N5)
          # simulate data from predicted trajectory
          mean_viral_loads <- unlist(lapply(time_points,
                                            exp_growth_decay,
                                            g_i = g,
                                            d_i = d,
                                            t_p = t_p,
                                            logN_0 = logN_0))
          log_viral_loads_og <- rnorm(length(time_points),
                                      mean = mean_viral_loads,
                                      sd = sigma_s)
          log10_viral_loads_og <- log_to_log10(log_viral_loads_og)
          log10_viral_loads <- log10_viral_loads_og
          log10_viral_loads[log10_viral_loads < LOD] <- -100
        }

        # build individual simulate data
        indiv_sim_dataset <- data.frame(dataset_id = data_id,
                                        individual = indiv,
                                        animal_species = data_indiv_filter$animal_species[1],
                                        study_name = name,
                                        viral_isolate = data_indiv_filter$viral_isolate[1],
                                        TCID_dose = dose,
                                        day_post_infection = time_points,
                                        data_type = data_indiv_filter$data_type,
                                        log10_viral_load_og = log10_viral_loads_og,
                                        log10_viral_load = log10_viral_loads,
                                        log10_LOD = LOD)

        # combine with all simulate data
        sim_dataset <- rbind(sim_dataset, indiv_sim_dataset)
      }

    }

  }

  write.csv(
    sim_dataset,
    paste0(
      "outputs/simulations/data/simulated-data-original-doses-balanced-distribution-fixed-log10ID50-",
      true_log10ID50s[j],
      ".csv"
    ),
    row.names = FALSE
  )

}

## simulate new doses, balanced distribution -----------------------------------

set.seed(1234)

true_log10ID50s <- c(2, 3, 4, 5, 6)
true_ID50s <- 10^true_log10ID50s
true_hit_probs <- log(2) / true_ID50s

study_names <- c("BL", "JL", "NN", "YC", "LL16", "LL19")
indiv_names <- unique(data$individual)
species_names <- unique(data$animal_species)
isolate_names <- unique(data$viral_isolate)

data_dose_scramble <- data

# 39 total individuals, doses: 10^2, 10^3, 10^4, 10^5, 10^6
# balanced design:
# 8 at 10^2
# 8 at 10^3
# 8 at 10^4
# 8 at 10^5
# 7 at 10^6
indiv_doses <- c(rep(10^2, 8), rep(10^3, 8), rep(10^4, 8), rep(10^5, 8), rep(10^6, 7))

for (j in 1:length(true_hit_probs)) {
  
  study_names <- c("BL", "JL", "NN", "YC", "LL16", "LL19")
  n_datasets <- 10
  
  sim_dataset <- setNames(
    data.frame(matrix(ncol = 11, nrow = 0)),
    c(
      "dataset_id",
      "individual",
      "animal_species",
      "study_name",
      "viral_isolate",
      "TCID_dose",
      "day_post_infection",
      "data_type",
      "log10_viral_load_og",
      "log10_viral_load",
      "log10_LOD"
    )
  )
  
  for (i in 1:n_datasets) {
    
    # scramble the doses
    dose_assignment <- data.frame(
      individual = unique(data_dose_scramble$individual),
      TCID_dose = sample(indiv_doses)
    )
    data_dose_scramble <- data_dose_scramble %>%
      select(-TCID_dose) %>%
      left_join(dose_assignment, by = "individual")
    
    data_id <- i
    
    # sample parameter values
    n_draws <- length(unique(draws$.draw))
    sampled_draw <- sample(1:n_draws, size = 1)
    draws_filter <- draws %>% filter(.draw == sampled_draw)
    hit_prob <- true_hit_probs[j]
    alpha_t <- draws_filter$alpha_t[1]
    alpha_N <- draws_filter$alpha_N[1]
    p_false <- draws_filter$p_false[1]
    
    for (name in study_names) {
      
      # get study-specific parameters
      draws_study_filter <- draws_filter %>% filter(study_name == name)
      g <- draws_study_filter$g
      d <- draws_study_filter$d
      log_tp5 <- draws_study_filter$log_tp5
      log_N5 <- draws_study_filter$log_N5
      sigma_s <- draws_study_filter$sigma_s
      
      # filter data by study
      data_study_filter <- data_dose_scramble %>% filter(study_name == name)
      
      # simulate data for each individual
      indiv_ids <- unique(data_study_filter$individual)
      for (indiv in indiv_ids) {
        
        # filter by individual
        data_indiv_filter <- data_study_filter %>% filter(individual == indiv)
        dose <- data_indiv_filter$TCID_dose[1]
        LOD <- data_indiv_filter$log10_LOD[1]
        
        # simulate infection, chance that monkey is not infected
        inf_prob <- DR_prob(hit_prob, dose)
        
        # simulate viral loads
        time_points <- data_indiv_filter$day_post_infection
        if (runif(1) > inf_prob) {
          # for failed infections
          log10_viral_loads_og <- rbinom(length(time_points), 1, prob=p_false)
          log10_viral_loads_og[log10_viral_loads_og == 0] <- -100 # true negative
          log10_viral_loads_og[log10_viral_loads_og == 1] <- 100 # false positive
          log10_viral_loads <- log10_viral_loads_og
        } else {
          # for successful infections
          # compute predicted trajectory
          t_p <- peak_time(hit_prob, dose, alpha_t, log_tp5)
          logN_0 <- log_init_RNA(hit_prob, dose, alpha_N, log_N5)
          # simulate data from predicted trajectory
          mean_viral_loads <- unlist(lapply(time_points,
                                            exp_growth_decay,
                                            g_i = g,
                                            d_i = d,
                                            t_p = t_p,
                                            logN_0 = logN_0))
          log_viral_loads_og <- rnorm(length(time_points),
                                      mean = mean_viral_loads,
                                      sd = sigma_s)
          log10_viral_loads_og <- log_to_log10(log_viral_loads_og)
          log10_viral_loads <- log10_viral_loads_og
          log10_viral_loads[log10_viral_loads < LOD] <- -100
        }
        
        # build individual simulate data
        indiv_sim_dataset <- data.frame(dataset_id = data_id,
                                        individual = indiv,
                                        animal_species = data_indiv_filter$animal_species[1],
                                        study_name = name,
                                        viral_isolate = data_indiv_filter$viral_isolate[1],
                                        TCID_dose = dose,
                                        day_post_infection = time_points,
                                        data_type = data_indiv_filter$data_type,
                                        log10_viral_load_og = log10_viral_loads_og,
                                        log10_viral_load = log10_viral_loads,
                                        log10_LOD = LOD)
        
        # combine with all simulate data
        sim_dataset <- rbind(sim_dataset, indiv_sim_dataset)
      }
      
    }
    
  }
  
  write.csv(
    sim_dataset,
    paste0(
      "outputs/simulations/data/simulated-data-new-doses-balanced-distribution-fixed-log10ID50-",
      true_log10ID50s[j],
      ".csv"
    ),
    row.names = FALSE
  )
  
}
