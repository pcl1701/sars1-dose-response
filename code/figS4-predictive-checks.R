rm(list = ls())

## Install and load necessary packages -----------------------------------------

req_pkgs <- c("dplyr", "tidybayes", "ggplot2", "scales", "patchwork")
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

prior_fit <- readRDS("outputs/fits/fit-prior.RDS")
prior_fit <- prior_fit %>%
  recover_types(data)

post_fit <- readRDS("outputs/fits/fit-all-data.RDS")
post_fit <- post_fit %>%
  recover_types(data)

## Fig S4A: prior viral load trajectories --------------------------------------

set.seed(1234)

theme_set(theme_bw())

# extract parameters
draws <- prior_fit %>%
  spread_draws(hit_prob,
               alpha_t, log_tp5[study_name],
               alpha_N, log_N5[study_name],
               g[study_name], d[study_name])

# sample a subset (100) of extracted parameters
n_draws <- length(unique(draws$.draw))
draws_per_study <- 100
sampled_draws <- sample(1:n_draws, size = draws_per_study)
draws <- draws %>%
  filter(.draw %in% sampled_draws)

# keep track of doses (TCID50) used in each study
study_doses <- list("BL" = c(10^5),
                    "LL16" = c(10^5),
                    "LL19" = c(10^5),
                    "NN" = c(10^3, 10^6),
                    "JL" = c(6*10^6/0.7),
                    "YC" = c(10^5))

# vectors to store successful trajectories
final_times <- c()
final_viral_loads <- c()
final_study_names <- c()
final_doses <- c()
final_ids <- c() # to keep track of individual rows (draws)

# vectors to store failed infections
failed_study_names <- c()
failed_doses <- c()

# times to compute trajectories
n_time_points <- 100
t_seq <- seq(from = 0, to = 35, length.out = n_time_points)

# compute and save trajectories for each draw
for (i in 1:nrow(draws)){
  
  # extract parameters from draw
  g <- draws[i,]$g
  d <- draws[i,]$d
  hit_prob <- draws[i,]$hit_prob
  alpha_t <- draws[i,]$alpha_t
  alpha_N <- draws[i,]$alpha_N
  log_tp5 <- draws[i,]$log_tp5
  log_N5 <- draws[i,]$log_N5
  study_name <- as.character(draws[i,]$study_name)
  
  for (dose in study_doses[[study_name]]){
    
    # simulate infection, chance that monkey is not infected
    inf_prob <- DR_prob(hit_prob, dose)
    # keep track of failed infections
    if (runif(1) > inf_prob) {
      failed_study_names <- c(failed_study_names, study_name)
      failed_doses <- c(failed_doses, dose)
      next
    }
    
    # if infected, compute predicted trajectory
    t_p <- peak_time(hit_prob, dose, alpha_t, log_tp5)
    logN_0 <- log_init_RNA(hit_prob, dose, alpha_N, log_N5)
    viral_loads <- unlist(lapply(t_seq,
                                 exp_growth_decay,
                                 g_i = g,
                                 d_i = d,
                                 t_p = t_p,
                                 logN_0 = logN_0))
    viral_loads <- exp(viral_loads)
    
    # save trajectory
    final_times <- c(final_times, t_seq)
    final_viral_loads <- c(final_viral_loads, viral_loads)
    final_study_names <- c(final_study_names, rep(study_name, n_time_points))
    final_doses <- c(final_doses, rep(dose, n_time_points))
    final_ids <- c(final_ids, rep(i, n_time_points))
  }
  
}

# data frame for plotting successful trajectories
df_inf <- data.frame(time = final_times,
                     viral_load = final_viral_loads,
                     study_name = final_study_names,
                     dose = final_doses,
                     id = final_ids) %>%
  filter(viral_load >= 1) %>% # remove points below 10^0
  mutate(facet = case_when(study_name == "JL" ~ "'JL -' ~ 10^6.9 ~ TCID[50]",
                           study_name == "BL" ~ "'BL -' ~ 10^5 ~ TCID[50]",
                           study_name == "LL16" ~ "'LL16 -' ~ 10^5 ~ TCID[50]",
                           study_name == "LL19" ~ "'LL19 -' ~ 10^5 ~ TCID[50]",
                           study_name == "YC" ~ "'YC -' ~ 10^5 ~ TCID[50]",
                           study_name == "NN" & dose == 1e+3 ~ "'NN -' ~ 10^3 ~ TCID[50]",
                           study_name == "NN" & dose == 1e+6 ~ "'NN -' ~ 10^6 ~ TCID[50]"))
# order facet variable for plotting
df_inf$facet <- factor(df_inf$facet, levels = c("'NN -' ~ 10^3 ~ TCID[50]",
                                                "'BL -' ~ 10^5 ~ TCID[50]",
                                                "'YC -' ~ 10^5 ~ TCID[50]",
                                                "'LL16 -' ~ 10^5 ~ TCID[50]",
                                                "'LL19 -' ~ 10^5 ~ TCID[50]",
                                                "'NN -' ~ 10^6 ~ TCID[50]",
                                                "'JL -' ~ 10^6.9 ~ TCID[50]"))

# create trajectories for failed infections
# lines that will be plotted below zero to indicate failed infections
uninf_study_name <- c()
uninf_dose <- c()
uninf_time <- c()
uninf_viral_load <- c()
uninf_id <- c()
for (i in 1:length(failed_study_names)) {
  uninf_study_name <- c(uninf_study_name, rep(failed_study_names[i], 2))
  uninf_dose <- c(uninf_dose, rep(failed_doses[i], 2))
  uninf_time <- c(uninf_time, c(0, 30))
  uninf_viral_load <- c(uninf_viral_load, runif(2, min = 10^-0.8, max = 10^-0.4))
  uninf_id <- c(uninf_id, rep(i, 2))
}

# data frame for plotting failed infections
df_uninf <- data.frame(study_name = uninf_study_name,
                       dose = uninf_dose,
                       time = uninf_time,
                       viral_load = uninf_viral_load,
                       id = uninf_id) %>%
  mutate(facet = case_when(study_name == "JL" ~ "'JL -' ~ 10^6.9 ~ TCID[50]",
                           study_name == "BL" ~ "'BL -' ~ 10^5 ~ TCID[50]",
                           study_name == "LL16" ~ "'LL16 -' ~ 10^5 ~ TCID[50]",
                           study_name == "LL19" ~ "'LL19 -' ~ 10^5 ~ TCID[50]",
                           study_name == "YC" ~ "'YC -' ~ 10^5 ~ TCID[50]",
                           study_name == "NN" & dose == 1e+3 ~ "'NN -' ~ 10^3 ~ TCID[50]",
                           study_name == "NN" & dose == 1e+6 ~ "'NN -' ~ 10^6 ~ TCID[50]"))
# order facet variable for plotting
df_uninf$facet <- factor(df_uninf$facet, levels = c("'NN -' ~ 10^3 ~ TCID[50]",
                                                    "'BL -' ~ 10^5 ~ TCID[50]",
                                                    "'YC -' ~ 10^5 ~ TCID[50]",
                                                    "'LL16 -' ~ 10^5 ~ TCID[50]",
                                                    "'LL19 -' ~ 10^5 ~ TCID[50]",
                                                    "'NN -' ~ 10^6 ~ TCID[50]",
                                                    "'JL -' ~ 10^6.9 ~ TCID[50]"))

# data frame for plotting LODs
df_LOD <- data.frame(facet = c("'NN -' ~ 10^3 ~ TCID[50]",
                               "'BL -' ~ 10^5 ~ TCID[50]",
                               "'YC -' ~ 10^5 ~ TCID[50]",
                               "'LL16 -' ~ 10^5 ~ TCID[50]",
                               "'LL19 -' ~ 10^5 ~ TCID[50]",
                               "'NN -' ~ 10^6 ~ TCID[50]",
                               "'JL -' ~ 10^6.9 ~ TCID[50]"),
                     val = c(10^3.14, 10^1.56, 10^1, 10^1, 10^1, 10^3.14, 10^2))
df_LOD$facet <- factor(df_LOD$facet, levels = c("'NN -' ~ 10^3 ~ TCID[50]",
                                                "'BL -' ~ 10^5 ~ TCID[50]",
                                                "'YC -' ~ 10^5 ~ TCID[50]",
                                                "'LL16 -' ~ 10^5 ~ TCID[50]",
                                                "'LL19 -' ~ 10^5 ~ TCID[50]",
                                                "'NN -' ~ 10^6 ~ TCID[50]",
                                                "'JL -' ~ 10^6.9 ~ TCID[50]"))

# color settings
pred_color <- "#007991"

# add shading, zero line, LOD lines
p1 <- ggplot() +
  geom_ribbon(data = data.frame(x = seq(-5, 40, 1), y = 1, lower = 1e-2, upper = 1e+0),
              aes(ymin = lower, ymax = upper, x = x, y = y),
              alpha = 0.2) +
  geom_hline(yintercept = 1e+0) +
  geom_hline(data = df_LOD,
             aes(yintercept = val),
             alpha = 0.5, linetype = "dashed")

# add predicted trajectories
p1 <- p1 + geom_line(data = df_inf,
                     aes(x = time, y = viral_load, group = id),
                     alpha = 0.1, color = pred_color)

# add failed infections
p1 <- p1 + geom_line(data = df_uninf,
                     aes(x = time, y = viral_load, group = id),
                     alpha = 0.1, color = pred_color)

# final plot settings
p1 <- p1 +
  facet_grid(~ facet, labeller = label_parsed) +
  scale_x_continuous(breaks = seq(0, 40, 5)) +
  scale_y_continuous(trans = "log10",
                     breaks = c(1e+0, 1e+1, 1e+2, 1e+3, 1e+4, 1e+5, 1e+6, 1e+7),
                     labels = trans_format("log10", math_format(10^.x))) +
  coord_cartesian(xlim = c(0, 30), ylim = c(10^-0.8, 10^6)) +
  xlab("Day Post Infection") +
  ylab("Viral Load (RNA Copies/Sample)") +
  theme(legend.position = "none")

plot(p1)

## Fig S4B: posterior viral load trajectories ----------------------------------

set.seed(1234)

theme_set(theme_bw())

# extract parameters
draws <- post_fit %>%
  spread_draws(hit_prob,
               alpha_t, log_tp5[study_name],
               alpha_N, log_N5[study_name],
               g[study_name], d[study_name])

# sample a subset (100) of extracted parameters
n_draws <- length(unique(draws$.draw))
draws_per_study <- 100
sampled_draws <- sample(1:n_draws, size = draws_per_study)
draws <- draws %>%
  filter(.draw %in% sampled_draws)

# keep track of doses (TCID50) used in each study
study_doses <- list("BL" = c(10^5),
                    "LL16" = c(10^5),
                    "LL19" = c(10^5),
                    "NN" = c(10^3, 10^6),
                    "JL" = c(6*10^6/0.7),
                    "YC" = c(10^5))

# vectors to store successful trajectories
final_times <- c()
final_viral_loads <- c()
final_study_names <- c()
final_doses <- c()
final_ids <- c() # to keep track of individual rows (draws)

# vectors to store failed infections
failed_study_names <- c()
failed_doses <- c()

# times to compute trajectories
n_time_points <- 100
t_seq <- seq(from = 0, to = 35, length.out = n_time_points)

# compute and save trajectories for each draw
for (i in 1:nrow(draws)){
  
  # extract parameters from draw
  g <- draws[i,]$g
  d <- draws[i,]$d
  hit_prob <- draws[i,]$hit_prob
  alpha_t <- draws[i,]$alpha_t
  alpha_N <- draws[i,]$alpha_N
  log_tp5 <- draws[i,]$log_tp5
  log_N5 <- draws[i,]$log_N5
  study_name <- as.character(draws[i,]$study_name)
  
  for (dose in study_doses[[study_name]]){
    
    # simulate infection, chance that monkey is not infected
    inf_prob <- DR_prob(hit_prob, dose)
    # keep track of failed infections
    if (runif(1) > inf_prob) {
      failed_study_names <- c(failed_study_names, study_name)
      failed_doses <- c(failed_doses, dose)
      next
    }
    
    # if infected, compute predicted trajectory
    t_p <- peak_time(hit_prob, dose, alpha_t, log_tp5)
    logN_0 <- log_init_RNA(hit_prob, dose, alpha_N, log_N5)
    viral_loads <- unlist(lapply(t_seq,
                                 exp_growth_decay,
                                 g_i = g,
                                 d_i = d,
                                 t_p = t_p,
                                 logN_0 = logN_0))
    viral_loads <- exp(viral_loads)
    
    # save trajectory
    final_times <- c(final_times, t_seq)
    final_viral_loads <- c(final_viral_loads, viral_loads)
    final_study_names <- c(final_study_names, rep(study_name, n_time_points))
    final_doses <- c(final_doses, rep(dose, n_time_points))
    final_ids <- c(final_ids, rep(i, n_time_points))
  }
  
}

# data frame for plotting successful trajectories
df_inf <- data.frame(time = final_times,
                     viral_load = final_viral_loads,
                     study_name = final_study_names,
                     dose = final_doses,
                     id = final_ids) %>%
  filter(viral_load >= 1) %>% # remove points below 10^0
  mutate(facet = case_when(study_name == "JL" ~ "'JL -' ~ 10^6.9 ~ TCID[50]",
                           study_name == "BL" ~ "'BL -' ~ 10^5 ~ TCID[50]",
                           study_name == "LL16" ~ "'LL16 -' ~ 10^5 ~ TCID[50]",
                           study_name == "LL19" ~ "'LL19 -' ~ 10^5 ~ TCID[50]",
                           study_name == "YC" ~ "'YC -' ~ 10^5 ~ TCID[50]",
                           study_name == "NN" & dose == 1e+3 ~ "'NN -' ~ 10^3 ~ TCID[50]",
                           study_name == "NN" & dose == 1e+6 ~ "'NN -' ~ 10^6 ~ TCID[50]"))
# order facet variable for plotting
df_inf$facet <- factor(df_inf$facet, levels = c("'NN -' ~ 10^3 ~ TCID[50]",
                                                "'BL -' ~ 10^5 ~ TCID[50]",
                                                "'YC -' ~ 10^5 ~ TCID[50]",
                                                "'LL16 -' ~ 10^5 ~ TCID[50]",
                                                "'LL19 -' ~ 10^5 ~ TCID[50]",
                                                "'NN -' ~ 10^6 ~ TCID[50]",
                                                "'JL -' ~ 10^6.9 ~ TCID[50]"))

# create trajectories for failed infections
# lines that will be plotted below zero to indicate failed infections
uninf_study_name <- c()
uninf_dose <- c()
uninf_time <- c()
uninf_viral_load <- c()
uninf_id <- c()
for (i in 1:length(failed_study_names)) {
  uninf_study_name <- c(uninf_study_name, rep(failed_study_names[i], 2))
  uninf_dose <- c(uninf_dose, rep(failed_doses[i], 2))
  uninf_time <- c(uninf_time, c(0, 30))
  uninf_viral_load <- c(uninf_viral_load, runif(2, min = 10^-0.8, max = 10^-0.4))
  uninf_id <- c(uninf_id, rep(i, 2))
}

# data frame for plotting failed infections
df_uninf <- data.frame(study_name = uninf_study_name,
                       dose = uninf_dose,
                       time = uninf_time,
                       viral_load = uninf_viral_load,
                       id = uninf_id) %>%
  mutate(facet = case_when(study_name == "JL" ~ "'JL -' ~ 10^6.9 ~ TCID[50]",
                           study_name == "BL" ~ "'BL -' ~ 10^5 ~ TCID[50]",
                           study_name == "LL16" ~ "'LL16 -' ~ 10^5 ~ TCID[50]",
                           study_name == "LL19" ~ "'LL19 -' ~ 10^5 ~ TCID[50]",
                           study_name == "YC" ~ "'YC -' ~ 10^5 ~ TCID[50]",
                           study_name == "NN" & dose == 1e+3 ~ "'NN -' ~ 10^3 ~ TCID[50]",
                           study_name == "NN" & dose == 1e+6 ~ "'NN -' ~ 10^6 ~ TCID[50]"))
# order facet variable for plotting
df_uninf$facet <- factor(df_uninf$facet, levels = c("'NN -' ~ 10^3 ~ TCID[50]",
                                                    "'BL -' ~ 10^5 ~ TCID[50]",
                                                    "'YC -' ~ 10^5 ~ TCID[50]",
                                                    "'LL16 -' ~ 10^5 ~ TCID[50]",
                                                    "'LL19 -' ~ 10^5 ~ TCID[50]",
                                                    "'NN -' ~ 10^6 ~ TCID[50]",
                                                    "'JL -' ~ 10^6.9 ~ TCID[50]"))

# data frame for plotting LODs
df_LOD <- data.frame(facet = c("'NN -' ~ 10^3 ~ TCID[50]",
                               "'BL -' ~ 10^5 ~ TCID[50]",
                               "'YC -' ~ 10^5 ~ TCID[50]",
                               "'LL16 -' ~ 10^5 ~ TCID[50]",
                               "'LL19 -' ~ 10^5 ~ TCID[50]",
                               "'NN -' ~ 10^6 ~ TCID[50]",
                               "'JL -' ~ 10^6.9 ~ TCID[50]"),
                     val = c(10^3.14, 10^1.56, 10^1, 10^1, 10^1, 10^3.14, 10^2))
df_LOD$facet <- factor(df_LOD$facet, levels = c("'NN -' ~ 10^3 ~ TCID[50]",
                                                "'BL -' ~ 10^5 ~ TCID[50]",
                                                "'YC -' ~ 10^5 ~ TCID[50]",
                                                "'LL16 -' ~ 10^5 ~ TCID[50]",
                                                "'LL19 -' ~ 10^5 ~ TCID[50]",
                                                "'NN -' ~ 10^6 ~ TCID[50]",
                                                "'JL -' ~ 10^6.9 ~ TCID[50]"))

# color settings
pred_color <- "#007991"

# add shading, zero line, LOD lines
p2 <- ggplot() +
  geom_ribbon(data = data.frame(x = seq(-5, 40, 1), y = 1, lower = 1e-2, upper = 1e+0),
              aes(ymin = lower, ymax = upper, x = x, y = y),
              alpha = 0.2) +
  geom_hline(yintercept = 1e+0) +
  geom_hline(data = df_LOD,
             aes(yintercept = val),
             alpha = 0.5, linetype = "dashed")

# add predicted trajectories
p2 <- p2 + geom_line(data = df_inf,
                     aes(x = time, y = viral_load, group = id),
                     alpha = 0.1, color = pred_color)

# add failed infections
p2 <- p2 + geom_line(data = df_uninf,
                     aes(x = time, y = viral_load, group = id),
                     alpha = 0.1, color = pred_color)

# final plot settings
p2 <- p2 +
  facet_grid(~ facet, labeller = label_parsed) +
  scale_x_continuous(breaks = seq(0, 40, 5)) +
  scale_y_continuous(trans = "log10",
                     breaks = c(1e+0, 1e+1, 1e+2, 1e+3, 1e+4, 1e+5, 1e+6, 1e+7),
                     labels = trans_format("log10", math_format(10^.x))) +
  coord_cartesian(xlim = c(0, 30), ylim = c(10^-0.8, 10^6)) +
  xlab("Day Post Infection") +
  ylab("Viral Load (RNA Copies/Sample)") +
  theme(legend.position = "none")

plot(p2)

## Fig S4C: prior predictive simulation ----------------------------------------

set.seed(1234)

theme_set(theme_bw())

# plot settings
draws_per_study <- 100
points_per_line <- 5

# extract parameters
draws <- prior_fit %>%
  spread_draws(hit_prob,
               alpha_t, log_tp5[study_name],
               alpha_N, log_N5[study_name],
               g[study_name], d[study_name],
               sigma_s[study_name])

# sample a subset (100) of extracted parameters
n_draws <- length(unique(draws$.draw))
sampled_draws <- sample(1:n_draws, size = draws_per_study)
draws <- draws %>%
  filter(.draw %in% sampled_draws)

# keep track of doses (TCID50) used in each study
study_doses <- list("BL" = c(10^5),
                    "LL16" = c(10^5),
                    "LL19" = c(10^5),
                    "NN" = c(10^3, 10^6),
                    "JL" = c(6*10^6/0.7),
                    "YC" = c(10^5))

# vectors to store successful trajectories
final_times <- c()
final_viral_loads <- c()
final_study_names <- c()
final_doses <- c()
final_ids <- c() # to keep track of individual rows (draws)

# vectors to store failed infections
failed_study_names <- c()
failed_doses <- c()

# compute and save simulated data for each draw
for (i in 1:nrow(draws)){
  
  # extract parameters from draw
  g <- draws[i,]$g
  d <- draws[i,]$d
  hit_prob <- draws[i,]$hit_prob
  alpha_t <- draws[i,]$alpha_t
  alpha_N <- draws[i,]$alpha_N
  log_tp5 <- draws[i,]$log_tp5
  log_N5 <- draws[i,]$log_N5
  study_name <- as.character(draws[i,]$study_name)
  sigma_s <- draws[i,]$sigma_s
  
  for (dose in study_doses[[study_name]]){
    
    # simulate infection, chance that monkey is not infected
    inf_prob <- DR_prob(hit_prob, dose)
    # keep track of failed infections
    if (runif(1) > inf_prob) {
      failed_study_names <- c(failed_study_names, study_name)
      failed_doses <- c(failed_doses, dose)
      next
    }
    
    # if infected, compute predicted trajectory
    t_p <- peak_time(hit_prob, dose, alpha_t, log_tp5)
    logN_0 <- log_init_RNA(hit_prob, dose, alpha_N, log_N5)
    # pick times randomly in range 0 to 30
    t_sample <- sample(1:30, size = points_per_line)
    # simulate data from predicted trajectory
    mean_viral_loads <- unlist(lapply(t_sample,
                                      exp_growth_decay,
                                      g_i = g,
                                      d_i = d,
                                      t_p = t_p,
                                      logN_0 = logN_0))
    viral_loads <- rnorm(points_per_line,
                         mean = mean_viral_loads,
                         sd = sigma_s)
    viral_loads <- exp(viral_loads)
    
    # save simulated data
    final_times <- c(final_times, t_sample)
    final_viral_loads <- c(final_viral_loads, viral_loads)
    final_study_names <- c(final_study_names, rep(study_name, points_per_line))
    final_doses <- c(final_doses, rep(dose, points_per_line))
    final_ids <- c(final_ids, rep(i, points_per_line))
  }
  
}

# data frame to plot successful infections
df_inf <- data.frame(time = final_times,
                     viral_load = final_viral_loads,
                     study_name = final_study_names,
                     dose = final_doses,
                     id = final_ids) %>%
  filter(viral_load > 1) %>%
  mutate(facet = case_when(study_name == "JL" ~ "'JL -' ~ 10^6.9 ~ TCID[50]",
                           study_name == "BL" ~ "'BL -' ~ 10^5 ~ TCID[50]",
                           study_name == "LL16" ~ "'LL16 -' ~ 10^5 ~ TCID[50]",
                           study_name == "LL19" ~ "'LL19 -' ~ 10^5 ~ TCID[50]",
                           study_name == "YC" ~ "'YC -' ~ 10^5 ~ TCID[50]",
                           study_name == "NN" & dose == 1e+3 ~ "'NN -' ~ 10^3 ~ TCID[50]",
                           study_name == "NN" & dose == 1e+6 ~ "'NN -' ~ 10^6 ~ TCID[50]"))
df_inf$facet <- factor(df_inf$facet, levels = c("'NN -' ~ 10^3 ~ TCID[50]",
                                                "'BL -' ~ 10^5 ~ TCID[50]",
                                                "'YC -' ~ 10^5 ~ TCID[50]",
                                                "'LL16 -' ~ 10^5 ~ TCID[50]",
                                                "'LL19 -' ~ 10^5 ~ TCID[50]",
                                                "'NN -' ~ 10^6 ~ TCID[50]",
                                                "'JL -' ~ 10^6.9 ~ TCID[50]"))

# create plotting data for failed infections
fail_time <- c()
fail_viral_load <- c()
fail_study_name <- c()
fail_dose <- c()
fail_id <- c()
for (i in 1:length(failed_study_names)) {
  fail_time <- c(fail_time, runif(5, min = 0, max = 30))
  fail_viral_load <- c(fail_viral_load, runif(5, min = 10^-0.9, max = 10^-0.3))
  fail_study_name <- c(fail_study_name, rep(failed_study_names[i], 5))
  fail_dose <- c(fail_dose, rep(failed_doses[i], 5))
  fail_id <- c(fail_id, rep(i, 5))
}

# data frame to plot failed infections
df_uninf <- data.frame(time = fail_time,
                       viral_load = fail_viral_load,
                       study_name = fail_study_name,
                       dose = fail_dose,
                       id = fail_id) %>%
  mutate(facet = case_when(study_name == "JL" ~ "'JL -' ~ 10^6.9 ~ TCID[50]",
                           study_name == "BL" ~ "'BL -' ~ 10^5 ~ TCID[50]",
                           study_name == "LL16" ~ "'LL16 -' ~ 10^5 ~ TCID[50]",
                           study_name == "LL19" ~ "'LL19 -' ~ 10^5 ~ TCID[50]",
                           study_name == "YC" ~ "'YC -' ~ 10^5 ~ TCID[50]",
                           study_name == "NN" & dose == 1e+3 ~ "'NN -' ~ 10^3 ~ TCID[50]",
                           study_name == "NN" & dose == 1e+6 ~ "'NN -' ~ 10^6 ~ TCID[50]"))
df_uninf$facet <- factor(df_uninf$facet, levels = c("'NN -' ~ 10^3 ~ TCID[50]",
                                                    "'BL -' ~ 10^5 ~ TCID[50]",
                                                    "'YC -' ~ 10^5 ~ TCID[50]",
                                                    "'LL16 -' ~ 10^5 ~ TCID[50]",
                                                    "'LL19 -' ~ 10^5 ~ TCID[50]",
                                                    "'NN -' ~ 10^6 ~ TCID[50]",
                                                    "'JL -' ~ 10^6.9 ~ TCID[50]"))

# create data frame for plotting observed data
indiv_data <- data %>%
  filter(log10_viral_load != 100 & data_type != "nasal swab") %>%
  mutate(below_LOD = case_when(log10_viral_load == -100 ~ 1,
                               TRUE ~ 0)) %>%
  mutate(facet = case_when(study_name == "JL" ~ "'JL -' ~ 10^6.9 ~ TCID[50]",
                           study_name == "BL" ~ "'BL -' ~ 10^5 ~ TCID[50]",
                           study_name == "LL16" ~ "'LL16 -' ~ 10^5 ~ TCID[50]",
                           study_name == "LL19" ~ "'LL19 -' ~ 10^5 ~ TCID[50]",
                           study_name == "YC" ~ "'YC -' ~ 10^5 ~ TCID[50]",
                           study_name == "NN" & TCID_dose == 1e+3 ~ "'NN -' ~ 10^3 ~ TCID[50]",
                           study_name == "NN" & TCID_dose == 1e+6 ~ "'NN -' ~ 10^6 ~ TCID[50]"))
indiv_data$facet <- factor(indiv_data$facet, levels = c("'NN -' ~ 10^3 ~ TCID[50]",
                                                        "'BL -' ~ 10^5 ~ TCID[50]",
                                                        "'YC -' ~ 10^5 ~ TCID[50]",
                                                        "'LL16 -' ~ 10^5 ~ TCID[50]",
                                                        "'LL19 -' ~ 10^5 ~ TCID[50]",
                                                        "'NN -' ~ 10^6 ~ TCID[50]",
                                                        "'JL -' ~ 10^6.9 ~ TCID[50]"))
# manually set points below LOD
indiv_data <- indiv_data %>%
  mutate(log10_viral_load = case_when(individual == "JL_CM-03" & log10_viral_load == -100 ~ -0.4,
                                      individual == "JL_CM-04" & log10_viral_load == -100 ~ -0.8,
                                      individual == "NN_CM-01" & log10_viral_load == -100 ~ -0.4,
                                      individual == "NN_CM-02" & log10_viral_load == -100 ~ -0.8,
                                      individual == "NN_CM-03" & log10_viral_load == -100 ~ -0.4,
                                      individual == "NN_CM-04" & log10_viral_load == -100 ~ -0.8,
                                      individual == "LL_RM-14" & log10_viral_load == -100 ~ -0.6,
                                      TRUE ~ log10_viral_load))

# data frame for plotting LODs
df_LOD <- data.frame(facet = c("'NN -' ~ 10^3 ~ TCID[50]",
                               "'BL -' ~ 10^5 ~ TCID[50]",
                               "'YC -' ~ 10^5 ~ TCID[50]",
                               "'LL16 -' ~ 10^5 ~ TCID[50]",
                               "'LL19 -' ~ 10^5 ~ TCID[50]",
                               "'NN -' ~ 10^6 ~ TCID[50]",
                               "'JL -' ~ 10^6.9 ~ TCID[50]"),
                     val = c(10^3.14, 10^1.56, 10^1, 10^1, 10^1, 10^3.14, 10^2))
df_LOD$facet <- factor(df_LOD$facet, levels = c("'NN -' ~ 10^3 ~ TCID[50]",
                                                "'BL -' ~ 10^5 ~ TCID[50]",
                                                "'YC -' ~ 10^5 ~ TCID[50]",
                                                "'LL16 -' ~ 10^5 ~ TCID[50]",
                                                "'LL19 -' ~ 10^5 ~ TCID[50]",
                                                "'NN -' ~ 10^6 ~ TCID[50]",
                                                "'JL -' ~ 10^6.9 ~ TCID[50]"))

# add shading, zero line, LOD lines
p3 <- ggplot() +
  geom_ribbon(data = data.frame(x = seq(-5, 40, 1), y = 1, lower = 1e-2, upper = 1e+0),
              aes(ymin = lower, ymax = upper, x = x, y = y),
              alpha = 0.2) +
  geom_hline(yintercept = 1e+0) +
  geom_hline(data = df_LOD,
             aes(yintercept = val),
             alpha = 0.5, linetype = "dashed")

# add simulated successful infections
p3 <- p3 +
  geom_line(data = df_inf,
            aes(x = time, y = viral_load, group = id),
            alpha = 0.1, color = "#007991", linewidth = 0.3) +
  geom_point(data = df_inf,
             aes(x = time, y = viral_load, group = id),
             alpha = 0.1, color = "#007991", stroke = 0, size = 2)

# add simulated failed infections
p3 <- p3 +
  geom_line(data = df_uninf,
            aes(x = time, y = viral_load, group = id),
            alpha = 0.1, color = "#007991", linewidth = 0.3) +
  geom_point(data = df_uninf,
             aes(x = time, y = viral_load, group = id),
             alpha = 0.1, color = "#007991", stroke = 0, size = 2)

# add observed data
p3 <- p3 +
  geom_line(data = indiv_data,
            aes(x = day_post_infection, y = 10^log10_viral_load, group = individual),
            alpha = 0.5,
            linewidth = 0.3,
            color = "#D7263D") +
  geom_point(data = indiv_data,
             aes(x = day_post_infection, y = 10^log10_viral_load, shape = as.factor(below_LOD)),
             alpha = 0.5,
             size = 2,
             stroke = 0,
             fill = "#D7263D") +
  scale_shape_manual(values = c(21, 25), guide = "none")

# final plot settings
p3 <- p3 +
  facet_grid(~ facet, labeller = label_parsed) +
  scale_x_continuous(breaks = seq(0, 40, 5)) +
  scale_y_continuous(trans = "log10",
                     breaks = c(1e+0, 1e+1, 1e+2, 1e+3, 1e+4, 1e+5, 1e+6, 1e+7),
                     labels = trans_format("log10", math_format(10^.x))) +
  coord_cartesian(xlim = c(0, 30), ylim = c(10^-0.8, 10^6)) +
  xlab("Day Post Infection") +
  ylab("Viral Load (RNA Copies/Sample)")

plot(p3)

## Fig S4D: posterior predictive simulation ------------------------------------

set.seed(1234)

theme_set(theme_bw())

# plot settings
draws_per_study <- 100
points_per_line <- 5

# extract parameters
draws <- post_fit %>%
  spread_draws(hit_prob,
               alpha_t, log_tp5[study_name],
               alpha_N, log_N5[study_name],
               g[study_name], d[study_name],
               sigma_s[study_name])

# sample a subset (100) of extracted parameters
n_draws <- length(unique(draws$.draw))
sampled_draws <- sample(1:n_draws, size = draws_per_study)
draws <- draws %>%
  filter(.draw %in% sampled_draws)

# keep track of doses (TCID50) used in each study
study_doses <- list("BL" = c(10^5),
                    "LL16" = c(10^5),
                    "LL19" = c(10^5),
                    "NN" = c(10^3, 10^6),
                    "JL" = c(6*10^6/0.7),
                    "YC" = c(10^5))

# vectors to store successful trajectories
final_times <- c()
final_viral_loads <- c()
final_study_names <- c()
final_doses <- c()
final_ids <- c() # to keep track of individual rows (draws)

# vectors to store failed infections
failed_study_names <- c()
failed_doses <- c()

# compute and save simulated data for each draw
for (i in 1:nrow(draws)){
  
  # extract parameters from draw
  g <- draws[i,]$g
  d <- draws[i,]$d
  hit_prob <- draws[i,]$hit_prob
  alpha_t <- draws[i,]$alpha_t
  alpha_N <- draws[i,]$alpha_N
  log_tp5 <- draws[i,]$log_tp5
  log_N5 <- draws[i,]$log_N5
  study_name <- as.character(draws[i,]$study_name)
  sigma_s <- draws[i,]$sigma_s
  
  for (dose in study_doses[[study_name]]){
    
    # simulate infection, chance that monkey is not infected
    inf_prob <- DR_prob(hit_prob, dose)
    # keep track of failed infections
    if (runif(1) > inf_prob) {
      failed_study_names <- c(failed_study_names, study_name)
      failed_doses <- c(failed_doses, dose)
      next
    }
    
    # if infected, compute predicted trajectory
    t_p <- peak_time(hit_prob, dose, alpha_t, log_tp5)
    logN_0 <- log_init_RNA(hit_prob, dose, alpha_N, log_N5)
    # pick times randomly in range 0 to 30
    t_sample <- sample(1:30, size = points_per_line)
    # simulate data from predicted trajectory
    mean_viral_loads <- unlist(lapply(t_sample,
                                      exp_growth_decay,
                                      g_i = g,
                                      d_i = d,
                                      t_p = t_p,
                                      logN_0 = logN_0))
    viral_loads <- rnorm(points_per_line,
                         mean = mean_viral_loads,
                         sd = sigma_s)
    viral_loads <- exp(viral_loads)
    
    # save simulated data
    final_times <- c(final_times, t_sample)
    final_viral_loads <- c(final_viral_loads, viral_loads)
    final_study_names <- c(final_study_names, rep(study_name, points_per_line))
    final_doses <- c(final_doses, rep(dose, points_per_line))
    final_ids <- c(final_ids, rep(i, points_per_line))
  }
  
}

# data frame to plot successful infections
df_inf <- data.frame(time = final_times,
                     viral_load = final_viral_loads,
                     study_name = final_study_names,
                     dose = final_doses,
                     id = final_ids) %>%
  filter(viral_load > 1) %>%
  mutate(facet = case_when(study_name == "JL" ~ "'JL -' ~ 10^6.9 ~ TCID[50]",
                           study_name == "BL" ~ "'BL -' ~ 10^5 ~ TCID[50]",
                           study_name == "LL16" ~ "'LL16 -' ~ 10^5 ~ TCID[50]",
                           study_name == "LL19" ~ "'LL19 -' ~ 10^5 ~ TCID[50]",
                           study_name == "YC" ~ "'YC -' ~ 10^5 ~ TCID[50]",
                           study_name == "NN" & dose == 1e+3 ~ "'NN -' ~ 10^3 ~ TCID[50]",
                           study_name == "NN" & dose == 1e+6 ~ "'NN -' ~ 10^6 ~ TCID[50]"))
df_inf$facet <- factor(df_inf$facet, levels = c("'NN -' ~ 10^3 ~ TCID[50]",
                                                "'BL -' ~ 10^5 ~ TCID[50]",
                                                "'YC -' ~ 10^5 ~ TCID[50]",
                                                "'LL16 -' ~ 10^5 ~ TCID[50]",
                                                "'LL19 -' ~ 10^5 ~ TCID[50]",
                                                "'NN -' ~ 10^6 ~ TCID[50]",
                                                "'JL -' ~ 10^6.9 ~ TCID[50]"))

# create plotting data for failed infections
fail_time <- c()
fail_viral_load <- c()
fail_study_name <- c()
fail_dose <- c()
fail_id <- c()
for (i in 1:length(failed_study_names)) {
  fail_time <- c(fail_time, runif(5, min = 0, max = 30))
  fail_viral_load <- c(fail_viral_load, runif(5, min = 10^-0.9, max = 10^-0.3))
  fail_study_name <- c(fail_study_name, rep(failed_study_names[i], 5))
  fail_dose <- c(fail_dose, rep(failed_doses[i], 5))
  fail_id <- c(fail_id, rep(i, 5))
}

# data frame to plot failed infections
df_uninf <- data.frame(time = fail_time,
                       viral_load = fail_viral_load,
                       study_name = fail_study_name,
                       dose = fail_dose,
                       id = fail_id) %>%
  mutate(facet = case_when(study_name == "JL" ~ "'JL -' ~ 10^6.9 ~ TCID[50]",
                           study_name == "BL" ~ "'BL -' ~ 10^5 ~ TCID[50]",
                           study_name == "LL16" ~ "'LL16 -' ~ 10^5 ~ TCID[50]",
                           study_name == "LL19" ~ "'LL19 -' ~ 10^5 ~ TCID[50]",
                           study_name == "YC" ~ "'YC -' ~ 10^5 ~ TCID[50]",
                           study_name == "NN" & dose == 1e+3 ~ "'NN -' ~ 10^3 ~ TCID[50]",
                           study_name == "NN" & dose == 1e+6 ~ "'NN -' ~ 10^6 ~ TCID[50]"))
df_uninf$facet <- factor(df_uninf$facet, levels = c("'NN -' ~ 10^3 ~ TCID[50]",
                                                    "'BL -' ~ 10^5 ~ TCID[50]",
                                                    "'YC -' ~ 10^5 ~ TCID[50]",
                                                    "'LL16 -' ~ 10^5 ~ TCID[50]",
                                                    "'LL19 -' ~ 10^5 ~ TCID[50]",
                                                    "'NN -' ~ 10^6 ~ TCID[50]",
                                                    "'JL -' ~ 10^6.9 ~ TCID[50]"))

# create data frame for plotting observed data
indiv_data <- data %>%
  filter(log10_viral_load != 100 & data_type != "nasal swab") %>%
  mutate(below_LOD = case_when(log10_viral_load == -100 ~ 1,
                               TRUE ~ 0)) %>%
  mutate(facet = case_when(study_name == "JL" ~ "'JL -' ~ 10^6.9 ~ TCID[50]",
                           study_name == "BL" ~ "'BL -' ~ 10^5 ~ TCID[50]",
                           study_name == "LL16" ~ "'LL16 -' ~ 10^5 ~ TCID[50]",
                           study_name == "LL19" ~ "'LL19 -' ~ 10^5 ~ TCID[50]",
                           study_name == "YC" ~ "'YC -' ~ 10^5 ~ TCID[50]",
                           study_name == "NN" & TCID_dose == 1e+3 ~ "'NN -' ~ 10^3 ~ TCID[50]",
                           study_name == "NN" & TCID_dose == 1e+6 ~ "'NN -' ~ 10^6 ~ TCID[50]"))
indiv_data$facet <- factor(indiv_data$facet, levels = c("'NN -' ~ 10^3 ~ TCID[50]",
                                                        "'BL -' ~ 10^5 ~ TCID[50]",
                                                        "'YC -' ~ 10^5 ~ TCID[50]",
                                                        "'LL16 -' ~ 10^5 ~ TCID[50]",
                                                        "'LL19 -' ~ 10^5 ~ TCID[50]",
                                                        "'NN -' ~ 10^6 ~ TCID[50]",
                                                        "'JL -' ~ 10^6.9 ~ TCID[50]"))
# manually set points below LOD
indiv_data <- indiv_data %>%
  mutate(log10_viral_load = case_when(individual == "JL_CM-03" & log10_viral_load == -100 ~ -0.4,
                                      individual == "JL_CM-04" & log10_viral_load == -100 ~ -0.8,
                                      individual == "NN_CM-01" & log10_viral_load == -100 ~ -0.4,
                                      individual == "NN_CM-02" & log10_viral_load == -100 ~ -0.8,
                                      individual == "NN_CM-03" & log10_viral_load == -100 ~ -0.4,
                                      individual == "NN_CM-04" & log10_viral_load == -100 ~ -0.8,
                                      individual == "LL_RM-14" & log10_viral_load == -100 ~ -0.6,
                                      TRUE ~ log10_viral_load))

# data frame for plotting LODs
df_LOD <- data.frame(facet = c("'NN -' ~ 10^3 ~ TCID[50]",
                               "'BL -' ~ 10^5 ~ TCID[50]",
                               "'YC -' ~ 10^5 ~ TCID[50]",
                               "'LL16 -' ~ 10^5 ~ TCID[50]",
                               "'LL19 -' ~ 10^5 ~ TCID[50]",
                               "'NN -' ~ 10^6 ~ TCID[50]",
                               "'JL -' ~ 10^6.9 ~ TCID[50]"),
                     val = c(10^3.14, 10^1.56, 10^1, 10^1, 10^1, 10^3.14, 10^2))
df_LOD$facet <- factor(df_LOD$facet, levels = c("'NN -' ~ 10^3 ~ TCID[50]",
                                                "'BL -' ~ 10^5 ~ TCID[50]",
                                                "'YC -' ~ 10^5 ~ TCID[50]",
                                                "'LL16 -' ~ 10^5 ~ TCID[50]",
                                                "'LL19 -' ~ 10^5 ~ TCID[50]",
                                                "'NN -' ~ 10^6 ~ TCID[50]",
                                                "'JL -' ~ 10^6.9 ~ TCID[50]"))

# add shading, zero line, LOD lines
p4 <- ggplot() +
  geom_ribbon(data = data.frame(x = seq(-5, 40, 1), y = 1, lower = 1e-2, upper = 1e+0),
              aes(ymin = lower, ymax = upper, x = x, y = y),
              alpha = 0.2) +
  geom_hline(yintercept = 1e+0) +
  geom_hline(data = df_LOD,
             aes(yintercept = val),
             alpha = 0.5, linetype = "dashed")

# add simulated successful infections
p4 <- p4 +
  geom_line(data = df_inf,
            aes(x = time, y = viral_load, group = id),
            alpha = 0.1, color = "#007991", linewidth = 0.3) +
  geom_point(data = df_inf,
             aes(x = time, y = viral_load, group = id),
             alpha = 0.1, color = "#007991", stroke = 0, size = 2)

# add simulated failed infections
p4 <- p4 +
  geom_line(data = df_uninf,
            aes(x = time, y = viral_load, group = id),
            alpha = 0.1, color = "#007991", linewidth = 0.3) +
  geom_point(data = df_uninf,
             aes(x = time, y = viral_load, group = id),
             alpha = 0.1, color = "#007991", stroke = 0, size = 2)

# add observed data
p4 <- p4 +
  geom_line(data = indiv_data,
            aes(x = day_post_infection, y = 10^log10_viral_load, group = individual),
            alpha = 0.5,
            linewidth = 0.3,
            color = "#D7263D") +
  geom_point(data = indiv_data,
             aes(x = day_post_infection, y = 10^log10_viral_load, shape = as.factor(below_LOD)),
             alpha = 0.5,
             size = 2,
             stroke = 0,
             fill = "#D7263D") +
  scale_shape_manual(values = c(21, 25), guide = "none")

# final plot settings
p4 <- p4 +
  facet_grid(~ facet, labeller = label_parsed) +
  scale_x_continuous(breaks = seq(0, 40, 5)) +
  scale_y_continuous(trans = "log10",
                     breaks = c(1e+0, 1e+1, 1e+2, 1e+3, 1e+4, 1e+5, 1e+6, 1e+7),
                     labels = trans_format("log10", math_format(10^.x))) +
  coord_cartesian(xlim = c(0, 30), ylim = c(10^-0.8, 10^6)) +
  xlab("Day Post Infection") +
  ylab("Viral Load (RNA Copies/Sample")

plot(p4)

## Fig S4: final combined plot -------------------------------------------------

fig <- (p1 / p2 / p3 / p4) + plot_annotation(tag_levels = "A")

## Save plot -------------------------------------------------------------------

ggsave('outputs/figures/S4 Fig.tiff', fig, width = 16, height = 12)
