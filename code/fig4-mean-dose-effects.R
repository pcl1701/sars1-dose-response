rm(list = ls())

## Install and load necessary packages -----------------------------------------

req_pkgs <- c("dplyr", "tidyr", "tidybayes", "ggplot2", "scales", "patchwork")
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

## Overall settings ------------------------------------------------------------

set.seed(1234)
theme_set(theme_bw())

## Fig 4A: dose-dependence of initial viral load -------------------------------

# extract parameters
draws <- fit %>%
  spread_draws(hit_prob, alpha_N, mean_log_N5)

# sample rows
n_draws <- 100
draw_ids <- sample(1:nrow(draws), size = n_draws, replace = FALSE)
slope_draws <- draws[draw_ids,]

# generate range of inoculation doses
log10_doses <- seq(from = 0, to = 8, length.out = 50)
doses <- 10^log10_doses
n_doses <- length(doses)

# vectors to store plotting values
final_doses <- c()
final_init_RNAs <- c()
final_ids <- c() # to keep track of individual rows (draws)

# for each slope draw (alpha_N, beta_N)
# generate all points in initial viral load vs. dose curve
for (i in 1:nrow(slope_draws)){
  
  # compute initial viral loads for all doses
  # for a single row (draw), (alpha_N, beta_N)
  log_init_RNAs <- log_init_RNA(slope_draws[i,]$hit_prob,
                                doses,
                                slope_draws[i,]$alpha_N,
                                slope_draws[i,]$mean_log_N5)
  init_RNAs <- exp(log_init_RNAs)
  
  # save doses, initial viral loads, row ids
  final_doses <- c(final_doses, doses)
  final_init_RNAs <- c(final_init_RNAs, init_RNAs)
  final_ids <- c(final_ids, rep(i, n_doses))
  
}

# data frame for plotting initial viral loads vs doses
df_N <- data.frame(dose = final_doses,
                   init_RNA = final_init_RNAs,
                   id = final_ids)

# plot results
pXa <- df_N %>% ggplot(aes(x = dose, y = init_RNA, group = id)) +
  geom_line(alpha = 0.2,
            color = "#488B67") +
  geom_vline(xintercept = c(10^3, 10^5, 10^6, 10^6.9),
             linetype = "dashed",
             color = "red") +
  xlab(bquote("Inoculation Dose (" * TCID[50] * ")")) +
  ylab("Initial Viral Load") +
  scale_x_continuous(trans = "log10",
                     breaks = c(1e+0, 1e+1, 1e+2, 1e+3, 1e+4, 1e+5, 1e+6, 1e+7),
                     labels = trans_format("log10", math_format(10^.x))) +
  scale_y_continuous(trans = "log10",
                     breaks = c(1e-2, 1e-1, 1e+0, 1e+1, 1e+2, 1e+3, 1e+4, 1e+5),
                     labels = trans_format("log10", math_format(10^.x))) +
  coord_cartesian(xlim = c(1e+2, 1e+7), ylim=c(1e-1, 1e+5))

plot(pXa)

## Fig 4B: dose-dependence of time of peak viral load --------------------------

# extract parameters
draws <- fit %>%
  spread_draws(hit_prob, alpha_t, mean_log_tp5)

# sample rows
slope_draws <- draws[draw_ids,]

# generate range of inoculation doses
log10_doses <- seq(from = 0, to = 8, length.out = 50)
doses <- 10^log10_doses
n_doses <- length(doses)

# vectors to store plotting values
final_doses <- c()
final_peak_times <- c()
final_ids <- c() # to keep track of individual rows

# for each slope draw (alpha_t, beta_t)
# generate all points in time of peak viral load vs. dose curve
for (i in 1:nrow(slope_draws)){
  
  # compute peak times for all doses
  # for a single row (draw), (alpha_t, beta_t)
  peak_times <- peak_time(slope_draws[i,]$hit_prob,
                          doses,
                          slope_draws[i,]$alpha_t,
                          slope_draws[i,]$mean_log_tp5)
  
  # save doses, times of peak viral load, row ids
  final_doses <- c(final_doses, doses)
  final_peak_times <- c(final_peak_times, peak_times)
  final_ids <- c(final_ids, rep(i, n_doses))
  
}

# data frame for plotting times of peak viral load vs doses
df_t <- data.frame(dose = final_doses,
                   peak_time = final_peak_times,
                   id = final_ids)

# plot results
pXb <- df_t %>% ggplot(aes(x = dose, y = peak_time, group = id)) +
  geom_line(alpha = 0.2,
            color = "#725AA6") +
  geom_vline(xintercept = c(10^3, 10^5, 10^6, 10^6.9),
             linetype = "dashed",
             color = "red") +
  xlab(bquote("Inoculation Dose (" * TCID[50] * ")")) +
  ylab("Peak Time (dpi)") +
  scale_x_continuous(trans = "log10",
                     breaks = c(1e+0, 1e+1, 1e+2, 1e+3, 1e+4, 1e+5, 1e+6, 1e+7),
                     labels = trans_format("log10", math_format(10^.x))) +
  scale_y_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7)) +
  coord_cartesian(xlim = c(1e+2, 1e+7), ylim = c(0, 5))

plot(pXb)

## Fig 4C: effect of 100-fold dose change --------------------------------------

# violin plot of 100-fold ratios
pXc <- fit %>%
  spread_draws(diff_logN7_logN5, diff_logtp7_logtp5) %>%
  pivot_longer(cols = c("diff_logN7_logN5", "diff_logtp7_logtp5"),
               names_to = "ratio",
               values_to = "val") %>%
  ggplot(aes(x = ratio, y = exp(val), fill = ratio)) +
  stat_eye(.width = c(0.66, 0.95), alpha = 0.8) +
  geom_hline(yintercept = 1e+0, linetype = "dashed", alpha = 0.6, linewidth = 1) +
  scale_x_discrete(labels = c("diff_logN7_logN5" = "Initial Viral Load\nFold Change",
                              "diff_logtp7_logtp5" = "Peak Time\nFold Change")) +
  scale_y_continuous(trans = "log10",
                     breaks = c(1e-1, 1e+0, 1e+1, 1e+2, 1e+3, 1e+4),
                     labels = c("0.1", "1", "10", "100", "1000", "10000")) +
  scale_fill_manual(values = c("#488B67", "#725AA6")) +
  xlab(element_blank()) +
  ylab("Effect of 100-Fold Dose Change") +
  theme(legend.position = "none")

plot(pXc)

## Fig 4D: dose-effect on trajectories -----------------------------------------

# extract parameters
draws <- fit %>%
  spread_draws(hit_prob,
               alpha_t, mean_log_tp5,
               alpha_N, mean_log_N5,
               mean_log_t_double, mean_log_t_half)

# sample a subset of extracted parameters
n_samples <- 1000
sampled_draws <- draws %>%
  sample_n(size = n_samples, replace = FALSE) %>%
  mutate(mean_g = log(2) / exp(mean_log_t_double),
         mean_d = log(2) / exp(mean_log_t_half))

# vectors to store successful trajectories
inf_times <- c()
inf_viral_loads <- c()
inf_doses <- c()
inf_ids <- c() # to keep track of individual rows (draws)

# vectors to store failed infections
failed_doses <- c()

# times to compute trajectories
n_time_points <- 1000
t_seq <- seq(from = 0, to = 35, length.out = n_time_points)

# compute and save trajectories for each draw
for (i in 1:nrow(sampled_draws)) {
  
  # sample a dose between 10^3 and 10^7 TCID50
  dose <- 10^runif(1, min=3, max=7)
  
  # extract parameters from draw
  g <- sampled_draws[i,]$mean_g
  d <- sampled_draws[i,]$mean_d
  hit_prob <- sampled_draws[i,]$hit_prob
  alpha_t <- sampled_draws[i,]$alpha_t
  alpha_N <- sampled_draws[i,]$alpha_N
  log_tp5 <- sampled_draws[i,]$mean_log_tp5
  log_N5 <- sampled_draws[i,]$mean_log_N5
  
  # simulate infection, chance that monkey is not infected
  inf_prob <- DR_prob(hit_prob, dose)
  # keep track of failed infections
  if (runif(1) > inf_prob) {
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
  inf_times <- c(inf_times, t_seq)
  inf_viral_loads <- c(inf_viral_loads, viral_loads)
  inf_doses <- c(inf_doses, rep(dose, n_time_points))
  inf_ids <- c(inf_ids, rep(i, n_time_points))

}


# data frame for plotting successful trajectories
df_inf <- data.frame(time = inf_times,
                     viral_load = inf_viral_loads,
                     TCID_dose = inf_doses,
                     id = inf_ids) %>%
  filter(viral_load >= 1) # remove points below 10^0

# create trajectories for failed infections
# lines that will be plotted below zero to indicate failed infections
uninf_dose <- c()
uninf_time <- c()
uninf_viral_load <- c()
uninf_id <- c()
for (i in 1:length(failed_doses)) {
  uninf_dose <- c(uninf_dose, rep(failed_doses[i], 2))
  uninf_time <- c(uninf_time, c(0, 30))
  uninf_viral_load <- c(uninf_viral_load, runif(2, min = 10^-0.8, max = 10^-0.4))
  uninf_id <- c(uninf_id, rep(i, 2))
}

# data frame for plotting failed infections
df_uninf <- data.frame(TCID_dose = uninf_dose,
                       time = uninf_time,
                       viral_load = uninf_viral_load,
                       id = uninf_id)

# add shading, zero line
pXd <- ggplot() +
  geom_ribbon(data = data.frame(x = seq(-5, 40, 1), y = 1, lower = 1e-2, upper = 1e+0),
              aes(ymin = lower, ymax = upper, x = x, y = y),
              alpha = 0.2) +
  geom_hline(yintercept = 1e+0)

# add predicted successful infection
pXd <- pXd + geom_line(data = df_inf,
                       aes(x = time, y = viral_load, group = id, color = log10(TCID_dose)),
                       alpha = 0.1, linewidth = 1)

# add predicted failed infections
pXd <- pXd + geom_line(data = df_uninf,
                       aes(x = time, y = viral_load, group = id, color = log10(TCID_dose)),
                       alpha = 0.1, linewidth = 1)

# final plot settings
pXd <- pXd +
  scale_x_continuous(breaks = seq(0, 40, 5)) +
  scale_y_continuous(trans = "log10",
                     breaks = c(1e+0, 1e+1, 1e+2, 1e+3, 1e+4, 1e+5, 1e+6, 1e+7),
                     labels = trans_format("log10", math_format(10^.x))) +
  coord_cartesian(xlim = c(0, 30), ylim = c(10^-0.8, 10^6)) +
  xlab("Day Post Infection") +
  ylab("Viral Load (RNA Copies/Sample)") +
  scale_color_gradient2(
    low = "#009DDC",
    mid = "#F7CB15",
    high = "#DD1C1A",
    midpoint = 5,
    name = "log10\ndose"
  ) +
  theme(legend.position = c(0.900, 0.730),
        legend.background = element_rect(color = 'black', linewidth = 0.2),
        legend.title = element_text(size = 7), 
        legend.text = element_text(size = 7))

plot(pXd)

## Fig 4: combined dose-dependence results plot --------------------------------

fig <- ((pXa | pXb) / (pXc | pXd)) + plot_annotation(tag_levels = "A")

## Save plot -------------------------------------------------------------------

ggsave("outputs/figures/Fig4.tiff", fig, width = 8, height = 8)
