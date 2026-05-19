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

fit <- readRDS("outputs/fits/fit-all-data.RDS")
fit <- fit %>%
  recover_types(data)

## Fig 5A: ID50 posterior draws ------------------------------------------------

theme_set(theme_bw())

p1 <- fit %>%
  spread_draws(log_ID50) %>%
  mutate(ID50 = exp(log_ID50)) %>%
  ggplot(aes(x = ID50)) +
  geom_density(color = NA,
               fill = "#ED254E",
               alpha = 0.6) +
  stat_pointinterval(.width = c(0.66, 0.95)) +
  scale_x_continuous(trans = "log10",
                     limits = c(1e-1, 1e+7),
                     breaks = c(1e-1, 1e+0, 1e+1, 1e+2, 1e+3, 1e+4, 1e+5, 1e+6, 1e+7),
                     labels = trans_format("log10", math_format(10^.x))) +
  scale_y_continuous(limits = c(0.0, 0.50),
                     breaks = c(0.0, 0.25, 0.50, 0.75, 1.0)) +
  xlab(bquote(ID[50] ~ "(" * TCID[50] * ")")) +
  ylab("Density")

plot(p1)

## Fig 5B: dose-response curve posterior draws ---------------------------------

set.seed(1234)
theme_set(theme_bw())

# sample a subset of posterior hit probabilities
n_samples <- 200
draws <- fit %>%
  spread_draws(hit_prob)
hit_prob_samples <- sample(draws$hit_prob, n_samples)

# set x-axis range
n_doses <- 200
log10_TCID_doses <- seq(-1, 8, length.out = n_doses)
TCID_doses <- 10^log10_TCID_doses

# calculate dose-response curve for each sampled hit probability
final_doses <- c()
final_probs <- c()
final_ids <- c()
for (i in 1:n_samples) {
  probs <- DR_prob(hit_prob_samples[i], TCID_doses)
  final_doses <- c(final_doses, TCID_doses)
  final_probs <- c(final_probs, probs)
  final_ids <- c(final_ids, rep(i, n_doses))
}
dose_response_df <- data.frame(dose = final_doses,
                               prob = final_probs,
                               id = final_ids)

# plot sampled dose-response curves
p2 <- ggplot() +
  geom_line(data = dose_response_df,
            aes(x = dose, y = prob, group = id),
            alpha = 0.1, color = "#ED254E")

# plot median dose-response curve
median_hit_prob <- median(draws$hit_prob)
median_probs <- DR_prob(median_hit_prob, TCID_doses)
df_median <- data.frame(dose = TCID_doses,
                        prob = median_probs)
p2 <- p2 + geom_line(data = df_median,
                   aes(x = dose, y = prob),
                   linewidth = 0.8)

# plot lines indicating doses used in studies
p2 <- p2 + geom_vline(xintercept = c(10^3, 10^5, 10^6, 10^6.9), linetype = "dashed", color = "red")

# settings for axes
p2 <- p2 + scale_x_continuous(trans = "log10",
                            breaks = c(1e-1, 1e+0, 1e+1, 1e+2, 1e+3, 1e+4, 1e+5, 1e+6, 1e+7),
                            labels = trans_format("log10", math_format(10^.x))) +
  coord_cartesian(xlim = c(1e-1, 1e+7)) +
  xlab(bquote("Dose (" * TCID[50] * ")")) + ylab("Probability of Infection") +
  ylab("Probability of Infection")

plot(p2)

## Fig 5: combined dose-response results plot ----------------------------------

fig <- p1 + p2 + plot_annotation(tag_levels = "A")

## Save plot -------------------------------------------------------------------

ggsave("outputs/figures/Fig5.tiff", fig, width = 7, height = 3)
