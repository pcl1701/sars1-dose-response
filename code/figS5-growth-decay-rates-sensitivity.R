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

## Load data -------------------------------------------------------------------

# all data
data_all <- read.csv("data/clean-data.csv")
data_all$individual <- as.factor(data_all$individual)
data_all$animal_species <- as.factor(data_all$animal_species)
data_all$study_name <- as.factor(data_all$study_name)
data_all$viral_isolate <- as.factor(data_all$viral_isolate)
data_all$data_type <- as.factor(data_all$data_type)

# only nasal swab data
data_nasal <- read.csv("data/clean-data.csv")
data_nasal <- data_nasal %>% filter(data_type == "nasal swab")
data_nasal$individual <- as.factor(data_nasal$individual)
data_nasal$animal_species <- as.factor(data_nasal$animal_species)
data_nasal$study_name <- as.factor(data_nasal$study_name)
data_nasal$viral_isolate <- as.factor(data_nasal$viral_isolate)
data_nasal$data_type <- as.factor(data_nasal$data_type)

# only oral swab data
data_oral <- read.csv("data/clean-data.csv")
data_oral <- data_oral %>% filter(data_type != "nasal swab")
data_oral$individual <- as.factor(data_oral$individual)
data_oral$animal_species <- as.factor(data_oral$animal_species)
data_oral$study_name <- as.factor(data_oral$study_name)
data_oral$viral_isolate <- as.factor(data_oral$viral_isolate)
data_oral$data_type <- as.factor(data_oral$data_type)

# only cynomolgus macaque data
data_CM <- read.csv("data/clean-data.csv")
data_CM <- data_CM %>% filter(animal_species == "CM")
data_CM$individual <- as.factor(data_CM$individual)
data_CM$animal_species <- as.factor(data_CM$animal_species)
data_CM$study_name <- as.factor(data_CM$study_name)
data_CM$viral_isolate <- as.factor(data_CM$viral_isolate)
data_CM$data_type <- as.factor(data_CM$data_type)

# only Nagata et al. data
data_NN <- read.csv("data/clean-data.csv")
data_NN <- data_NN %>% filter(study_name == "NN")
data_NN$individual <- as.factor(data_NN$individual)
data_NN$animal_species <- as.factor(data_NN$animal_species)
data_NN$study_name <- as.factor(data_NN$study_name)
data_NN$viral_isolate <- as.factor(data_NN$viral_isolate)
data_NN$data_type <- as.factor(data_NN$data_type)

## Load model fits -------------------------------------------------------------

# all data fit
fit_all <- readRDS("outputs/fits/fit-all-data.RDS")
fit_all <- fit_all %>%
  recover_types(data_all)

# only nasal data fit
fit_nasal <- readRDS("outputs/fits/fit-nasal-data.RDS")
fit_nasal <- fit_nasal %>%
  recover_types(data_nasal)

# only oral data fit
fit_oral <- readRDS("outputs/fits/fit-oral-data.RDS")
fit_oral <- fit_oral %>%
  recover_types(data_oral)

# only cynomolgus macaque fit
fit_CM <- readRDS("outputs/fits/fit-cynomolgus-data.RDS")
fit_CM <- fit_CM %>%
  recover_types(data_CM)

# only Nagata et al. fit
fit_NN <- readRDS("outputs/fits/fit-nagata-et-al-data.RDS")
fit_NN <- fit_NN %>%
  recover_types(data_NN)

## compare growth and decay rates ----------------------------------------------

theme_set(theme_bw())

get_g_ests <- function(fit) {
  fit %>%
    spread_draws(g[study_name]) %>%
    select(study_name, g)
}

get_d_ests <- function(fit) {
  fit %>%
    spread_draws(d[study_name]) %>%
    select(study_name, d)
}

get_t_double_ests <- function(fit) {
  fit %>%
    spread_draws(log_t_double[study_name]) %>%
    mutate(t_double = exp(log_t_double)) %>%
    select(study_name, t_double)
}

get_t_half_ests <- function(fit) {
  fit %>%
    spread_draws(log_t_half[study_name]) %>%
    mutate(t_half = exp(log_t_half)) %>%
    select(study_name, t_half)
}

# extract growth rate estimates
g_all <- get_g_ests(fit_all)
g_oral <- get_g_ests(fit_oral)
g_nasal <- get_g_ests(fit_nasal)
g_CM <- get_g_ests(fit_CM)
g_NN <- get_g_ests(fit_NN)

# extract decay rate estimates
d_all <- get_d_ests(fit_all)
d_oral <- get_d_ests(fit_oral)
d_nasal <- get_d_ests(fit_nasal)
d_CM <- get_d_ests(fit_CM)
d_NN <- get_d_ests(fit_NN)

# extract doubling time estimates
t_double_all <- get_t_double_ests(fit_all)
t_double_oral <- get_t_double_ests(fit_oral)
t_double_nasal <- get_t_double_ests(fit_nasal)
t_double_CM <- get_t_double_ests(fit_CM)
t_double_NN <- get_t_double_ests(fit_NN)

# extract halving time estimates
t_half_all <- get_t_half_ests(fit_all)
t_half_oral <- get_t_half_ests(fit_oral)
t_half_nasal <- get_t_half_ests(fit_nasal)
t_half_CM <- get_t_half_ests(fit_CM)
t_half_NN <- get_t_half_ests(fit_NN)

# create data frame for plotting growth rates
df_g <- rbind(g_all, g_oral, g_nasal, g_CM, g_NN)
df_g$fit <- c(rep("All Data", nrow(g_all)),
              rep("Only Oral", nrow(g_oral)),
              rep("Only Nasal", nrow(g_nasal)),
              rep("Only Cynomolgus", nrow(g_CM)),
              rep("Only Nagata", nrow(g_NN)))
df_g$fit <- factor(df_g$fit, levels = c("All Data", "Only Oral", "Only Nasal", "Only Cynomolgus", "Only Nagata"))

# create data frame for plotting decay rates
df_d <- rbind(d_all, d_oral, d_nasal, d_CM, d_NN)
df_d$fit <- c(rep("All Data", nrow(d_all)),
              rep("Only Oral", nrow(d_oral)),
              rep("Only Nasal", nrow(d_nasal)),
              rep("Only Cynomolgus", nrow(d_CM)),
              rep("Only Nagata", nrow(d_NN)))
df_d$fit <- factor(df_d$fit, levels = c("All Data", "Only Oral", "Only Nasal", "Only Cynomolgus", "Only Nagata"))

# create data frame for plotting doubling times
df_t_double <- rbind(t_double_all, t_double_oral, t_double_nasal, t_double_CM, t_double_NN)
df_t_double$fit <- c(rep("All Data", nrow(t_double_all)),
                     rep("Only Oral", nrow(t_double_oral)),
                     rep("Only Nasal", nrow(t_double_nasal)),
                     rep("Only Cynomolgus", nrow(t_double_CM)),
                     rep("Only Nagata", nrow(t_double_NN)))
df_t_double$fit <- factor(df_t_double$fit, levels = c("All Data", "Only Oral", "Only Nasal", "Only Cynomolgus", "Only Nagata"))

# create data frame for plotting halving times
df_t_half <- rbind(t_half_all, t_half_oral, t_half_nasal, t_half_CM, t_half_NN)
df_t_half$fit <- c(rep("All Data", nrow(t_half_all)),
                   rep("Only Oral", nrow(t_half_oral)),
                   rep("Only Nasal", nrow(t_half_nasal)),
                   rep("Only Cynomolgus", nrow(t_half_CM)),
                   rep("Only Nagata", nrow(t_half_NN)))
df_t_half$fit <- factor(df_t_half$fit, levels = c("All Data", "Only Oral", "Only Nasal", "Only Cynomolgus", "Only Nagata"))

# plot growth rates
pg <- df_g %>%
  ggplot(aes(x = study_name, y = g, color = study_name)) +
  stat_pointinterval(.width = c(0.66, 0.95)) +
  xlab(element_blank()) +
  ylab("Growth Rate") +
  # scale_y_continuous(breaks = seq(0, 10, 2)) +
  # coord_cartesian(ylim = c(0, 6)) +
  theme(legend.position = "none") +
  facet_grid(~ fit, scales = "free_x")

plot(pg)

# plot decay rates
pd <- df_d %>%
  ggplot(aes(x = study_name, y = d, color = study_name)) +
  stat_pointinterval(.width = c(0.66, 0.95)) +
  xlab(element_blank()) +
  ylab("Decay Rate") +
  # scale_y_continuous(breaks = seq(0, 2, 0.2)) +
  # coord_cartesian(ylim = c(0, 2)) +
  theme(legend.position = "none") +
  facet_grid(~ fit, scales = "free_x")

plot(pd)

# plot doubling times
pt_double <- df_t_double %>%
  ggplot(aes(x = study_name, y = t_double, color = study_name)) +
  stat_pointinterval(.width = c(0.66, 0.95)) +
  xlab(element_blank()) +
  ylab("Doubling Time (Days)") +
  scale_y_continuous(trans = "log10", breaks = c(1e-1, 1e+0, 1e+1)) +
  coord_cartesian(ylim = c(10^-1.5, 10^1.5)) +
  theme(legend.position = "none") +
  facet_grid(~ fit, scales = "free_x")

plot(pt_double)

# plot halving times
pt_half <- df_t_half %>%
  ggplot(aes(x = study_name, y = t_half, color = study_name)) +
  stat_pointinterval(.width = c(0.66, 0.95)) +
  xlab(element_blank()) +
  ylab("Halving Time (Days)") +
  scale_y_continuous(trans = "log10", breaks = c(1e-1, 1e+0, 1e+1)) +
  theme(legend.position = "none") +
  facet_grid(~ fit, scales = "free_x")

plot(pt_half)

## Combined plot ---------------------------------------------------------------

fig <- (pg / pd / pt_double / pt_half) + plot_annotation(tag_levels = "A")

## Save plot -------------------------------------------------------------------

ggsave("outputs/figures/S5 Fig.tiff", fig, width = 12, height = 10)
