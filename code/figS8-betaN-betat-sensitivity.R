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

## Compare dose-dependence intercept parameter estimates -----------------------

theme_set(theme_bw())

get_betaN_ests <- function(fit) {
  fit %>%
    spread_draws(log_N5[study_name]) %>%
    select(study_name, log_N5)
}

get_betat_ests <- function(fit) {
  fit %>%
    spread_draws(log_tp5[study_name]) %>%
    select(study_name, log_tp5)
}

# extract betaN estimates
betaN_all <- get_betaN_ests(fit_all)
betaN_oral <- get_betaN_ests(fit_oral)
betaN_nasal <- get_betaN_ests(fit_nasal)
betaN_CM <- get_betaN_ests(fit_CM)
betaN_NN <- get_betaN_ests(fit_NN)

# extract betat estimates
betat_all <- get_betat_ests(fit_all)
betat_oral <- get_betat_ests(fit_oral)
betat_nasal <- get_betat_ests(fit_nasal)
betat_CM <- get_betat_ests(fit_CM)
betat_NN <- get_betat_ests(fit_NN)

# create data frame for plotting betaN
df_betaN <- rbind(betaN_all, betaN_oral, betaN_nasal, betaN_CM, betaN_NN)
df_betaN$fit <- c(rep("All Data", nrow(betaN_all)),
                  rep("Only Oral", nrow(betaN_oral)),
                  rep("Only Nasal", nrow(betaN_nasal)),
                  rep("Only Cynomolgus", nrow(betaN_CM)),
                  rep("Only Nagata", nrow(betaN_NN)))
df_betaN$fit <- factor(df_betaN$fit, levels = c("All Data", "Only Oral", "Only Nasal", "Only Cynomolgus", "Only Nagata"))

# create data frame for plotting betat
df_betat <- rbind(betat_all, betat_oral, betat_nasal, betat_CM, betat_NN)
df_betat$fit <- c(rep("All Data", nrow(betat_all)),
                  rep("Only Oral", nrow(betat_oral)),
                  rep("Only Nasal", nrow(betat_nasal)),
                  rep("Only Cynomolgus", nrow(betat_CM)),
                  rep("Only Nagata", nrow(betat_NN)))
df_betat$fit <- factor(df_betat$fit, levels = c("All Data", "Only Oral", "Only Nasal", "Only Cynomolgus", "Only Nagata"))

# plot betaN
pbetaN <- df_betaN %>%
  ggplot(aes(x = study_name, y = exp(log_N5), color = study_name)) +
  stat_pointinterval(.width = c(0.66, 0.95)) +
  xlab(element_blank()) +
  ylab(bquote("exp(" * beta[N] * ")")) +
  scale_y_continuous(trans = "log10",
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  coord_cartesian(ylim = c(1e+1, 10^4.5)) +
  theme(legend.position = "none") +
  facet_grid(~ fit, scales = "free_x")

plot(pbetaN)

# plot betat
pbetat <- df_betat %>%
  ggplot(aes(x = study_name, y = exp(log_tp5), color = study_name)) +
  stat_pointinterval(.width = c(0.66, 0.95)) +
  xlab(element_blank()) +
  ylab(bquote("exp(" * beta[t] * ")")) +
  scale_y_continuous(trans = "log10",
                     breaks = c(1e-1, 1e+0, 1e+1),
                     labels = trans_format("log10", math_format(10^.x))) +
  coord_cartesian(ylim = c(10^-1, 10^1)) +
  theme(legend.position = "none") +
  facet_grid(~ fit, scales = "free_x")

plot(pbetat)

## Combined plot ---------------------------------------------------------------

fig <- (pbetaN / pbetat) + plot_annotation(tag_levels = "A")

## Save plot -------------------------------------------------------------------

ggsave("outputs/figures/S8 Fig.tiff", fig, width = 12, height = 6)
