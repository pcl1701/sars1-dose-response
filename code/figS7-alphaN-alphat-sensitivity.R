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

## Compare dose-dependence slope parameter estimates ---------------------------

theme_set(theme_bw())

get_alphaN_ests <- function(fit) {
  fit %>%
    spread_draws(alpha_N) %>%
    pull(alpha_N)
}

get_alphat_ests <- function(fit) {
  fit %>%
    spread_draws(alpha_t) %>%
    pull(alpha_t)
}

# extract alphaN estimates
alphaN_all <- get_alphaN_ests(fit_all)
alphaN_oral <- get_alphaN_ests(fit_oral)
alphaN_nasal <- get_alphaN_ests(fit_nasal)
alphaN_CM <- get_alphaN_ests(fit_CM)
alphaN_NN <- get_alphaN_ests(fit_NN)

# extract alphat estimates
alphat_all <- get_alphat_ests(fit_all)
alphat_oral <- get_alphat_ests(fit_oral)
alphat_nasal <- get_alphat_ests(fit_nasal)
alphat_CM <- get_alphat_ests(fit_CM)
alphat_NN <- get_alphat_ests(fit_NN)

# create data frame for plotting alphaN
df_alphaN <- data.frame(alphaN = c(alphaN_all, alphaN_oral, alphaN_nasal, alphaN_CM, alphaN_NN),
                        fit = rep(c("All Data", "Only Oral", "Only Nasal", "Only Cynomolgus", "Only Nagata"),
                                  each = length(alphaN_all)))
df_alphaN$fit <- factor(df_alphaN$fit, levels = c("All Data", "Only Oral", "Only Nasal", "Only Cynomolgus", "Only Nagata"))

# create data frame for plotting alphat
df_alphat <- data.frame(alphat = c(alphat_all, alphat_oral, alphat_nasal, alphat_CM, alphat_NN),
                        fit = rep(c("All Data", "Only Oral", "Only Nasal", "Only Cynomolgus", "Only Nagata"),
                                  each = length(alphat_all)))
df_alphat$fit <- factor(df_alphat$fit, levels = c("All Data", "Only Oral", "Only Nasal", "Only Cynomolgus", "Only Nagata"))

# plot alphaN estimates
palphaN <- df_alphaN %>%
  ggplot(aes(x = fit, y = alphaN, color = fit)) +
  stat_pointinterval(.width = c(0.66, 0.95)) +
  xlab(element_blank()) +
  ylab(bquote(alpha[N])) +
  # coord_cartesian(ylim = c(0, 4)) +
  theme(legend.position = "none")

plot(palphaN)

# plot alphat estimates
palphat <- df_alphat %>%
  ggplot(aes(x = fit, y = alphat, color = fit)) +
  stat_pointinterval(.width = c(0.66, 0.95)) +
  xlab(element_blank()) +
  ylab(bquote(alpha[t])) +
  # coord_cartesian(ylim = c(-1, 1)) +
  theme(legend.position = "none")

plot(palphat)

## Combined plot ---------------------------------------------------------------

fig <- (palphaN / palphat) + plot_annotation(tag_levels = "A")

## Save plot -------------------------------------------------------------------

ggsave("outputs/figures/S7 Fig.tiff", fig, width = 6, height = 5)
