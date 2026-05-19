rm(list = ls())

## Install and load necessary packages -----------------------------------------

req_pkgs <- c("dplyr", "tidybayes", "ggplot2", "scales", "patchwork")
new_pkgs <- req_pkgs[!(req_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) {
  install.packages(new_pkgs, dependences = TRUE)
}
invisible(lapply(req_pkgs, library, character.only = TRUE))

## Load data and model fit -----------------------------------------------------

data <- read.csv("data/clean-data.csv")
data$individual <- as.factor(data$individual)
data$animal_species <- as.factor(data$animal_species)
data$study_name <- as.factor(data$study_name)
data$viral_isolate <- as.factor(data$viral_isolate)
data$data_type <- as.factor(data$data_type)

fit_original <- readRDS("outputs/fits/fit-all-data.RDS")
fit_original <- fit_original %>%
  recover_types(data)

fit_alphaN1 <- readRDS("outputs/fits/fit-sensitivity-alphaN-2x.RDS")
fit_alphaN1 <- fit_alphaN1 %>%
  recover_types(data)

fit_alphaN2 <- readRDS("outputs/fits/fit-sensitivity-alphaN-3x.RDS")
fit_alphaN2 <- fit_alphaN2 %>%
  recover_types(data)

fit_alphat1 <- readRDS("outputs/fits/fit-sensitivity-alphat-2x.RDS")
fit_alphat1 <- fit_alphat1 %>%
  recover_types(data)

fit_alphat2 <- readRDS("outputs/fits/fit-sensitivity-alphat-3x.RDS")
fit_alphat2 <- fit_alphat2 %>%
  recover_types(data)

fit_original_prior <- readRDS("outputs/fits/fit-prior.RDS")
fit_original_prior <- fit_original_prior %>%
  recover_types(data)

fit_alphaN1_prior <- readRDS("outputs/fits/fit-sensitivity-alphaN-2x-prior.RDS")
fit_alphaN1_prior <- fit_alphaN1_prior %>%
  recover_types(data)

fit_alphaN2_prior <- readRDS("outputs/fits/fit-sensitivity-alphaN-3x-prior.RDS")
fit_alphaN2_prior <- fit_alphaN2_prior %>%
  recover_types(data)

fit_alphat1_prior <- readRDS("outputs/fits/fit-sensitivity-alphat-2x-prior.RDS")
fit_alphat1_prior <- fit_alphat1_prior %>%
  recover_types(data)

fit_alphat2_prior <- readRDS("outputs/fits/fit-sensitivity-alphat-3x-prior.RDS")
fit_alphat2_prior <- fit_alphat2_prior %>%
  recover_types(data)

## Fig S3A: initial viral load dose-dependence sensitivity analysis ------------

theme_set(theme_bw())

# extract posterior samples
df_original <- fit_original %>%
  spread_draws(diff_logN7_logN5)
df_alphaN1 <- fit_alphaN1 %>%
  spread_draws(diff_logN7_logN5)
df_alphaN2 <- fit_alphaN2 %>%
  spread_draws(diff_logN7_logN5)
df_post <- data.frame(diff_N = c(df_original$diff_logN7_logN5,
                                 df_alphaN1$diff_logN7_logN5,
                                 df_alphaN2$diff_logN7_logN5),
                      sd = c(rep("1x SD", nrow(df_original)),
                             rep("2x SD", nrow(df_alphaN1)),
                             rep("3x SD", nrow(df_alphaN2))),
                      side = "right")

# extract prior samples
df_original_prior <- fit_original_prior %>%
  spread_draws(diff_logN7_logN5)
df_alphaN1_prior <- fit_alphaN1_prior %>%
  spread_draws(diff_logN7_logN5)
df_alphaN2_prior <- fit_alphaN2_prior %>%
  spread_draws(diff_logN7_logN5)
df_prior <- data.frame(diff_N = c(df_original_prior$diff_logN7_logN5,
                                  df_alphaN1_prior$diff_logN7_logN5,
                                  df_alphaN2_prior$diff_logN7_logN5),
                       sd = c(rep("1x SD", nrow(df_original_prior)),
                              rep("2x SD", nrow(df_alphaN1_prior)),
                              rep("3x SD", nrow(df_alphaN2_prior))),
                       side = "left")

# plot all samples
df_prior_post <- rbind(df_prior, df_post)

p1 <- df_prior_post %>%
  ggplot(aes(x = as.factor(sd), y = exp(diff_N), side = side, fill = side)) +
  stat_slab(alpha = 0.8, scale = 0.65) +
  scale_fill_manual(values = c("left" = "#B1677A",
                               "right" = "#92374D"),
                    guide = "none") +
  xlab(element_blank()) +
  ylab(bquote(10^7 * ":" * 10^5 ~ "Dose Ratio")) +
  scale_y_continuous(trans = "log10",
                     # limits = c(1e-4, 1e+7),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  coord_cartesian(ylim = c(1e-4, 1e+7)) +
  ggtitle("Initial Viral Load")

plot(p1)

## Fig S3B: peak viral load time dose-dependence sensitivity analysis ----------

theme_set(theme_bw())

# extract posterior samples
df_original <- fit_original %>%
  spread_draws(diff_logtp7_logtp5)
df_alphat1 <- fit_alphat1 %>%
  spread_draws(diff_logtp7_logtp5)
df_alphat2 <- fit_alphat2 %>%
  spread_draws(diff_logtp7_logtp5)
df_post <- data.frame(diff_tp = c(df_original$diff_logtp7_logtp5,
                                  df_alphat1$diff_logtp7_logtp5,
                                  df_alphat2$diff_logtp7_logtp5),
                      sd = c(rep("1x SD", nrow(df_original)),
                             rep("2x SD", nrow(df_alphat1)),
                             rep("3x SD", nrow(df_alphat2))),
                      side = "right")

# extract prior samples
df_original_prior <- fit_original_prior %>%
  spread_draws(diff_logtp7_logtp5)
df_alphat1_prior <- fit_alphat1_prior %>%
  spread_draws(diff_logtp7_logtp5)
df_alphat2_prior <- fit_alphat2_prior %>%
  spread_draws(diff_logtp7_logtp5)
df_prior <- data.frame(diff_tp = c(df_original_prior$diff_logtp7_logtp5,
                                   df_alphat1_prior$diff_logtp7_logtp5,
                                   df_alphat2_prior$diff_logtp7_logtp5),
                       sd = c(rep("1x SD", nrow(df_original_prior)),
                              rep("2x SD", nrow(df_alphat1_prior)),
                              rep("3x SD", nrow(df_alphat2_prior))),
                       side = "left")

# plot all samples
df_prior_post <- rbind(df_prior, df_post)

p2 <- df_prior_post %>%
  ggplot(aes(x = as.factor(sd), y = exp(diff_tp), side = side, fill = side)) +
  stat_slab(alpha = 0.8, scale = 0.75) +
  scale_fill_manual(values = c("left" = "#92BEDF",
                               "right" = "#5AA9E6"),
                    guide = "none") +
  xlab(element_blank()) +
  ylab(bquote(10^7 * ":" * 10^5 ~ "Dose Ratio")) +
  scale_y_continuous(trans = "log10",
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  coord_cartesian(ylim = c(1e-4, 1e+4)) +
  ggtitle("Peak Viral Load Time")

plot(p2)

## Fig S3: combined dose-dependence sensitivity plot ---------------------------

fig <- p1 / p2 + plot_annotation(tag_levels = "A")

## Save plot -------------------------------------------------------------------

ggsave("outputs/figures/S3 Fig.tiff", fig, width = 6, height = 7)
