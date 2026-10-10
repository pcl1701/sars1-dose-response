rm(list = ls())

## Install and load necessary packages -----------------------------------------

req_pkgs <- c("dplyr", "tidybayes", "ggplot2", "scales", "patchwork")
new_pkgs <- req_pkgs[!(req_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) {
  install.packages(new_pkgs, dependences = TRUE)
}
invisible(lapply(req_pkgs, library, character.only = TRUE))

theme_set(theme_bw())

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

fit_original_prior <- readRDS("outputs/fits/fit-prior.RDS")
fit_original_prior <- fit_original_prior %>%
  recover_types(data)

fit_wider <- readRDS("outputs/fits/fit-sensitivity-pfalse-wider.RDS")
fit_wider <- fit_wider %>%
  recover_types(data)

fit_wider_prior <- readRDS("outputs/fits/fit-sensitivity-pfalse-wider-prior.RDS")
fit_wider_prior <- fit_wider_prior %>%
  recover_types(data)

## ID50 priors and posteriors --------------------------------------------------

# extract ID50 estimates
ID50s_original <- fit_original %>%
  spread_draws(log_ID50) %>%
  mutate(ID50 = exp(log_ID50)) %>%
  pull(ID50)

ID50s_original_prior <- fit_original_prior %>%
  spread_draws(log_ID50) %>%
  mutate(ID50 = exp(log_ID50)) %>%
  pull(ID50)

ID50s_wider <- fit_wider %>%
  spread_draws(log_ID50) %>%
  mutate(ID50 = exp(log_ID50)) %>%
  pull(ID50)

ID50s_wider_prior <- fit_wider_prior %>%
  spread_draws(log_ID50) %>%
  mutate(ID50 = exp(log_ID50)) %>%
  pull(ID50)

# prepare plotting data frames
df_ID50_original <- data.frame(ID50 = c(ID50s_original, ID50s_original_prior),
                               type = c(rep("Posterior", length(ID50s_original)),
                                        rep("Prior", length(ID50s_original_prior))))

df_ID50_wider <- data.frame(ID50 = c(ID50s_wider, ID50s_wider_prior),
                            type = c(rep("Posterior", length(ID50s_wider)),
                                     rep("Prior", length(ID50s_wider_prior))))

p1 <- ggplot() +
  # Posterior
  geom_density(
    data = subset(df_ID50_original, type == "Posterior"),
    aes(x = ID50),
    color = NA,
    fill = "#ED254E",
    alpha = 0.6
  ) +
  # Prior
  geom_density(
    data = subset(df_ID50_original, type == "Prior"),
    aes(x = ID50),
    color = "black",
    fill = NA,
    linetype = "dashed",
    linewidth = 0.8
  ) +
  scale_x_continuous(
    trans = "log10",
    limits = c(1e-1, 1e+7),
    breaks = c(
      1e-1, 1e+0, 1e+1, 1e+2,
      1e+3, 1e+4, 1e+5, 1e+6, 1e+7
    ),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  scale_y_continuous(
    limits = c(0.0, 0.50),
    breaks = c(0.0, 0.25, 0.50, 0.75, 1.0)
  ) +
  xlab(bquote(ID[50] ~ "(" * TCID[50] * ")")) +
  ylab("Density")

p2 <- ggplot() +
  # Posterior
  geom_density(
    data = subset(df_ID50_wider, type == "Posterior"),
    aes(x = ID50),
    color = NA,
    fill = "#ED254E",
    alpha = 0.6
  ) +
  # Prior
  geom_density(
    data = subset(df_ID50_wider, type == "Prior"),
    aes(x = ID50),
    color = "black",
    fill = NA,
    linetype = "dashed",
    linewidth = 0.8
  ) +
  scale_x_continuous(
    trans = "log10",
    limits = c(1e-1, 1e+7),
    breaks = c(
      1e-1, 1e+0, 1e+1, 1e+2,
      1e+3, 1e+4, 1e+5, 1e+6, 1e+7
    ),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  scale_y_continuous(
    limits = c(0.0, 0.50),
    breaks = c(0.0, 0.25, 0.50, 0.75, 1.0)
  ) +
  xlab(bquote(ID[50] ~ "(" * TCID[50] * ")")) +
  ylab("Density")


## pfalse priors and posteriors ------------------------------------------------

# extract pfalse estimates
pfalse_original <- fit_original %>%
  spread_draws(p_false) %>%
  pull(p_false)

pfalse_original_prior <- fit_original_prior %>%
  spread_draws(p_false) %>%
  pull(p_false)

pfalse_wider <- fit_wider %>%
  spread_draws(p_false) %>%
  pull(p_false)

pfalse_wider_prior <- fit_wider_prior %>%
  spread_draws(p_false) %>%
  pull(p_false)

# prepare plotting data frames
df_pfalse_original <- data.frame(pfalse = c(pfalse_original, pfalse_original_prior),
                                 type = c(rep("Posterior", length(pfalse_original)),
                                          rep("Prior", length(pfalse_original_prior))))

df_pfalse_wider <- data.frame(pfalse = c(pfalse_wider, pfalse_wider_prior),
                              type = c(rep("Posterior", length(pfalse_wider)),
                                       rep("Prior", length(pfalse_wider_prior))))

p3 <- ggplot() +
  # Posterior
  geom_density(
    data = subset(df_pfalse_original, type == "Posterior"),
    aes(x = pfalse),
    color = NA,
    fill = "#07BEB8",
    alpha = 0.6
  ) +
  # Prior
  geom_density(
    data = subset(df_pfalse_original, type == "Prior"),
    aes(x = pfalse),
    color = "black",
    fill = NA,
    linetype = "dashed",
    linewidth = 0.8
  ) +
  scale_x_continuous(
    trans = "log10",
    limits = c(1e-5, 1e+0),
    breaks = c(
      1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e+0
    ),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  xlab("p_false") +
  ylab("Density")


p4 <- ggplot() +
  # Posterior
  geom_density(
    data = subset(df_pfalse_wider, type == "Posterior"),
    aes(x = pfalse),
    color = NA,
    fill = "#07BEB8",
    alpha = 0.6
  ) +
  # Prior
  geom_density(
    data = subset(df_pfalse_wider, type == "Prior"),
    aes(x = pfalse),
    color = "black",
    fill = NA,
    linetype = "dashed",
    linewidth = 0.8
  ) +
  scale_x_continuous(
    trans = "log10",
    limits = c(1e-5, 1e+0),
    breaks = c(
      1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e+0
    ),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  xlab("p_false") +
  ylab("Density")

## combined plot ---------------------------------------------------------------

p3 <- p3 +
  ggtitle("Original p_false Prior") +
  theme(
    plot.title
    = element_text(hjust = 0.5)
  )

p4 <- p4 +
  ggtitle("Wider p_false Prior") +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

fig <- (p3 + p4) / (p1 + p2)

## Save plot -------------------------------------------------------------------

ggsave("outputs/figures/S5 Fig.tiff", fig, width = 8, height = 6)
