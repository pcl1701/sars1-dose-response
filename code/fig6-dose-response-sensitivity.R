rm(list = ls())

## Install and load necessary packages -----------------------------------------

req_pkgs <- c("dplyr", "tidybayes", "ggplot2", "scales", "patchwork")
new_pkgs <- req_pkgs[!(req_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) {
  install.packages(new_pkgs, dependences = TRUE)
}
invisible(lapply(req_pkgs, library, character.only = TRUE))

## Load data and model fit -----------------------------------------------------

data_all <- read.csv("data/clean-data.csv")
data_all$individual <- as.factor(data_all$individual)
data_all$animal_species <- as.factor(data_all$animal_species)
data_all$study_name <- as.factor(data_all$study_name)
data_all$viral_isolate <- as.factor(data_all$viral_isolate)
data_all$data_type <- as.factor(data_all$data_type)

data_nasal <- read.csv("data/clean-data.csv")
data_nasal <- data_nasal %>% filter(data_type == "nasal swab")
data_nasal$individual <- as.factor(data_nasal$individual)
data_nasal$animal_species <- as.factor(data_nasal$animal_species)
data_nasal$study_name <- as.factor(data_nasal$study_name)
data_nasal$viral_isolate <- as.factor(data_nasal$viral_isolate)
data_nasal$data_type <- as.factor(data_nasal$data_type)

data_oral <- read.csv("data/clean-data.csv")
data_oral <- data_oral %>% filter(data_type != "nasal swab")
data_oral$individual <- as.factor(data_oral$individual)
data_oral$animal_species <- as.factor(data_oral$animal_species)
data_oral$study_name <- as.factor(data_oral$study_name)
data_oral$viral_isolate <- as.factor(data_oral$viral_isolate)
data_oral$data_type <- as.factor(data_oral$data_type)

data_RM <- read.csv("data/clean-data.csv")
data_RM <- data_RM %>% filter(animal_species == "RM")
data_RM$individual <- as.factor(data_RM$individual)
data_RM$animal_species <- as.factor(data_RM$animal_species)
data_RM$study_name <- as.factor(data_RM$study_name)
data_RM$viral_isolate <- as.factor(data_RM$viral_isolate)
data_RM$data_type <- as.factor(data_RM$data_type)

data_CM <- read.csv("data/clean-data.csv")
data_CM <- data_CM %>% filter(animal_species == "CM")
data_CM$individual <- as.factor(data_CM$individual)
data_CM$animal_species <- as.factor(data_CM$animal_species)
data_CM$study_name <- as.factor(data_CM$study_name)
data_CM$viral_isolate <- as.factor(data_CM$viral_isolate)
data_CM$data_type <- as.factor(data_CM$data_type)

data_NN <- read.csv("data/clean-data.csv")
data_NN <- data_NN %>% filter(study_name == "NN")
data_NN$individual <- as.factor(data_NN$individual)
data_NN$animal_species <- as.factor(data_NN$animal_species)
data_NN$study_name <- as.factor(data_NN$study_name)
data_NN$viral_isolate <- as.factor(data_NN$viral_isolate)
data_NN$data_type <- as.factor(data_NN$data_type)

fit_all <- readRDS("outputs/fits/fit-all-data.RDS")
fit_all <- fit_all %>%
  recover_types(data_all)

fit_nasal <- readRDS("outputs/fits/fit-nasal-data.RDS")
fit_nasal <- fit_nasal %>%
  recover_types(data_nasal)

fit_oral <- readRDS("outputs/fits/fit-oral-data.RDS")
fit_oral <- fit_oral %>%
  recover_types(data_oral)

fit_CM <- readRDS("outputs/fits/fit-cynomolgus-data.RDS")
fit_CM <- fit_CM %>%
  recover_types(data_CM)

fit_NN <- readRDS("outputs/fits/fit-nagata-et-al-data.RDS")
fit_NN <- fit_NN %>%
  recover_types(data_NN)

## Fig 6A: ID50 posterior by swab type -----------------------------------------

theme_set(theme_bw())

# extract ID50 estimates
ID50s_all <- fit_all %>%
  spread_draws(log_ID50) %>%
  mutate(ID50 = exp(log_ID50)) %>%
  pull(ID50)

ID50s_nasal <- fit_nasal %>%
  spread_draws(log_ID50) %>%
  mutate(ID50 = exp(log_ID50)) %>%
  pull(ID50)

ID50s_oral <- fit_oral %>%
  spread_draws(log_ID50) %>%
  mutate(ID50 = exp(log_ID50)) %>%
  pull(ID50)

# prepare plotting data frames
df_ID50_all <- data.frame(ID50 = ID50s_all,
                          swab = "All Data")
df_ID50_oral_nasal <- data.frame(ID50 = c(ID50s_oral, ID50s_nasal),
                           swab = c(rep("Oral Swabs", length(ID50s_oral)),
                                    rep("Nasal Swabs", length(ID50s_nasal))))
df_ID50_oral_nasal$swab <- factor(df_ID50_oral_nasal$swab, levels = c("Oral Swabs", "Nasal Swabs"))

# plot results

# oral swab doses
df_dose_oral <- data_oral %>%
  select(individual, TCID_dose) %>%
  distinct()
df_dose_oral <- as.data.frame(table(df_dose_oral$TCID_dose))
colnames(df_dose_oral)[1] <- "dose"
df_dose_oral$y <- 0.65
df_dose_oral$dose <- as.numeric(as.character(df_dose_oral$dose))

# nasal swab doses
df_dose_nasal <- data_nasal %>%
  select(individual, TCID_dose) %>%
  distinct()
df_dose_nasal <- as.data.frame(table(df_dose_nasal$TCID_dose))
colnames(df_dose_nasal)[1] <- "dose"
df_dose_nasal$y <- 0.55
df_dose_nasal$dose <- as.numeric(as.character(df_dose_nasal$dose))

p1 <- ggplot() +
  geom_density(data = df_ID50_oral_nasal, aes(x = ID50, fill = swab),
               alpha = 0.7, color = NA) +
  geom_density(data = df_ID50_all, aes(x = ID50),
               alpha = 0.7, color = "black", fill = NA, linetype = "dashed") +
  geom_point(data = df_dose_oral, aes(x = dose, y = y),
             alpha = 0.7, color = "#F46036", size = 5) +
  geom_text(data = df_dose_oral, aes(x = dose, y = y, label = Freq),
            alpha = 1, size = 3) +
  geom_point(data = df_dose_nasal, aes(x = dose, y = y),
             alpha = 0.7, color = "#75B9BE", size = 5) +
  geom_text(data = df_dose_nasal, aes(x = dose, y = y, label = Freq),
            alpha = 1, size = 3)

# final plot settings
p1 <- p1 +
  scale_x_continuous(trans = "log10",
                     limits = c(1e-1, 1e+7),
                     breaks = c(1e-1, 1e+0, 1e+1, 1e+2, 1e+3, 1e+4, 1e+5, 1e+6, 1e+7),
                     labels = trans_format("log10", math_format(10^.x))) +
  scale_y_continuous(limits = c(0.0, 0.75),
                     breaks = c(0.0, 0.25, 0.50, 0.75, 1.00)) +
  xlab(bquote(ID[50] ~ "(" * TCID[50] * ")")) +
  ylab("Density") +
  scale_fill_manual(name = element_blank(),
                    values = c("Oral Swabs" = "#F46036",
                               "Nasal Swabs" = "#75B9BE")) +
  theme(legend.position = "bottom")

plot(p1)

## Fig 6B: ID50 posterior by CM and NN -----------------------------------------

theme_set(theme_bw())

# extract ID50 estimates
ID50s_all <- fit_all %>%
  spread_draws(log_ID50) %>%
  mutate(ID50 = exp(log_ID50)) %>%
  pull(ID50)

ID50s_CM <- fit_CM %>%
  spread_draws(log_ID50) %>%
  mutate(ID50 = exp(log_ID50)) %>%
  pull(ID50)

ID50s_NN <- fit_NN %>%
  spread_draws(log_ID50) %>%
  mutate(ID50 = exp(log_ID50)) %>%
  pull(ID50)

# prepare plotting data frames
df_ID50_all <- data.frame(ID50 = ID50s_all,
                          study = "All Data")
df_ID50_CM_NN <- data.frame(ID50 = c(ID50s_CM, ID50s_NN),
                            study = c(rep("Lawler et al. and Nagata et al.", length(ID50s_CM)),
                                      rep("Only Nagata et al.", length(ID50s_NN))))
df_ID50_CM_NN$study <- factor(df_ID50_CM_NN$study,
                              levels = c("Lawler et al. and Nagata et al.", "Only Nagata et al."))

# plot results

# CM doses
df_dose_CM <- data_CM %>%
  select(individual, TCID_dose) %>%
  distinct()
df_dose_CM <- as.data.frame(table(df_dose_CM$TCID_dose))
colnames(df_dose_CM)[1] <- "dose"
df_dose_CM$y <- 0.65
df_dose_CM$dose <- as.numeric(as.character(df_dose_CM$dose))

# NN doses
df_dose_NN <- data_NN %>%
  select(individual, TCID_dose) %>%
  distinct()
df_dose_NN <- as.data.frame(table(df_dose_NN$TCID_dose))
colnames(df_dose_NN)[1] <- "dose"
df_dose_NN$y <- 0.55
df_dose_NN$dose <- as.numeric(as.character(df_dose_NN$dose))

p2 <- ggplot() +
  geom_density(data = df_ID50_CM_NN, aes(x = ID50, fill = study),
               alpha = 0.7, color = NA) +
  geom_density(data = df_ID50_all, aes(x = ID50),
               alpha = 0.7, color = "black", fill = NA, linetype = "dashed") +
  geom_point(data = df_dose_CM, aes(x = dose, y = y),
             alpha = 0.7, color = "#3A4F41", size = 5) +
  geom_text(data = df_dose_CM, aes(x = dose, y = y, label = Freq),
            alpha = 1, size = 3) +
  geom_point(data = df_dose_NN, aes(x = dose, y = y),
             alpha = 0.7, color = "#B9314F", size = 5) +
  geom_text(data = df_dose_NN, aes(x = dose, y = y, label = Freq),
            alpha = 1, size = 3)

# final plot settings
p2 <- p2 +
  scale_x_continuous(trans = "log10",
                     limits = c(1e-1, 1e+7),
                     breaks = c(1e-1, 1e+0, 1e+1, 1e+2, 1e+3, 1e+4, 1e+5, 1e+6, 1e+7),
                     labels = trans_format("log10", math_format(10^.x))) +
  scale_y_continuous(limits = c(0.0, 0.75),
                     breaks = c(0.0, 0.25, 0.50, 0.75, 1.00)) +
  xlab(bquote(ID[50] ~ "(" * TCID[50] * ")")) +
  ylab("Density") +
  scale_fill_manual(name = element_blank(),
                    values = c("Lawler et al. and Nagata et al." = "#3A4F41",
                               "Only Nagata et al." = "#B9314F")) +
  theme(legend.position = "bottom")

plot(p2)

## Fig 6: combined ID50 plot ---------------------------------------------------

fig <- p1 + p2 + plot_annotation(tag_levels = "A")

## Save plot -------------------------------------------------------------------

ggsave("outputs/figures/Fig6.tiff", fig, width = 8, height = 4)
