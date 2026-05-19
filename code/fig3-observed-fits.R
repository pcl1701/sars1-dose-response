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

## Overall settings ------------------------------------------------------------

set.seed(1234)
theme_set(theme_bw())

## Fig 3A: fit to observed data ------------------------------------------------

# extract parameters
draws <- fit %>%
  spread_draws(hit_prob,
               alpha_t, log_tp5[study_name],
               alpha_N, log_N5[study_name],
               g[study_name], d[study_name])

# study names + doses
study_names_doses <- list(BL = c(1e+5),
                          JL = c(8571429),
                          NN = c(1e+3, 1e+6),
                          YC = c(1e+5),
                          LL16 = c(1e+5),
                          LL19 = c(1e+5))

# randomly select draws indices
n_draws <- 100
draw_ids <- sample(1:length(unique(draws$.draw)), size = n_draws, replace = FALSE)

# vectors to store successful trajectories
inf_times <- c()
inf_viral_loads <- c()
inf_study_names <- c()
inf_doses <- c()
inf_ids <- c() # to keep track of individual rows (draws)

# vectors to store failed infections
failed_study_names <- c()
failed_doses <- c()

# keep track of which rows were successful infections for each study
success_ids <- matrix(TRUE, ncol = 7, nrow = n_draws)

# times to compute trajectories
n_time_points <- 1000
t_seq <- seq(from = 0, to = 35, length.out = n_time_points)

col_idx <- 1 # to keep track of successful infections
for (study_name in names(study_names_doses)) {

  # extract dose levels for current study
  cur_doses <- study_names_doses[[study_name]]
  
  # filter draws by study name
  cur_study <- study_name
  cur_draws <- draws %>% filter(study_name == cur_study)
  
  for (dose in cur_doses) {
    
    row_idx <- 1 # to keep track of successful infections
    # compute trajectory for each row (draw)
    for (i in draw_ids) {
      
      # extract parameters from draw
      g <- cur_draws[i,]$g
      d <- cur_draws[i,]$d
      hit_prob <- cur_draws[i,]$hit_prob
      alpha_t <- cur_draws[i,]$alpha_t
      alpha_N <- cur_draws[i,]$alpha_N
      log_tp5 <- cur_draws[i,]$log_tp5
      log_N5 <- cur_draws[i,]$log_N5
      
      # simulate infection, chance that monkey is not infected
      inf_prob <- DR_prob(hit_prob, dose)
      # keep track of failed infections
      if (runif(1) > inf_prob) {
        failed_study_names <- c(failed_study_names, study_name)
        failed_doses <- c(failed_doses, dose)
        success_ids[row_idx,col_idx] <- FALSE
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
      inf_study_names <- c(inf_study_names, rep(study_name, n_time_points))
      inf_doses <- c(inf_doses, rep(dose, n_time_points))
      inf_ids <- c(inf_ids, rep(i, n_time_points))
      
      row_idx <- row_idx + 1 # to keep track of successful infections
    }
    
    col_idx <- col_idx + 1 # to keep track of successful infections
    
  }

}

# data frame for plotting successful trajectories
df_inf <- data.frame(time = inf_times,
                     viral_load = inf_viral_loads,
                     study_name = inf_study_names,
                     TCID_dose = inf_doses,
                     id = inf_ids) %>%
  filter(viral_load >= 1) # remove points below 10^0

# create trajectories for failed infections
# lines that will be plotted below zero to indicate failed infections
uninf_times <- c()
uninf_viral_loads <- c()
uninf_study_names <- c()
uninf_doses <- c()
uninf_ids <- c() # to keep track of individual rows (draws)
for (i in 1:length(failed_study_names)) {
  uninf_times <- c(uninf_times, c(0, 30))
  uninf_viral_loads <- c(uninf_viral_loads, runif(2, min = 10^-0.8, max = 10^-0.4))
  uninf_study_names <- c(uninf_study_names, rep(failed_study_names[i], 2))
  uninf_doses <- c(uninf_doses, rep(failed_doses[i], 2))
  uninf_ids <- c(uninf_ids, rep(i, 2))
}

# data frame for plotting failed infections
df_uninf <- data.frame(time = uninf_times,
                       viral_load = uninf_viral_loads,
                       study_name = uninf_study_names,
                       TCID_dose = uninf_doses,
                       id = uninf_ids)

# data frame for plotting LODs
df_LOD <- data.frame(study_name = c("BL", "JL", "NN", "NN", "YC", "LL16", "LL19"),
                     TCID_dose = c(1e+5, 8571429, 1e+3, 1e+6, 1e+5, 1e+5, 1e+5),
                     val = c(10^1.56, 10^2, 10^3.14, 10^3.14, 10^1, 10^1, 10^1))

# data frame for plotting observed data
df_obs <- data %>%
  filter(log10_viral_load != 100) %>% # remove "+" measurements
  filter(data_type != "nasal swab") %>% # remove nasal swabs
  mutate(below_LOD = case_when(log10_viral_load == -100 ~ 1,
                               TRUE ~ 0)) %>%
  select(individual, TCID_dose, day_post_infection, log10_viral_load, study_name, below_LOD)
# manually spread out negative measurements below the LOD
df_obs <- df_obs %>%
  mutate(log10_viral_load = case_when(individual == "JL_CM-03" & below_LOD == 1 ~ -0.4,
                                      individual == "JL_CM-04" & below_LOD == 1 ~ -0.8,
                                      individual == "NN_CM-01" & below_LOD == 1 ~ -0.4,
                                      individual == "NN_CM-02" & below_LOD == 1 ~ -0.8,
                                      individual == "NN_CM-03" & below_LOD == 1 ~ -0.4,
                                      individual == "NN_CM-04" & below_LOD == 1 ~ -0.8,
                                      individual == "LL_RM-14" & below_LOD == 1 ~ -0.6,
                                      TRUE ~ log10_viral_load))

# set facet values
facet_labels <- c("'Li et al. (2005) -' ~ 10^5 ~ TCID[50]",
                  "'Lawler et al. (2006) -' ~ 10^6.9 ~ TCID[50]",
                  "'Liu et al. (2016) -' ~ 10^5 ~ TCID[50]",
                  "'Liu et al. (2019) -' ~ 10^5 ~ TCID[50]",
                  "'Nagata et al. (2007) -' ~ 10^3 ~ TCID[50]",
                  "'Nagata et al. (2007) -' ~ 10^6 ~ TCID[50]",
                  "'Chen et al. (2008) -' ~ 10^5 ~ TCID[50]")
df_inf <- df_inf %>%
  mutate(facet = case_when(
    (study_name == "BL") ~ facet_labels[1],
    (study_name == "JL") ~ facet_labels[2],
    (study_name == "LL16") ~ facet_labels[3],
    (study_name == "LL19") ~ facet_labels[4],
    (study_name == "NN") & (TCID_dose == 1e+3) ~ facet_labels[5],
    (study_name == "NN") & (TCID_dose == 1e+6) ~ facet_labels[6],
    (study_name == "YC") ~ facet_labels[7]
  ))
df_uninf <- df_uninf %>%
  mutate(facet = case_when(
    (study_name == "BL") ~ facet_labels[1],
    (study_name == "JL") ~ facet_labels[2],
    (study_name == "LL16") ~ facet_labels[3],
    (study_name == "LL19") ~ facet_labels[4],
    (study_name == "NN") & (TCID_dose == 1e+3) ~ facet_labels[5],
    (study_name == "NN") & (TCID_dose == 1e+6) ~ facet_labels[6],
    (study_name == "YC") ~ facet_labels[7]
  ))
df_obs <- df_obs %>%
  mutate(facet = case_when(
    (study_name == "BL") ~ facet_labels[1],
    (study_name == "JL") ~ facet_labels[2],
    (study_name == "LL16") ~ facet_labels[3],
    (study_name == "LL19") ~ facet_labels[4],
    (study_name == "NN") & (TCID_dose == 1e+3) ~ facet_labels[5],
    (study_name == "NN") & (TCID_dose == 1e+6) ~ facet_labels[6],
    (study_name == "YC") ~ facet_labels[7]
  ))
df_LOD <- df_LOD %>%
  mutate(facet = case_when(
    (study_name == "BL") ~ facet_labels[1],
    (study_name == "JL") ~ facet_labels[2],
    (study_name == "LL16") ~ facet_labels[3],
    (study_name == "LL19") ~ facet_labels[4],
    (study_name == "NN") & (TCID_dose == 1e+3) ~ facet_labels[5],
    (study_name == "NN") & (TCID_dose == 1e+6) ~ facet_labels[6],
    (study_name == "YC") ~ facet_labels[7]
  ))
# set facet order
facet_order <- c("'Nagata et al. (2007) -' ~ 10^3 ~ TCID[50]",
                 "'Li et al. (2005) -' ~ 10^5 ~ TCID[50]",
                 "'Chen et al. (2008) -' ~ 10^5 ~ TCID[50]",
                 "'Liu et al. (2016) -' ~ 10^5 ~ TCID[50]",
                 "'Liu et al. (2019) -' ~ 10^5 ~ TCID[50]",
                 "'Nagata et al. (2007) -' ~ 10^6 ~ TCID[50]",
                 "'Lawler et al. (2006) -' ~ 10^6.9 ~ TCID[50]")
df_inf$facet <- factor(df_inf$facet, levels = facet_order)
df_uninf$facet <- factor(df_uninf$facet, levels = facet_order)
df_obs$facet <- factor(df_obs$facet, levels = facet_order)
df_LOD$facet <- factor(df_LOD$facet, levels = facet_order)

df_dose <- data.frame(facet = df_LOD$facet) %>%
  mutate(dose = case_when(facet == "'Li et al. (2005) -' ~ 10^5 ~ TCID[50]" ~ "10^5 ~ TCID[50]",
                          facet == "'Lawler et al. (2006) -' ~ 10^6.9 ~ TCID[50]" ~ "10^6.9 ~ TCID[50]",
                          facet == "'Nagata et al. (2007) -' ~ 10^3 ~ TCID[50]" ~ "10^3 ~ TCID[50]",
                          facet == "'Nagata et al. (2007) -' ~ 10^6 ~ TCID[50]" ~ "10^6 ~ TCID[50]",
                          facet == "'Chen et al. (2008) -' ~ 10^5 ~ TCID[50]" ~ "10^5 ~ TCID[50]",
                          facet == "'Liu et al. (2016) -' ~ 10^5 ~ TCID[50]" ~ "10^5 ~ TCID[50]",
                          facet == "'Liu et al. (2019) -' ~ 10^5 ~ TCID[50]" ~ "10^5 ~ TCID[50]"))

facet_labels <- c("'Nagata et al. (2007) -' ~ 10^3 ~ TCID[50]" = "Nagata et al. (2007)",
                  "'Li et al. (2005) -' ~ 10^5 ~ TCID[50]" = "Li et al. (2005)",
                  "'Chen et al. (2008) -' ~ 10^5 ~ TCID[50]" = "Chen et al. (2008)",
                  "'Liu et al. (2016) -' ~ 10^5 ~ TCID[50]" = "Liu et al. (2016)",
                  "'Liu et al. (2019) -' ~ 10^5 ~ TCID[50]" = "Liu et al. (2019)",
                  "'Nagata et al. (2007) -' ~ 10^6 ~ TCID[50]" = "Nagata et al. (2007)",
                  "'Lawler et al. (2006) -' ~ 10^6.9 ~ TCID[50]" = "Lawler et al. (2006)")

# color settings
obs_color <- "#222E50"
colors <- c("#A41623",
            "#FF5714",
            "#F3A712",
            "#2C6E49",
            "#31AFD4",
            "#454E9E",
            "#9C86FF")

# add shading, zero line, LOD lines
pXa <- ggplot() +
  geom_ribbon(data = data.frame(x = seq(-5, 40, 1), y = 1, lower = 1e-2, upper = 1e+0),
              aes(ymin = lower, ymax = upper, x = x, y = y),
              alpha = 0.2) +
  geom_hline(yintercept = 1e+0) +
  geom_hline(data = df_LOD,
             aes(yintercept = val),
             alpha = 0.5, linetype = "dashed")

# add predicted successful infection
pXa <- pXa + geom_line(data = df_inf,
                       aes(x = time, y = viral_load, group = id, color = facet),
                       alpha = 0.15)

# add predicted failed infections
pXa <- pXa + geom_line(data = df_uninf,
                       aes(x = time, y = viral_load, group = id, color = facet),
                       alpha = 0.15)

# add observed data
pXa <- pXa +
  geom_line(data = df_obs,
            aes(x = day_post_infection, y = 10^log10_viral_load, group = individual),
            alpha = 0.4,
            linewidth = 0.3,
            color = obs_color) +
  geom_point(data = df_obs,
             aes(x = day_post_infection, y = 10^log10_viral_load, shape = as.factor(below_LOD)),
             alpha = 0.4,
             size = 2,
             stroke = 0,
             fill = obs_color) +
  scale_shape_manual(values = c(21, 25))

# add dose labels
pXa <- pXa +
  geom_label(data = df_dose,
             aes(x = Inf, y  = Inf, label = dose),
             hjust = "right", vjust = "top", parse = TRUE,
             fill = NA, label.size = NA)

# final plot settings
pXa <- pXa +
  facet_wrap(~ facet, nrow = 4, labeller = as_labeller(facet_labels)) +
  scale_x_continuous(breaks = seq(0, 40, 5)) +
  scale_y_continuous(trans = "log10",
                     breaks = c(1e+0, 1e+1, 1e+2, 1e+3, 1e+4, 1e+5, 1e+6, 1e+7),
                     labels = trans_format("log10", math_format(10^.x))) +
  scale_color_manual(values = colors) +
  # scale_color_brewer(palette = "Dark2") +
  coord_cartesian(xlim = c(0, 30), ylim = c(10^-0.8, 10^6)) +
  xlab("Day Post Infection") +
  ylab("Viral Load (RNA Copies/Sample)") +
  theme(legend.position = "none")

plot(pXa)

## Fig 3B: initial viral load distributions ------------------------------------

final_studies <- c()
final_doses <- c()
final_init_N <- c()

col_idx <- 1 # to keep track of study+dose combination

for (study_name in names(study_names_doses)) {
  
  # extract dose levels for current study
  cur_doses <- study_names_doses[[study_name]]
  
  # filter draws by study name
  cur_study <- study_name
  cur_draws <- draws %>% filter(study_name == cur_study)
  cur_draws <- cur_draws[draw_ids,]
  
  for (dose in cur_doses) {
    
    # extract draws for successful infections only
    success_draws <- cur_draws[success_ids[,col_idx],]
    
    # compute initial viral loads
    init_N <- success_draws %>%
      mutate(init_N = exp(log_init_RNA(hit_prob, dose, alpha_N, log_N5))) %>%
      pull(init_N)
    
    # save initial viral loads
    final_studies <- c(final_studies, rep(study_name, length(init_N)))
    final_doses <- c(final_doses, rep(dose, length(init_N)))
    final_init_N <- c(final_init_N, init_N)
    
    col_idx <- col_idx + 1
    
  }
  
}

# combine and plot
df_init_N <- data.frame(study_name = final_studies,
                        TCID_dose = final_doses,
                        init_N = final_init_N)

# set study values
study_labels <- c("Li et al. (2005)",
                  "Lawler et al. (2006)",
                  "Nagata et al. (2007)",
                  "Chen et al. (2008)",
                  "Liu et al. (2016)",
                  "Liu et al. (2019)")
df_init_N <- df_init_N %>%
  mutate(study = case_when(
    (study_name == "BL") ~ study_labels[1],
    (study_name == "JL") ~ study_labels[2],
    (study_name == "NN") ~ study_labels[3],
    (study_name == "YC") ~ study_labels[4],
    (study_name == "LL16") ~ study_labels[5],
    (study_name == "LL19") ~ study_labels[6]
  ))
# set study order
study_order <- c("Li et al. (2005)",
                 "Lawler et al. (2006)",
                 "Nagata et al. (2007)",
                 "Chen et al. (2008)",
                 "Liu et al. (2016)",
                 "Liu et al. (2019)")
df_init_N$study <- factor(df_init_N$study, levels = study_order)

# set facet values
facet_labels <- c("10^3 ~ TCID[50]",
                  "10^5 ~ TCID[50]",
                  "10^6 ~ TCID[50]",
                  "10^6.9 ~ TCID[50]")
df_init_N <- df_init_N %>%
  mutate(facet = case_when(
    (study_name == "BL") ~ facet_labels[2],
    (study_name == "JL") ~ facet_labels[4],
    (study_name == "LL16") ~ facet_labels[2],
    (study_name == "LL19") ~ facet_labels[2],
    (study_name == "NN") & (TCID_dose == 1e+3) ~ facet_labels[1],
    (study_name == "NN") & (TCID_dose == 1e+6) ~ facet_labels[3],
    (study_name == "YC") ~ facet_labels[2]
  ))
# set facet order
facet_order <- c("10^3 ~ TCID[50]",
                 "10^5 ~ TCID[50]",
                 "10^6 ~ TCID[50]",
                 "10^6.9 ~ TCID[50]")
df_init_N$facet <- factor(df_init_N$facet, levels = facet_order)

pXb <- ggplot(df_init_N, aes(x = study, y = init_N, color = interaction(study, facet))) +
  stat_pointinterval(.width = c(0.66, 0.95)) +
  scale_y_continuous(trans = "log10",
                     breaks = c(1e-2, 1e-1, 1e+0, 1e+1, 1e+2, 1e+3, 1e+4, 1e+5, 1e+6, 1e+7),
                     labels = trans_format("log10", math_format(10^.x)),
                     limits = c(1e-1, 1e+6)) +
  scale_color_manual(values = colors) +
  # scale_color_brewer(palette = "Dark2") +
  xlab(element_blank()) +
  ylab("Initial Viral Load (RNA Copies/Sample)") +
  facet_grid(.~facet, labeller = label_parsed,
             scales = "free", space = "free_y", switch = "y") +
  theme(legend.position = "none",
        panel.spacing = unit(0, "in"),
        axis.text.x = element_blank())

plot(pXb)

## Fig 3C: peak time distributions ---------------------------------------------

final_studies <- c()
final_doses <- c()
final_peak_t <- c()

col_idx <- 1 # to keep track of study+dose combination

for (study_name in names(study_names_doses)) {
  
  # extract dose levels for current study
  cur_doses <- study_names_doses[[study_name]]
  
  # filter draws by study name
  cur_study <- study_name
  cur_draws <- draws %>% filter(study_name == cur_study)
  cur_draws <- cur_draws[draw_ids,]
  
  for (dose in cur_doses) {
    
    # extract draws for successful infections only
    success_draws <- cur_draws[success_ids[,col_idx],]
    
    # compute peak times
    peak_t <- success_draws %>%
      mutate(peak_t = peak_time(hit_prob, dose, alpha_t, log_tp5)) %>%
      pull(peak_t)
    
    # save peak times
    final_studies <- c(final_studies, rep(study_name, length(peak_t)))
    final_doses <- c(final_doses, rep(dose, length(peak_t)))
    final_peak_t <- c(final_peak_t, peak_t)
    
    col_idx <- col_idx + 1
    
  }
  
}

# combine and plot
df_peak_t <- data.frame(study_name = final_studies,
                        TCID_dose = final_doses,
                        peak_t = final_peak_t)

# set study values
study_labels <- c("Li et al. (2005)",
                  "Lawler et al. (2006)",
                  "Nagata et al. (2007)",
                  "Chen et al. (2008)",
                  "Liu et al. (2016)",
                  "Liu et al. (2019)")
df_peak_t <- df_peak_t %>%
  mutate(study = case_when(
    (study_name == "BL") ~ study_labels[1],
    (study_name == "JL") ~ study_labels[2],
    (study_name == "NN") ~ study_labels[3],
    (study_name == "YC") ~ study_labels[4],
    (study_name == "LL16") ~ study_labels[5],
    (study_name == "LL19") ~ study_labels[6]
  ))
# set study order
study_order <- c("Li et al. (2005)",
                 "Lawler et al. (2006)",
                 "Nagata et al. (2007)",
                 "Chen et al. (2008)",
                 "Liu et al. (2016)",
                 "Liu et al. (2019)")
df_peak_t$study <- factor(df_peak_t$study, levels = study_order)

# set facet values
facet_labels <- c("10^3 ~ TCID[50]",
                  "10^5 ~ TCID[50]",
                  "10^6 ~ TCID[50]",
                  "10^6.9 ~ TCID[50]")
df_peak_t <- df_peak_t %>%
  mutate(facet = case_when(
    (study_name == "BL") ~ facet_labels[2],
    (study_name == "JL") ~ facet_labels[4],
    (study_name == "LL16") ~ facet_labels[2],
    (study_name == "LL19") ~ facet_labels[2],
    (study_name == "NN") & (TCID_dose == 1e+3) ~ facet_labels[1],
    (study_name == "NN") & (TCID_dose == 1e+6) ~ facet_labels[3],
    (study_name == "YC") ~ facet_labels[2]
  ))
# set facet order
facet_order <- c("10^3 ~ TCID[50]",
                 "10^5 ~ TCID[50]",
                 "10^6 ~ TCID[50]",
                 "10^6.9 ~ TCID[50]")
df_peak_t$facet <- factor(df_peak_t$facet, levels = facet_order)

pXc <- ggplot(df_peak_t, aes(x = study, y = peak_t, color = interaction(study, facet))) +
  stat_pointinterval(.width = c(0.66, 0.95)) +
  scale_y_continuous(trans = "log10",
                     breaks = c(1e-2, 1e-1, 1e+0, 1e+1, 1e+2, 1e+3, 1e+4, 1e+5, 1e+6, 1e+7),
                     labels = trans_format("log10", math_format(10^.x)),
                     limits = c(1e-1, 1e1)) +
  scale_color_manual(values = colors) +
  # scale_color_brewer(palette = "Dark2") +
  xlab(element_blank()) +
  ylab("Time Peak Viral Load (DPI)") +
  facet_grid(.~facet, labeller = label_parsed,
             scales = "free", space = "free_y", switch = "y") +
  theme(legend.position = "none",
        panel.spacing = unit(0, "in"),
        axis.text.x = element_blank())

plot(pXc)

## Fig 3D: infection duration distributions ------------------------------------

# infection duration (time spent with viral load > 1)
# INPUTS:
# (log_peak_N): log peak viral load
# (peak_t): peak time
# (g): growth rate
# (d): decay rate
inf_dur <- function(log_peak_N, peak_t, g, d) {
  
  # viral load never greater than 1
  if (log_peak_N < 0) {
    return(0)
  }
  
  # time before peak viral load
  # given log peak viral load, peak time, and growth rate
  # equation is: log(N) = g * t + (log_peak_N - g * peak_t)
  # time when log(N) = 0 ==> t = -(log_peak_N - g * peak_t) / g
  
  # time after peak viral load
  # given log peak viral load, peak time, and decay rate
  # equation is: log(N) = -d * t + (log_peak_N + d * peak_t)
  # time when log(N) = 0 ==> t = (log_peak_N + d * peak_t) / d
  
  pre_zero_t <- max(-(log_peak_N - g * peak_t) / g, 0)
  post_zero_t <- (log_peak_N + d * peak_t) / d
  
  return((peak_t - pre_zero_t) + (post_zero_t - peak_t))
}

final_studies <- c()
final_doses <- c()
final_duration <- c()

col_idx <- 1 # to keep track of study+dose combination

for (study_name in names(study_names_doses)) {
  
  # extract dose levels for current study
  cur_doses <- study_names_doses[[study_name]]
  
  # filter draws by study name
  cur_study <- study_name
  cur_draws <- draws %>% filter(study_name == cur_study)
  cur_draws <- cur_draws[draw_ids,]
  
  for (dose in cur_doses) {
    
    # extract draws for successful infections only
    success_draws <- cur_draws[success_ids[,col_idx],]

    # compute infection duration    
    duration <- success_draws %>%
      mutate(init_N = exp(log_init_RNA(hit_prob, dose, alpha_N, log_N5))) %>%
      mutate(peak_t = peak_time(hit_prob, dose, alpha_t, log_tp5)) %>%
      mutate(log_peak_N = log(init_N * exp(g * peak_t))) %>%
      rowwise() %>%
      mutate(duration = inf_dur(log_peak_N, peak_t, g, d)) %>%
      pull(duration)
    
    # save infection durations
    final_studies <- c(final_studies, rep(study_name, length(duration)))
    final_doses <- c(final_doses, rep(dose, length(duration)))
    final_duration <- c(final_duration, duration)
    
    col_idx <- col_idx + 1
    
  }
  
}

# combine and plot
df_duration <- data.frame(study_name = final_studies,
                          TCID_dose = final_doses,
                          duration = final_duration)

# set study values
study_labels <- c("Li et al. (2005)",
                  "Lawler et al. (2006)",
                  "Nagata et al. (2007)",
                  "Chen et al. (2008)",
                  "Liu et al. (2016)",
                  "Liu et al. (2019)")
df_duration <- df_duration %>%
  mutate(study = case_when(
    (study_name == "BL") ~ study_labels[1],
    (study_name == "JL") ~ study_labels[2],
    (study_name == "NN") ~ study_labels[3],
    (study_name == "YC") ~ study_labels[4],
    (study_name == "LL16") ~ study_labels[5],
    (study_name == "LL19") ~ study_labels[6]
  ))
# set study order
study_order <- c("Li et al. (2005)",
                 "Lawler et al. (2006)",
                 "Nagata et al. (2007)",
                 "Chen et al. (2008)",
                 "Liu et al. (2016)",
                 "Liu et al. (2019)")
df_duration$study <- factor(df_duration$study, levels = study_order)

# set facet values
facet_labels <- c("10^3 ~ TCID[50]",
                  "10^5 ~ TCID[50]",
                  "10^6 ~ TCID[50]",
                  "10^6.9 ~ TCID[50]")
df_duration <- df_duration %>%
  mutate(facet = case_when(
    (study_name == "BL") ~ facet_labels[2],
    (study_name == "JL") ~ facet_labels[4],
    (study_name == "LL16") ~ facet_labels[2],
    (study_name == "LL19") ~ facet_labels[2],
    (study_name == "NN") & (TCID_dose == 1e+3) ~ facet_labels[1],
    (study_name == "NN") & (TCID_dose == 1e+6) ~ facet_labels[3],
    (study_name == "YC") ~ facet_labels[2]
  ))
# set facet order
facet_order <- c("10^3 ~ TCID[50]",
                 "10^5 ~ TCID[50]",
                 "10^6 ~ TCID[50]",
                 "10^6.9 ~ TCID[50]")
df_duration$facet <- factor(df_duration$facet, levels = facet_order)

pXd <- ggplot(df_duration %>% filter(duration != 0), aes(x = study, y = duration, color = interaction(study, facet))) +
  stat_pointinterval(.width = c(0.66, 0.95)) +
  scale_y_continuous(trans = "log10",
                     breaks = c(1e-2, 1e-1, 1e+0, 1e+1, 1e+2, 1e+3, 1e+4, 1e+5, 1e+6, 1e+7),
                     labels = trans_format("log10", math_format(10^.x)),
                     limits = c(1e+0, 1e+2)) +
  scale_color_manual(values = colors) +
  # scale_color_brewer(palette = "Dark2") +
  xlab(element_blank()) +
  ylab("Duration of Infection (Days)") +
  facet_grid(.~facet, labeller = label_parsed,
             scales = "free", space = "free_y", switch = "y") +
  theme(legend.position = "none",
        panel.spacing = unit(0, "in"),
        axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))

plot(pXd)

## Fig 3E: area under the curve distributions ----------------------------------

# area under curve
# INPUTS:
# (N_0): initial viral load
# (t_p): peak time
# (g): growth rate
# (d): decay rate
AUC <- function(N_0, t_p, g, d) {
  
  # area before peak time
  a1 <- (N_0 / g) * (exp(g * t_p) - 1)
  
  # area after peak time
  a2 <- (N_0 * exp(g * t_p)) / d
  
  return(a1 + a2)
}

final_studies <- c()
final_doses <- c()
final_AUCs <- c()

col_idx <- 1 # to keep track of study+dose combination

for (study_name in names(study_names_doses)) {
  
  # extract dose levels for current study
  cur_doses <- study_names_doses[[study_name]]
  
  # filter draws by study name
  cur_study <- study_name
  cur_draws <- draws %>% filter(study_name == cur_study)
  cur_draws <- cur_draws[draw_ids,]
  
  for (dose in cur_doses) {
    
    # extract draws for successful infections only
    success_draws <- cur_draws[success_ids[,col_idx],]
    
    # compute AUCs
    AUCs <- success_draws %>%
      mutate(init_N = exp(log_init_RNA(hit_prob, dose, alpha_N, log_N5))) %>%
      mutate(peak_t = peak_time(hit_prob, dose, alpha_t, log_tp5)) %>%
      rowwise() %>%
      mutate(AUC = AUC(init_N, peak_t, g, d)) %>%
      pull(AUC)
    
    # save AUCs
    final_studies <- c(final_studies, rep(study_name, length(AUCs)))
    final_doses <- c(final_doses, rep(dose, length(AUCs)))
    final_AUCs <- c(final_AUCs, AUCs)
    
    col_idx <- col_idx + 1
    
  }
  
}

# combine and plot
df_AUCs <- data.frame(study_name = final_studies,
                      TCID_dose = final_doses,
                      AUC = final_AUCs)

# set study values
study_labels <- c("Li et al. (2005)",
                  "Lawler et al. (2006)",
                  "Nagata et al. (2007)",
                  "Chen et al. (2008)",
                  "Liu et al. (2016)",
                  "Liu et al. (2019)")
df_AUCs <- df_AUCs %>%
  mutate(study = case_when(
    (study_name == "BL") ~ study_labels[1],
    (study_name == "JL") ~ study_labels[2],
    (study_name == "NN") ~ study_labels[3],
    (study_name == "YC") ~ study_labels[4],
    (study_name == "LL16") ~ study_labels[5],
    (study_name == "LL19") ~ study_labels[6]
  ))
# set study order
study_order <- c("Li et al. (2005)",
                 "Lawler et al. (2006)",
                 "Nagata et al. (2007)",
                 "Chen et al. (2008)",
                 "Liu et al. (2016)",
                 "Liu et al. (2019)")
df_AUCs$study <- factor(df_AUCs$study, levels = study_order)

# set facet values
facet_labels <- c("10^3 ~ TCID[50]",
                  "10^5 ~ TCID[50]",
                  "10^6 ~ TCID[50]",
                  "10^6.9 ~ TCID[50]")
df_AUCs <- df_AUCs %>%
  mutate(facet = case_when(
    (study_name == "BL") ~ facet_labels[2],
    (study_name == "JL") ~ facet_labels[4],
    (study_name == "LL16") ~ facet_labels[2],
    (study_name == "LL19") ~ facet_labels[2],
    (study_name == "NN") & (TCID_dose == 1e+3) ~ facet_labels[1],
    (study_name == "NN") & (TCID_dose == 1e+6) ~ facet_labels[3],
    (study_name == "YC") ~ facet_labels[2]
  ))
# set facet order
facet_order <- c("10^3 ~ TCID[50]",
                 "10^5 ~ TCID[50]",
                 "10^6 ~ TCID[50]",
                 "10^6.9 ~ TCID[50]")
df_AUCs$facet <- factor(df_AUCs$facet, levels = facet_order)

pXe <- ggplot(df_AUCs, aes(x = study, y = AUC, color = interaction(study, facet))) +
  stat_pointinterval(.width = c(0.66, 0.95)) +
  scale_y_continuous(trans = "log10",
                     limits = c(1e+0, 1e+6),
                     labels = trans_format("log10", math_format(10^.x))) +
  scale_color_manual(values = colors) +
  # scale_color_brewer(palette = "Dark2") +
  xlab(element_blank()) +
  ylab("Area Under the Curve") +
  facet_grid(.~facet, labeller = label_parsed,
             scales = "free", space = "free_y", switch = "y") +
  theme(legend.position = "none",
        panel.spacing = unit(0, "in"),
        axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))

plot(pXe)

## Fig 3: combined kinetics results plot ---------------------------------------

fig <- (free(pXa) | ((pXb | pXc) / (pXd | pXe))) + plot_layout(widths = c(1.5, 2)) + plot_annotation(tag_levels="A")

## Save plot -------------------------------------------------------------------

ggsave("outputs/figures/Fig3.tiff", fig, width = 16, height = 8)
