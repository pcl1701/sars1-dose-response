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


## Overall settings ------------------------------------------------------------

theme_set(theme_bw())

n_datasets <- 10
right_label_size <- 4

true_log10ID50s <- c(2, 3, 4, 5, 6)


## Fig 7A - original doses, uneven distribution --------------------------------

color <- "#F08700"

figs <- vector("list", length(true_log10ID50s))

for (j in 1:length(true_log10ID50s)) {
  
  df_all <- read.csv(paste0(
    "outputs/simulations/data/simulated-data-original-doses-uneven-distribution-fixed-log10ID50-",
    true_log10ID50s[j],
    ".csv"
  ))
  
  
  p1 <- ggplot() +
    geom_vline(
      xintercept = 10^true_log10ID50s[j],
      linetype = "dashed",
      color = "black",
      linewidth = 0.8
    )
  
  for (i in 1:n_datasets) {
    
    # filter by dataset id
    data_id <- i
    data <- df_all %>% filter(dataset_id == data_id)
    data$individual <- as.factor(data$individual)
    data$animal_species <- as.factor(data$animal_species)
    data$study_name <- as.factor(data$study_name)
    data$viral_isolate <- as.factor(data$viral_isolate)
    data$data_type <- as.factor(data$data_type)
    
    fit <- readRDS(paste0(
      "outputs/simulations/fits/fit-simulated-data-original-doses-uneven-distribution-fixed-log10ID50-",
      true_log10ID50s[j],
      "-dataset-",
      data_id,
      ".RDS"
    ))
    fit <- fit %>%
      recover_types(data)
    
    ID50s <- fit %>%
      spread_draws(log_ID50) %>%
      mutate(ID50 = exp(log_ID50)) %>%
      pull(ID50)
    
    df_ID50 <- data.frame(ID50 = ID50s)
    
    p1 <- p1 +
      geom_density(data = df_ID50, aes(x = ID50, color = "Simulated Data"), fill = NA, linewidth = 0.4, key_glyph = "path")
    
  }
  
  p1 <- p1 +
    scale_color_manual(
      name = NULL,
      values = c("Simulated Data" = color)
    ) +
    scale_x_continuous(trans = "log10",
                       limits = c(1e-1, 1e+7),
                       breaks = c(1e-1, 1e+0, 1e+1, 1e+2, 1e+3, 1e+4, 1e+5, 1e+6, 1e+7),
                       labels = trans_format("log10", math_format(10^.x))) +
    xlab(bquote(ID[50] ~ "(" * TCID[50] * ")")) +
    ylab("Density") +
    theme(legend.position = "none")
  
  
  # simulated datasets
  df_sim <- data.frame()
  
  df_all <- read.csv(paste0(
    "outputs/simulations/data/simulated-data-original-doses-uneven-distribution-fixed-log10ID50-",
    true_log10ID50s[j],
    ".csv"
  ))
  
  for (i in 1:n_datasets) {
    
    data_id <- i
    data <- df_all %>% filter(dataset_id == data_id)
    
    data$individual <- as.factor(data$individual)
    data$animal_species <- as.factor(data$animal_species)
    data$study_name <- as.factor(data$study_name)
    data$viral_isolate <- as.factor(data$viral_isolate)
    data$data_type <- as.factor(data$data_type)
    
    fit <- readRDS(paste0(
      "outputs/simulations/fits/fit-simulated-data-original-doses-uneven-distribution-fixed-log10ID50-",
      true_log10ID50s[j],
      "-dataset-",
      data_id,
      ".RDS"
    ))
    
    fit <- fit %>%
      recover_types(data)
    
    ID50s <- fit %>%
      spread_draws(log_ID50) %>%
      mutate(ID50 = exp(log_ID50)) %>%
      pull(ID50)
    
    df_sim <- rbind(
      df_sim,
      data.frame(
        dataset = paste0("Simulated ", i),
        ID50 = ID50s
      )
    )
  }
  
  df_all_ID50 <- df_sim
  
  df_CI <- df_all_ID50 %>%
    group_by(dataset) %>%
    summarize(
      median = median(ID50),
      lower_66 = quantile(ID50, 0.17),
      upper_66 = quantile(ID50, 0.83),
      lower_95 = quantile(ID50, 0.025),
      upper_95 = quantile(ID50, 0.975)
    )
  
  df_CI <- df_CI %>%
    mutate(
      x = as.numeric(factor(dataset))
    )
  
  df_CI$group <- "Simulated Data"

  p2 <- ggplot(df_CI, aes(x = x, y = median)) +
    
    # true ID50
    geom_hline(
      yintercept = 10^true_log10ID50s[j],
      linetype = "dashed",
      color = "black",
      linewidth = 0.8
    ) +
    
    # 95% interval
    geom_errorbar(
      aes(ymin = lower_95, ymax = upper_95, color = group),
      linewidth = 0.5,
      width = 0
    ) +
    
    # 66% interval
    geom_errorbar(
      aes(ymin = lower_66, ymax = upper_66, color = group),
      linewidth = 2,
      width = 0
    ) +
    
    # median point
    geom_point(
      aes(color = group),
      size = 3
    ) +
    
    scale_color_manual(
      values = c(
        "Simulated Data" = color
      ),
      guide = "none"
    ) +
    
    scale_x_continuous(
      breaks = c(mean(df_CI$x)),
      labels = c("Simulated Data")
    ) +
    
    scale_y_log10(
      limits = c(1e-1, 1e7),
      breaks = c(1e-1, 1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7),
      labels = trans_format("log10", math_format(10^.x))
    ) +
    
    xlab(NULL) +
    ylab(bquote(ID[50] ~ "(" * TCID[50] * ")")) +
    
    theme_bw() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
  
  p3 <- ggplot() +
    annotate(
      "text",
      x = 0.5,
      y = 0.5,
      label = deparse(
        bquote(True~ID[50] == 10^.(true_log10ID50s[j]))
      ),
      parse = TRUE,
      angle = 270,
      size = right_label_size
    ) +
    theme_void()
  
  fig <- p1 + p2 + p3 +
    plot_layout(widths = c(3, 2, 0.3))
  
  figs[[j]] <- fig

}

# dose distribution plot
indiv_doses <- c(rep(10^3, 2), rep(10^5, 33), rep(10^6, 2), rep(8571429, 2))
df_doses <- data.frame(
  dose = indiv_doses
) %>% count(dose)

# dose plot
fig_top <- ggplot(df_doses, aes(x = factor(dose), y = n)) +
  geom_col(fill = color) +
  scale_x_discrete(
    labels = c(
      expression(10^3),
      expression(10^5),
      expression(10^6),
      expression(10^6.9)
    )
  ) +
  labs(
    x = expression("Inoculation Dose"~"(TCID"[50]*")"),
    y = "Number of individuals",
    fill = NULL
  ) + theme(legend.position = "none")

fig_bot <- wrap_plots(figs, ncol = 1)

fig_A <- wrap_elements((fig_top / fig_bot) + plot_layout(heights = c(1, 5)))


## Fig 7B - original doses, balanced distribution ------------------------------

color <- "#8390FA"

figs <- vector("list", length(true_log10ID50s))

for (j in 1:length(true_log10ID50s)) {
  
  df_all <- read.csv(paste0(
    "outputs/simulations/data/simulated-data-original-doses-balanced-distribution-fixed-log10ID50-",
    true_log10ID50s[j],
    ".csv"
  ))
  
  
  p1 <- ggplot() +
    geom_vline(
      xintercept = 10^true_log10ID50s[j],
      linetype = "dashed",
      color = "black",
      linewidth = 0.8
    )
  
  for (i in 1:n_datasets) {
    
    # filter by dataset id
    data_id <- i
    data <- df_all %>% filter(dataset_id == data_id)
    data$individual <- as.factor(data$individual)
    data$animal_species <- as.factor(data$animal_species)
    data$study_name <- as.factor(data$study_name)
    data$viral_isolate <- as.factor(data$viral_isolate)
    data$data_type <- as.factor(data$data_type)
    
    fit <- readRDS(paste0(
      "outputs/simulations/fits/fit-simulated-data-original-doses-balanced-distribution-fixed-log10ID50-",
      true_log10ID50s[j],
      "-dataset-",
      data_id,
      ".RDS"
    ))
    fit <- fit %>%
      recover_types(data)
    
    ID50s <- fit %>%
      spread_draws(log_ID50) %>%
      mutate(ID50 = exp(log_ID50)) %>%
      pull(ID50)
    
    df_ID50 <- data.frame(ID50 = ID50s)
    
    p1 <- p1 +
      geom_density(data = df_ID50, aes(x = ID50, color = "Simulated Data"), fill = NA, linewidth = 0.4, key_glyph = "path")
    
  }
  
  p1 <- p1 +
    scale_color_manual(
      name = NULL,
      values = c("Simulated Data" = color)
    ) +
    scale_x_continuous(trans = "log10",
                       limits = c(1e-1, 1e+7),
                       breaks = c(1e-1, 1e+0, 1e+1, 1e+2, 1e+3, 1e+4, 1e+5, 1e+6, 1e+7),
                       labels = trans_format("log10", math_format(10^.x))) +
    xlab(bquote(ID[50] ~ "(" * TCID[50] * ")")) +
    ylab("Density") +
    theme(legend.position = "none")
  
  
  # simulated datasets
  df_sim <- data.frame()
  
  df_all <- read.csv(paste0(
    "outputs/simulations/data/simulated-data-original-doses-balanced-distribution-fixed-log10ID50-",
    true_log10ID50s[j],
    ".csv"
  ))
  
  for (i in 1:n_datasets) {
    
    data_id <- i
    data <- df_all %>% filter(dataset_id == data_id)
    
    data$individual <- as.factor(data$individual)
    data$animal_species <- as.factor(data$animal_species)
    data$study_name <- as.factor(data$study_name)
    data$viral_isolate <- as.factor(data$viral_isolate)
    data$data_type <- as.factor(data$data_type)
    
    fit <- readRDS(paste0(
      "outputs/simulations/fits/fit-simulated-data-original-doses-balanced-distribution-fixed-log10ID50-",
      true_log10ID50s[j],
      "-dataset-",
      data_id,
      ".RDS"
    ))
    
    fit <- fit %>%
      recover_types(data)
    
    ID50s <- fit %>%
      spread_draws(log_ID50) %>%
      mutate(ID50 = exp(log_ID50)) %>%
      pull(ID50)
    
    df_sim <- rbind(
      df_sim,
      data.frame(
        dataset = paste0("Simulated ", i),
        ID50 = ID50s
      )
    )
  }
  
  df_all_ID50 <- df_sim
  
  df_CI <- df_all_ID50 %>%
    group_by(dataset) %>%
    summarize(
      median = median(ID50),
      lower_66 = quantile(ID50, 0.17),
      upper_66 = quantile(ID50, 0.83),
      lower_95 = quantile(ID50, 0.025),
      upper_95 = quantile(ID50, 0.975)
    )
  
  df_CI <- df_CI %>%
    mutate(
      x = as.numeric(factor(dataset))
    )
  
  df_CI$group <- "Simulated Data"
  
  p2 <- ggplot(df_CI, aes(x = x, y = median)) +
    
    # true ID50
    geom_hline(
      yintercept = 10^true_log10ID50s[j],
      linetype = "dashed",
      color = "black",
      linewidth = 0.8
    ) +
    
    # 95% interval
    geom_errorbar(
      aes(ymin = lower_95, ymax = upper_95, color = group),
      linewidth = 0.5,
      width = 0
    ) +
    
    # 66% interval
    geom_errorbar(
      aes(ymin = lower_66, ymax = upper_66, color = group),
      linewidth = 2,
      width = 0
    ) +
    
    # median point
    geom_point(
      aes(color = group),
      size = 3
    ) +
    
    scale_color_manual(
      values = c(
        "Simulated Data" = color
      ),
      guide = "none"
    ) +
    
    scale_x_continuous(
      breaks = c(mean(df_CI$x)),
      labels = c("Simulated Data")
    ) +
    
    scale_y_log10(
      limits = c(1e-1, 1e7),
      breaks = c(1e-1, 1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7),
      labels = trans_format("log10", math_format(10^.x))
    ) +
    
    xlab(NULL) +
    ylab(bquote(ID[50] ~ "(" * TCID[50] * ")")) +
    
    theme_bw() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
  
  p3 <- ggplot() +
    annotate(
      "text",
      x = 0.5,
      y = 0.5,
      label = deparse(
        bquote(True~ID[50] == 10^.(true_log10ID50s[j]))
      ),
      parse = TRUE,
      angle = 270,
      size = right_label_size
    ) +
    theme_void()
  
  fig <- p1 + p2 + p3 +
    plot_layout(widths = c(3, 2, 0.3))
  
  figs[[j]] <- fig
  
}

# dose distribution plot
indiv_doses <- c(rep(10^3, 10), rep(10^5, 10), rep(10^6, 10), rep(8571429, 9))
df_doses <- data.frame(
  dose = indiv_doses
) %>% count(dose)

# dose plot
fig_top <- ggplot(df_doses, aes(x = factor(dose), y = n)) +
  geom_col(fill = color) +
  scale_x_discrete(
    labels = c(
      expression(10^3),
      expression(10^5),
      expression(10^6),
      expression(10^6.9)
    )
  ) +
  scale_y_continuous(
    breaks = c(0, 2, 4, 6, 8, 10),
    limits = c(0, 10)
  ) +
  labs(
    x = expression("Inoculation Dose"~"(TCID"[50]*")"),
    y = "Number of individuals",
    fill = NULL
  ) + theme(legend.position = "none")

fig_bot <- wrap_plots(figs, ncol = 1)

fig_B <- wrap_elements((fig_top / fig_bot) + plot_layout(heights = c(1, 5)))


## Fig 7C - new doses, balanced distribution -----------------------------------

color <- "#7BC950"

figs <- vector("list", length(true_log10ID50s))

for (j in 1:length(true_log10ID50s)) {
  
  df_all <- read.csv(paste0(
    "outputs/simulations/data/simulated-data-new-doses-balanced-distribution-fixed-log10ID50-",
    true_log10ID50s[j],
    ".csv"
  ))
  
  
  p1 <- ggplot() +
    geom_vline(
      xintercept = 10^true_log10ID50s[j],
      linetype = "dashed",
      color = "black",
      linewidth = 0.8
    )
  
  for (i in 1:n_datasets) {
    
    # filter by dataset id
    data_id <- i
    data <- df_all %>% filter(dataset_id == data_id)
    data$individual <- as.factor(data$individual)
    data$animal_species <- as.factor(data$animal_species)
    data$study_name <- as.factor(data$study_name)
    data$viral_isolate <- as.factor(data$viral_isolate)
    data$data_type <- as.factor(data$data_type)
    
    fit <- readRDS(paste0(
      "outputs/simulations/fits/fit-simulated-data-new-doses-balanced-distribution-fixed-log10ID50-",
      true_log10ID50s[j],
      "-dataset-",
      data_id,
      ".RDS"
    ))
    fit <- fit %>%
      recover_types(data)
    
    ID50s <- fit %>%
      spread_draws(log_ID50) %>%
      mutate(ID50 = exp(log_ID50)) %>%
      pull(ID50)
    
    df_ID50 <- data.frame(ID50 = ID50s)
    
    p1 <- p1 +
      geom_density(data = df_ID50, aes(x = ID50, color = "Simulated Data"), fill = NA, linewidth = 0.4, key_glyph = "path")
    
  }
  
  p1 <- p1 +
    scale_color_manual(
      name = NULL,
      values = c("Simulated Data" = color)
    ) +
    scale_x_continuous(trans = "log10",
                       limits = c(1e-1, 1e+7),
                       breaks = c(1e-1, 1e+0, 1e+1, 1e+2, 1e+3, 1e+4, 1e+5, 1e+6, 1e+7),
                       labels = trans_format("log10", math_format(10^.x))) +
    xlab(bquote(ID[50] ~ "(" * TCID[50] * ")")) +
    ylab("Density") +
    theme(legend.position = "none")
  
  
  # simulated datasets
  df_sim <- data.frame()
  
  df_all <- read.csv(paste0(
    "outputs/simulations/data/simulated-data-new-doses-balanced-distribution-fixed-log10ID50-",
    true_log10ID50s[j],
    ".csv"
  ))
  
  for (i in 1:n_datasets) {
    
    data_id <- i
    data <- df_all %>% filter(dataset_id == data_id)
    
    data$individual <- as.factor(data$individual)
    data$animal_species <- as.factor(data$animal_species)
    data$study_name <- as.factor(data$study_name)
    data$viral_isolate <- as.factor(data$viral_isolate)
    data$data_type <- as.factor(data$data_type)
    
    fit <- readRDS(paste0(
      "outputs/simulations/fits/fit-simulated-data-new-doses-balanced-distribution-fixed-log10ID50-",
      true_log10ID50s[j],
      "-dataset-",
      data_id,
      ".RDS"
    ))
    
    fit <- fit %>%
      recover_types(data)
    
    ID50s <- fit %>%
      spread_draws(log_ID50) %>%
      mutate(ID50 = exp(log_ID50)) %>%
      pull(ID50)
    
    df_sim <- rbind(
      df_sim,
      data.frame(
        dataset = paste0("Simulated ", i),
        ID50 = ID50s
      )
    )
  }
  
  df_all_ID50 <- df_sim
  
  df_CI <- df_all_ID50 %>%
    group_by(dataset) %>%
    summarize(
      median = median(ID50),
      lower_66 = quantile(ID50, 0.17),
      upper_66 = quantile(ID50, 0.83),
      lower_95 = quantile(ID50, 0.025),
      upper_95 = quantile(ID50, 0.975)
    )
  
  df_CI <- df_CI %>%
    mutate(
      x = as.numeric(factor(dataset))
    )
  
  df_CI$group <- "Simulated Data"
  
  p2 <- ggplot(df_CI, aes(x = x, y = median)) +
    
    # true ID50
    geom_hline(
      yintercept = 10^true_log10ID50s[j],
      linetype = "dashed",
      color = "black",
      linewidth = 0.8
    ) +
    
    # 95% interval
    geom_errorbar(
      aes(ymin = lower_95, ymax = upper_95, color = group),
      linewidth = 0.5,
      width = 0
    ) +
    
    # 66% interval
    geom_errorbar(
      aes(ymin = lower_66, ymax = upper_66, color = group),
      linewidth = 2,
      width = 0
    ) +
    
    # median point
    geom_point(
      aes(color = group),
      size = 3
    ) +
    
    scale_color_manual(
      values = c(
        "Simulated Data" = color
      ),
      guide = "none"
    ) +
    
    scale_x_continuous(
      breaks = c(mean(df_CI$x)),
      labels = c("Simulated Data")
    ) +
    
    scale_y_log10(
      limits = c(1e-1, 1e7),
      breaks = c(1e-1, 1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7),
      labels = trans_format("log10", math_format(10^.x))
    ) +
    
    xlab(NULL) +
    ylab(bquote(ID[50] ~ "(" * TCID[50] * ")")) +
    
    theme_bw() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
  
  p3 <- ggplot() +
    annotate(
      "text",
      x = 0.5,
      y = 0.5,
      label = deparse(
        bquote(True~ID[50] == 10^.(true_log10ID50s[j]))
      ),
      parse = TRUE,
      angle = 270,
      size = right_label_size
    ) +
    theme_void()
  
  fig <- p1 + p2 + p3 +
    plot_layout(widths = c(3, 2, 0.3))
  
  figs[[j]] <- fig
  
}

# dose distribution plot
indiv_doses <- c(rep(10^2, 8), rep(10^3, 8), rep(10^4, 8), rep(10^5, 8), rep(10^6, 7))
df_doses <- data.frame(
  dose = indiv_doses
) %>% count(dose)

# dose plot
fig_top <- ggplot(df_doses, aes(x = factor(dose), y = n)) +
  geom_col(fill = color) +
  scale_x_discrete(
    labels = c(
      expression(10^2),
      expression(10^3),
      expression(10^4),
      expression(10^5),
      expression(10^6)
    )
  ) +
  labs(
    x = expression("Inoculation Dose"~"(TCID"[50]*")"),
    y = "Number of individuals",
    fill = NULL
  ) + theme(legend.position = "none")

fig_bot <- wrap_plots(figs, ncol = 1)

fig_C <- wrap_elements((fig_top / fig_bot) + plot_layout(heights = c(1, 5)))

## Combined plot ---------------------------------------------------------------

figs <- list(fig_A, fig_B, fig_C)

final_fig <- wrap_plots(figs, ncol=3) +
  plot_annotation(tag_levels = "A")


ggsave("outputs/figures/Fig7.tiff", final_fig, width = 15, height = 12)

