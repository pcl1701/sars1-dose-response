rm(list = ls())

## Install and load necessary packages -----------------------------------------

req_pkgs <- c("dplyr", "ggplot2", "scales")
new_pkgs <- req_pkgs[!(req_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) {
  install.packages(new_pkgs, dependences = TRUE)
}
invisible(lapply(req_pkgs, library, character.only = TRUE))

## Load data -------------------------------------------------------------------

data <- read.csv("data/clean-data.csv")
data$individual <- as.factor(data$individual)
data$animal_species <- as.factor(data$animal_species)
data$study_name <- as.factor(data$study_name)
data$viral_isolate <- as.factor(data$viral_isolate)
data$data_type <- as.factor(data$data_type)

## Fig S2: individual observed trajectories ------------------------------------

theme_set(theme_bw())

# plot settings
oral_color <- "#222E50"
nasal_color <- "#EC0B43"

# prepare data for plotting
indiv_data <- data %>%
  filter(log10_viral_load != 100) %>%
  mutate(dose = case_when(study_name == "JL" ~ "10^6.9 ~ TCID[50]",
                          study_name == "YC" ~ "10^5 ~ TCID[50]",
                          study_name == "BL" ~ "10^5 ~ TCID[50]",
                          study_name == "LL16" ~ "10^5 ~ TCID[50]",
                          study_name == "LL19" ~ "10^5 ~ TCID[50]",
                          study_name == "NN" & TCID_dose == 1e+3 ~ "10^3 ~ TCID[50]",
                          study_name == "NN" & TCID_dose == 1e+6 ~ "10^6 ~ TCID[50]")) %>%
  mutate(data_type = case_when(data_type == "pharyngeal swab" ~ "oral",
                               data_type == "nasal swab" ~ "nasal",
                               data_type == "oropharyngeal swab" ~ "oral",
                               data_type == "oral swab" ~ "oral",
                               data_type == "throat swab" ~ "oral")) %>%
  mutate(below_LOD = case_when(log10_viral_load == -100 ~ 1,
                               TRUE ~ 0)) %>%
  select(individual, day_post_infection, data_type, log10_viral_load,
         log10_LOD, dose, below_LOD, study_name) %>%
  mutate(facet = paste("'", individual, ",'", " ~ ", dose, sep = ""))

# manually set points below LOD
indiv_data <- indiv_data %>%
  mutate(log10_viral_load = case_when(study_name == "JL" & data_type == "oral" & below_LOD == 1 ~ -0.4,
                                      study_name == "JL" & data_type == "nasal" & below_LOD == 1 ~ -0.8,
                                      study_name == "NN" & data_type == "oral" & below_LOD == 1 ~ -0.4,
                                      study_name == "NN" & data_type == "nasal" & below_LOD == 1 ~ -0.8,
                                      individual == "LL_RM-14" & below_LOD == 1 ~ -0.6,
                                      TRUE ~ log10_viral_load))

# set up plot
# add shading, zero line, LOD lines
p <- ggplot() +
  geom_ribbon(data = data.frame(x = seq(-5, 40, 1), y = 1, lower = 1e-2, upper = 1e+0),
              aes(ymin = lower, ymax = upper, x = x, y = y),
              alpha = 0.2) +
  geom_hline(yintercept = 1e+0) +
  geom_hline(data = indiv_data,
             aes(yintercept = 10^log10_LOD),
             linetype = "dashed",
             alpha = 0.5)

# add individual data
p <- p +
  geom_line(data = indiv_data,
            aes(x = day_post_infection, y = 10^log10_viral_load,
                group = data_type, color = data_type),
            alpha = 0.5,
            linewidth = 0.3) +
  geom_point(data = indiv_data,
             aes(x = day_post_infection, y = 10^log10_viral_load,
                 fill = data_type, shape = as.factor(below_LOD)),
             alpha = 0.5,
             size = 2,
             stroke = 0)

# final plot settings
p <- p +
  facet_wrap(~ facet, labeller = label_parsed) +
  scale_x_continuous(breaks = seq(0, 40, 5)) +
  scale_y_continuous(trans = "log10",
                     breaks = c(1e+0, 1e+1, 1e+2, 1e+3, 1e+4, 1e+5, 1e+6, 1e+7),
                     labels = trans_format("log10", math_format(10^.x))) +
  coord_cartesian(xlim = c(0, 30), ylim = c(10^-0.8, 10^7)) +
  xlab("Day Post Infection") +
  ylab("Viral Load (RNA Copies/Sample)") +
  theme(legend.position = "bottom") +
  scale_shape_manual(values = c(21, 25),
                     guide = "none") +
  scale_color_manual(values = c("oral" = oral_color,
                                "nasal" = nasal_color),
                     name = element_blank(),
                     labels = c("oral" = "Oral",
                                "nasal" = "Nasal")) +
  scale_fill_manual(values = c("oral" = oral_color,
                               "nasal" = nasal_color),
                    name = element_blank(),
                    labels = c("oral" = "Oral",
                               "nasal" = "Nasal")) +
  guides(fill = guide_legend(override.aes = list(shape=21)))

plot(p)

## Save plot -------------------------------------------------------------------

ggsave('outputs/figures/S2 Fig.tiff', p, width = 10, height = 10.5)
