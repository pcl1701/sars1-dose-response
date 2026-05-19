rm(list = ls())

## Install and load necessary packages -----------------------------------------

req_pkgs <- c("dplyr", "tidyr", "ggplot2", "scales", "patchwork")
new_pkgs <- req_pkgs[!(req_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) {
  install.packages(new_pkgs, dependences = TRUE)
}
invisible(lapply(req_pkgs, library, character.only = TRUE))

## Load data -------------------------------------------------------------------

data <- read.csv("data/additional-data.csv")
data$individual <- as.factor(data$individual)
data$animal_species <- as.factor(data$animal_species)
data$data_type <- as.factor(data$data_type)
data$sex <- as.factor(data$sex)
data$age <- as.factor(data$age)

## Fig 1A: data summary bar plots ----------------------------------------------

theme_set(theme_bw())

data_final <- NA
indivs <- unique(data$individual)
for (indiv in indivs) {
  data_indiv <- data %>% filter(individual == indiv) %>% distinct()
  if (nrow(data_indiv) == 1) {
    data_indiv$data_type <- "Oral ~ Only"
    data_final <- rbind(data_final, data_indiv)
  }
  else {
    data_indiv$data_type <- "Oral ~ '&' ~ Nasal"
    data_final <- rbind(data_final, data_indiv[1,])
  }
}
data_final <- data_final[2:nrow(data_final),]

# species ordering
data_final <- data_final %>%
  mutate(animal_species = case_when(animal_species == "RM" ~ "Rhesus",
                                    animal_species == "CM" ~ "Cynomolgus"))
data_final$animal_species <- factor(data_final$animal_species, levels = c("Rhesus", "Cynomolgus"))

# dose ordering
data_final <- data_final %>%
  mutate(TCID_dose = case_when(TCID_dose == 1000 ~ "10^3",
                               TCID_dose == 100000 ~ "10^5",
                               TCID_dose == 1000000 ~ "10^6",
                               TCID_dose == 8571429 ~ "10^6.9"))
data_final$TCID_dose <- factor(data_final$TCID_dose,
                               levels = c("10^3",
                                          "10^5",
                                          "10^6",
                                          "10^6.9"))

# data type ordering
data_final$data_type <- factor(data_final$data_type, levels = c("Oral ~ Only", "Oral ~ '&' ~ Nasal"))

# sex ordering
data_final <- data_final %>%
  mutate(sex = case_when(sex == "F" ~ "Female",
                         sex == "M" ~ "Male"))
data_final$sex <- factor(data_final$sex, levels = c("Female", "Male"))

# age ordering
data_final <- data_final %>%
  mutate(age = case_when(age == "adult" ~ "Adult",
                         TRUE ~ "Juvenile"))
data_final$age <- factor(data_final$age, levels = c("Juvenile", "Adult"))

# facet ordering
data_final <- data_final %>%
  pivot_longer(!individual, names_to = "facet", values_to = "value")
data_final <- data_final %>%
  mutate(facet = case_when(facet == "age" ~ "Age",
                           facet == "animal_species" ~ "Species",
                           facet == "data_type" ~ "Swab Type",
                           facet == "sex" ~ "Sex",
                           facet == "TCID_dose" ~ "Inoculation Dose"))
data_final$facet <- factor(data_final$facet,
                           levels = c("Inoculation Dose",
                                      "Swab Type",
                                      "Species",
                                      "Sex",
                                      "Age"))

# plot bar plots
p1 <- ggplot(data_final, aes(x = value)) +
  geom_bar() +
  geom_text(aes(label = after_stat(count)), stat = "count", nudge_y = 1, size = 2.5) +
  scale_x_discrete(labels = str2expression) +
  scale_y_continuous(breaks = c(0*7.8, 1*7.8, 2*7.8, 3*7.8, 4*7.8, 5*7.8),
                     labels = c("0%", "20%", "40%", "60%", "80%", "100%")) +
  xlab(element_blank()) +
  ylab("Individuals") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  facet_grid(~facet, scales = "free_x", space = "free")

plot(p1)

## Fig 1B: data summary bar plots ----------------------------------------------

theme_set(theme_bw())

# species ordering
data_final <- data %>%
  mutate(animal_species = case_when(animal_species == "RM" ~ "Rhesus",
                                   animal_species == "CM" ~ "Cynomolgus"))
data_final$animal_species <- factor(data_final$animal_species, levels = c("Rhesus", "Cynomolgus"))

# dose ordering
data_final <- data_final %>%
  mutate(TCID_dose = case_when(TCID_dose == 1000 ~ "10^3",
                               TCID_dose == 100000 ~ "10^5",
                               TCID_dose == 1000000 ~ "10^6",
                               TCID_dose == 8571429 ~ "10^6.9"))
data_final$TCID_dose <- factor(data_final$TCID_dose,
                               levels = c("10^3",
                                          "10^5",
                                          "10^6",
                                          "10^6.9"))

# data type ordering
data_final <- data_final %>%
  mutate(data_type = case_when(data_type == "nasal swab" ~ "Nasal",
                               TRUE ~ "Oral"))
data_final$data_type <- factor(data_final$data_type, levels = c("Oral", "Nasal"))

# sex ordering
data_final <- data_final %>%
  mutate(sex = case_when(sex == "F" ~ "Female",
                         sex == "M" ~ "Male"))
data_final$sex <- factor(data_final$sex, levels = c("Female", "Male"))

# age ordering
data_final <- data_final %>%
  mutate(age = case_when(age == "adult" ~ "Adult",
                         TRUE ~ "Juvenile"))
data_final$age <- factor(data_final$age, levels = c("Juvenile", "Adult"))

# facet ordering
data_final <- data_final %>%
  pivot_longer(!individual, names_to = "facet", values_to = "value")
data_final <- data_final %>%
  mutate(facet = case_when(facet == "age" ~ "Age",
                           facet == "animal_species" ~ "Species",
                           facet == "data_type" ~ "Swab Type",
                           facet == "sex" ~ "Sex",
                           facet == "TCID_dose" ~ "Inoculation Dose"))
data_final$facet <- factor(data_final$facet,
                           levels = c("Inoculation Dose",
                                      "Swab Type",
                                      "Species",
                                      "Sex",
                                      "Age"))

# plot bar plots
p2 <- ggplot(data_final, aes(x = value)) +
  geom_bar() +
  geom_text(aes(label = after_stat(count)), stat = "count", nudge_y = 4.5, size = 2.5) +
  scale_x_discrete(labels = str2expression) +
  scale_y_continuous(breaks = c(0*36, 1*36, 2*36, 3*36, 4*36, 5*36),
                     labels = c("0%", "20%", "40%", "60%", "80%", "100%")) +
  xlab(element_blank()) +
  ylab("Observations") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  facet_grid(~facet, scales = "free_x", space = "free")

plot(p2)

## Fig 1: combined data summary plot -------------------------------------------

fig <- (p1 / p2) + plot_annotation(tag_levels = "A")

## Save plot -------------------------------------------------------------------

ggsave("outputs/figures/Fig1.tiff", fig, width = 9, height = 6)
