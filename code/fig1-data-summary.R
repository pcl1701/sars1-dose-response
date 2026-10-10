rm(list = ls())

## Install and load necessary packages -----------------------------------------

req_pkgs <- c(
  "dplyr",
  "ggplot2",
  "ComplexUpset",
  "patchwork"
)

new_pkgs <- req_pkgs[!(req_pkgs %in% installed.packages()[, "Package"])]

if (length(new_pkgs)) {
  install.packages(new_pkgs, dependencies = TRUE)
}

invisible(lapply(req_pkgs, library, character.only = TRUE))


## Load data -------------------------------------------------------------------

data <- read.csv("data/additional-data.csv")

data$individual <- as.factor(data$individual)
data$animal_species <- as.factor(data$animal_species)
data$data_type <- as.factor(data$data_type)
data$sex <- as.factor(data$sex)
data$age <- as.factor(data$age)


## Recode variables ------------------------------------------------------------

data <- data %>%
  mutate(
    
    # Species
    animal_species = case_when(
      animal_species == "RM" ~ "Rhesus",
      animal_species == "CM" ~ "Cynomolgus"
    ),
    
    # Dose
    TCID_dose = case_when(
      TCID_dose == 1000 ~ "10³ Dose",
      TCID_dose == 100000 ~ "10⁵ Dose",
      TCID_dose == 1000000 ~ "10⁶ Dose",
      TCID_dose == 8571429 ~ "10⁶·⁹ Dose"
    ),
    
    # Swab type
    data_type = case_when(
      data_type == "nasal swab" ~ "Nasal",
      TRUE ~ "Oral"
    ),
    
    # Sex
    sex = case_when(
      sex == "F" ~ "Female",
      sex == "M" ~ "Male"
    ),
    
    # Age
    age = case_when(
      age == "adult" ~ "Adult",
      TRUE ~ "Juvenile"
    )
  )


## Fig 1A: individual-level experimental design -------------------------------

# Collapse to one row per individual.
#
# Oral  = individual has at least one oral sample
# Nasal = individual has at least one nasal sample
#
# Thus:
#   Oral only       = Oral TRUE,  Nasal FALSE
#   Oral + Nasal     = Oral TRUE,  Nasal TRUE
#   Nasal only       = Oral FALSE, Nasal TRUE

data_individual <- data %>%
  group_by(individual) %>%
  summarise(
    animal_species = first(animal_species),
    TCID_dose = first(TCID_dose),
    sex = first(sex),
    age = first(age),
    
    Oral = any(data_type == "Oral"),
    Nasal = any(data_type == "Nasal"),
    
    .groups = "drop"
  ) %>%
  mutate(
    
    # Species sets
    Rhesus = animal_species == "Rhesus",
    Cynomolgus = animal_species == "Cynomolgus",
    
    # Dose sets
    `10³ Dose` = TCID_dose == "10³ Dose",
    `10⁵ Dose` = TCID_dose == "10⁵ Dose",
    `10⁶ Dose` = TCID_dose == "10⁶ Dose",
    `10⁶·⁹ Dose` = TCID_dose == "10⁶·⁹ Dose",
    
    # Sex sets
    Female = sex == "Female",
    Male = sex == "Male",
    
    # Age sets
    Juvenile = age == "Juvenile",
    Adult = age == "Adult"
  )


# Sets to display
sets <- c(
  "10³ Dose",
  "10⁵ Dose",
  "10⁶ Dose",
  "10⁶·⁹ Dose",
  "Oral",
  "Nasal",
  "Rhesus",
  "Cynomolgus",
  "Female",
  "Male",
  "Juvenile",
  "Adult"
)


p1 <- upset(
  data_individual,
  intersect = sets,
  min_size = 1,
  sort_intersections = "descending",
  set_sizes = upset_set_size(),
  width_ratio = 0.25,
  name = "Experimental characteristics"
)


## Fig 1B: observation-level experimental design -------------------------------

# Each row in the original data is an observation.
#
# Unlike Fig 1A, Oral and Nasal are mutually exclusive here because each
# observation is one swab type.

data_observation <- data %>%
  mutate(
    
    # Species sets
    Rhesus = animal_species == "Rhesus",
    Cynomolgus = animal_species == "Cynomolgus",
    
    # Dose sets
    `10³ Dose` = TCID_dose == "10³ Dose",
    `10⁵ Dose` = TCID_dose == "10⁵ Dose",
    `10⁶ Dose` = TCID_dose == "10⁶ Dose",
    `10⁶·⁹ Dose` = TCID_dose == "10⁶·⁹ Dose",
    
    # Swab sets
    Oral = data_type == "Oral",
    Nasal = data_type == "Nasal",
    
    # Sex sets
    Female = sex == "Female",
    Male = sex == "Male",
    
    # Age sets
    Juvenile = age == "Juvenile",
    Adult = age == "Adult"
  )


p2 <- upset(
  data_observation,
  intersect = sets,
  min_size = 1,
  sort_intersections = "descending",
  set_sizes = upset_set_size(),
  width_ratio = 0.25,
  name = "Experimental characteristics"
)


## Combine ---------------------------------------------------------------------

fig <- (wrap_elements(p1) + wrap_elements(p2)) + plot_annotation(tag_levels = "A")

plot(fig)

## Save ------------------------------------------------------------------------

ggsave("outputs/figures/Fig1.tiff", fig, width = 12, height = 6)
