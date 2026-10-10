# Data aggregation and mechanistic modeling enable dose-response analysis of SARS-CoV-1 in non-human primates

Philip Lee(1), Celine Snedden(2), Dylan Morris(2), James Lloyd-Smith(2, 3)

1. Center for Computational Biology, University of California Berkeley, Berkeley, CA, USA
2. Department of Ecology and Evolutionary Biology, University of California Los Angeles, Los Angeles, CA, USA
3. Department of Computational Medicine, University of California Los Angeles, Los Angeles, CA, USA

## Repository information

This repository provides the data and code to reproduce the analyses presented in the associated paper (_Paper information to come_). 

## Article abstract

Dose-response modeling provides estimates of infectious and lethal doses, which can be used to inform control and prevention measures. Unfortunately, data from experimental challenge studies, which are needed to perform dose-response modeling, are often sparse. For example, non-human primate (NHP) challenge studies tend to have small samples sizes and little dose variation, often with only one or two dose levels per study. Thus, it is unfeasible to apply traditional dose-response modeling approaches to data from single NHP studies. To address this challenge, we developed a mechanistic Bayesian model that aggregates and analyzes NHP pathogen load data across multiple studies. Our model links dose-infectivity to pathogen kinetics, which allows us to estimate the infectious dose and evaluate dose effects on within-host viral kinetics simultaneously. With this model, we performed a dose-response analysis for SARS-CoV-1 in NHPs using data compiled from six NHP challenge studies. Our work demonstrates the value in reusing previous data from animal experiments. Our modeling framework can be applied to other pathogens, enabling robust dose-response inference when individual challenge studies are inconclusive.

## Citation information

If you use any of the code or data associated with this manuscript, please cite our work. (_Citation information to come_).

## File structure

All content is organized into the following folders:

- `code`: Contains all code, including the Stan models, model fitting, and figure generation.
  - Analyses are ordered by increasing numbers (e.g., `01-fit-all-data.R`, `02-fit-nasal-data.R`). These files should be run in ascending numerical order.
  - Files that generate figures are preceded with `fig` and the corresponding number from the paper.

- `data`: Contains data as `.csv` files, including the aggregated NHP viral load data (`clean-data.csv`) and additional demographic information about each individual NHP (`additional-data.csv`).

- `outputs`: Contains all output files, including model fits and figures.
