# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(flextable) # install.packages("flextable")
library(officer) # install.packages("officer")
library(janitor)
# Bioconductor libraries
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
# source plot function
source("C:/R/CD38-effect-of-treatment/code/functions/plot.gg_violin_interaction.r")
source("C:/R/CD38-effect-of-treatment/code/functions/plot.gg_patient_pairs_interaction.r")
# load reference set
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/data_scores_k1208.RData")



# DEFINE SET ####
baseline_K1208 <- data_scores_k1208 %>%
    dplyr::filter(Group == "Index", Patient %nin% c(15, 18)) %>%
    mutate(
        MMDx = case_when(
            Patient == 1 ~ "Severe fully-developed ABMR",
            Patient == 2 ~ "Severe fully-developed ABMR",
            Patient == 3 ~ "Mixed rejection",
            Patient == 4 ~ "Moderate fully-developed ABMR",
            Patient == 5 ~ "Severe fully-developed ABMR",
            Patient == 6 ~ "Moderate fully-developed ABMR",
            Patient == 7 ~ "Severe fully-developed ABMR",
            Patient == 8 ~ "Late-stage ABMR",
            Patient == 9 ~ "Late-stage ABMR",
            Patient == 10 ~ "Severe fully-developed ABMR",
            Patient == 11 ~ "Moderate to severe early-stage ABMR",
            Patient == 12 ~ "Moderate fully-developed ABMR",
            Patient == 13 ~ "Severe early-stage ABMR",
            Patient == 14 ~ "Severe early-stage ABMR",
            Patient == 16 ~ "Severe fully-developed ABMR",
            Patient == 17 ~ "Severe fully-developed ABMR",
            Patient == 19 ~ "Late-stage ABMR",
            Patient == 20 ~ "Severe fully-developed ABMR",
            Patient == 21 ~ "Moderate fully-developed ABMR",
            Patient == 22 ~ "Moderate to severe early-stage ABMR"
        ) %>% factor(levels = c(
            "Moderate to severe early-stage ABMR",
            "Severe early-stage ABMR",
            "Moderate fully-developed ABMR",
            "Severe fully-developed ABMR",
            "Late-stage ABMR",
            "Mixed rejection"
        )),
        .after = "Patient"
    )



baseline_K1208 %>%
    dplyr::select(Patient, CEL, MMDx) %>%
    arrange(CEL)


# SUMMARISE MMDx ####
baseline_K1208 %>%
    dplyr::select(Felzartamab, MMDx)

baseline_K1208 %>%
    tabyl(MMDx, Felzartamab) %>%
    janitor::adorn_totals()
