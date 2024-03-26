# HOUSEKEEPING ####
library(tidyverse)
library(readr)
# Custom operators, functions, and datasets
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
# load data
data <- read_csv("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0Georg Felz CD38 Vienna/G_Rstuff/data/felzartamab_study_data.csv")


# WRANGLE THE DATA ####
data_histo <- data %>%
    dplyr::select(
        Trial_Center, STUDY_EVALUATION_ID,
        Felzartamab,  contains("ABMRactivity")
    ) %>%
    mutate_at(vars(contains("ABMRactivity")), ~ as.numeric(.) %>% suppressWarnings()) %>%
    pivot_longer(
        cols = contains("ABMRactivity"),
        names_to = "Group", values_to = "ABMRactivity"
    ) %>%
    mutate(
        Group = Group %>%
            str_remove_all("Bx_Morphologic_ABMRActivity_Banff19") %>%
            str_remove_all("Bx_Morphologic_ABMRactivity_Banff19"),
        Followup = case_when(
            Group == "Index" ~ "Day0",
            Group == "FU0" ~ "Week12",
            Group == "FU1" ~ "Week24",
            Group == "FU2" ~ "Week52",
        )
    ) %>%
    dplyr::filter(Followup %>% str_detect("FU1b", negate = TRUE))


# SAVE THE DATA ####
save(data_histo, file = "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0Georg Felz CD38 Vienna/G_Rstuff/data/data_histo.RData")



data_histo %>% dplyr::filter(STUDY_EVALUATION_ID %in% c(2, 9, 13, 19))
