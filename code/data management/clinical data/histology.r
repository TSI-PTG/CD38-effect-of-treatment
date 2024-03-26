# HOUSEKEEPING ####
library(tidyverse)
library(readr)
library(haven)
# Custom operators, functions, and datasets
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
# load data
data <- read_csv("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0Georg Felz CD38 Vienna/G_Rstuff/data/felzartamab_study_data.csv")





# VARIABLES OF INTEREST ####
vars <- c(
    # contains("ABMRactivity"),
    "Trial_Center",
    "STUDY_EVALUATION_ID",
    "Felzartamab",
    "Age_at_Tx",
    "Donor_Age"
    #     "FU1Bx_i",
    #     "FU1Bx_ti",
    #     "FU1Bx_t",
    #     "FU1Bx_v",
    #     "FU1Bx_ah",
    #     "FU1Bx_cg",
    #     "FU1Bx_ci",
    #     "FU1Bx_ct",
    #     "FU1Bx_cv",
    #     "FU1Bx_mm",
    #     "FU1Bx_ptc",
    #     "FU1Bx_C4d"
)


# WRANGLE THE DATA ####
data_histo <- data %>%
    dplyr::select(
        # contains("ABMRactivity"),
        ends_with("CEL"),
        # ends_with("_g"),

        # all_of(vars)
    ) %>%
    mutate_at(vars(contains("ABMRactivity"), ends_with("_g")), ~ as.numeric(.) %>% suppressWarnings()) %>%
    pivot_longer(
        cols = contains("ABMRactivity"),
        names_to = c("Group"), 
        names_pattern = "^(\\w+)_.*",
        values_to = "Bx_Morphologic_ABMRActivity_Banff19"
    ) %>%
    pivot_longer(
        cols = ends_with("CEL"),
        names_to = NULL, values_to = "CEL"
    ) %>%
    distinct(CEL, .keep_all = TRUE) %>%
    mutate(
        Group = Group %>%
            str_remove_all("Bx_Morphologic_ABMRActivity_Banff19"),
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
