# HOUSEKEEPING ####
library(tidyverse)
library(readr)
# Custom operators, functions, and datasets
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
# load data
data <- read_csv("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0Georg Felz CD38 Vienna/G_Rstuff/data/felzartamab_study_data.csv")


# WRANGLE THE DATA ####
data_cfdna <- data %>%
    dplyr::select(
        Trial_Center, STUDY_EVALUATION_ID,
        Felzartamab,  contains("GcfDNA_cp_ml")
    ) %>%
    mutate_at(vars(contains("cfDNA")), ~ as.numeric(.)%>% suppressWarnings()) %>%
    pivot_longer(
        cols = contains("GcfDNA_cp_ml"),
        names_to = "Followup", values_to = "cfDNA"
    ) %>%
    mutate(
        Followup = Followup %>% str_remove("GcfDNA_cp_ml_"),
        Group = case_when(
            Followup == "Day0" ~ "Index",
            Followup == "Week12" ~ "FU0",
            Followup == "Week24" ~ "FU1",
            Followup == "Week52" ~ "FU2",
        )
    )
    

# SAVE THE DATA ####
save(data_cfdna, file = "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0Georg Felz CD38 Vienna/G_Rstuff/data/data_cfDNA.RData")
