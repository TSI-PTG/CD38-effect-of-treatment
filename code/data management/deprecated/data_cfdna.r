# HOUSEKEEPING ####
library(tidyverse) # install.packages("tidyverse")
library(haven) # install.packages("haven")
# Custom operators, functions, and datasets
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
# # load data
data <- read_spss("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0Georg Felz CD38 Vienna/G_Rstuff/data/Generalfile_Felzartamab SPSS.sav")


# WRANGLE THE DATA ####
data_cfdna <- data %>%
    dplyr::select(
        Trial_Center,
        STUDY_EVALUATION_ID,
        Felzartamab,
        contains("GcfDNA_cp_ml")
    ) %>%
    dplyr::rename(Center = Trial_Center, Patient = STUDY_EVALUATION_ID) %>%
    mutate_at(vars(contains("cfDNA")), ~ as.numeric(.) %>% suppressWarnings()) %>%
    pivot_longer(
        cols = contains("cfDNA"),
        names_to = "Followup",
        values_to = "cfDNA"
    ) %>%
    mutate(
        Patient = Patient %>% factor(),
        Followup = Followup %>% str_remove("GcfDNA_cp_ml_"),
        Followup = ifelse(Followup == "Day0", "Baseline", Followup) %>%
            factor(levels = c("Baseline", "Week12", "Week24", "Week52")),
        Group = case_when(
            Followup == "Baseline" ~ "Index",
            Followup == "Week12" ~ "FU1b",
            Followup == "Week24" ~ "FU1",
            Followup == "Week52" ~ "FU2",
        ) %>% factor(levels = c("Index", "FU1b", "FU1", "FU2")),
        Felzartamab = Felzartamab %>% factor(labels = c("Placebo", "Felzartamab")),
        Felzartamab_Group = paste(Group, Felzartamab, sep = "_") %>%
            factor(levels = c(
                "Index_Placebo", "FU1b_Placebo", "FU1_Placebo", "FU2_Placebo",
                "Index_Felzartamab", "FU1b_Felzartamab", "FU1_Felzartamab", "FU2_Felzartamab"
            )),
        Felzartamab_Followup = paste(Followup, Felzartamab, sep = "_") %>%
            factor(levels = c(
                "Baseline_Placebo", "Week12_Placebo", "Week24_Placebo", "Week52_Placebo",
                "Baseline_Felzartamab", "Week12_Felzartamab", "Week24_Felzartamab", "Week52_Felzartamab"
            ))
    ) %>%
    dplyr::select(Center, Patient, Felzartamab, Group, Followup, Felzartamab_Group, Felzartamab_Followup, cfDNA) %>%
    arrange(Felzartamab, Patient, Followup)


# SAVE THE DATA ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/"
save(data_cfdna, file = paste(saveDir, "data_cfdna.RData", sep = ""))
