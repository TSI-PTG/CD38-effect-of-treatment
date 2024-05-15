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
    dplyr::rename(Center = Trial_Center, Patient = STUDY_EVALUATION_ID)  %>% 
    mutate_at(vars(contains("cfDNA")), ~ as.numeric(.) %>% suppressWarnings()) %>%
    pivot_longer(
        cols = contains("cfDNA"),
        names_to = "Followup",
        values_to = "cfDNA"
    ) %>%
    mutate(
        Patient = Patient %>% factor(),
        Followup = Followup %>% str_remove("GcfDNA_cp_ml_") %>%
            factor(levels = c("Day0", "Week12", "Week24", "Week52")),
        Group = case_when(
            Followup == "Day0" ~ "Index",
            Followup == "Week12" ~ "FU0",
            Followup == "Week24" ~ "FU1",
            Followup == "Week52" ~ "FU2",
        ) %>% factor(levels = c("Index", "FU0", "FU1", "FU2")),
        Felzartamab = Felzartamab %>% factor(labels = c("Placebo", "Felzartamab")),
        Felzartamab_Group = paste(Group, Felzartamab, sep = ":") %>%
            factor(levels = c(
                "Index:Placebo", "FU1b:Placebo", "FU1:Placebo", "FU2:Placebo",
                "Index:Felzartamab", "FU1b:Felzartamab", "FU1:Felzartamab", "FU2:Felzartamab"
            )),
        Felzartamab_Followup = paste(Followup, Felzartamab, sep = ":") %>%
            factor(levels = c(
                "Day0:Placebo", "Week12:Placebo", "Week24:Placebo", "Week52:Placebo",
                "Day0:Felzartamab", "Week12:Felzartamab", "Week24:Felzartamab", "Week52:Felzartamab"
            ))
    ) %>%
    dplyr::select(Center, Patient, Felzartamab, Group, Followup, Felzartamab_Group, Felzartamab_Followup, cfDNA)


# SAVE THE DATA ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/"
save(data_cfdna, file = paste(saveDir, "data_cfDNA.RData", sep = ""))
