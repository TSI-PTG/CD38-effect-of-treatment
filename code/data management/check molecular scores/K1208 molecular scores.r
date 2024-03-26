# HOUSEKEEPING ####
library(tidyverse) # install.packages("tidyverse")
library(readr) # install.packages("readr")
library(haven) # install.packages("haven")
# Custom operators, functions, and datasets
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
# load data
data <- read_spss("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0Georg Felz CD38 Vienna/G_Rstuff/data/Generalfile_Felzartamab SPSS.sav")
# load reference set
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0Georg Felz CD38 Vienna/G_Rstuff/data/vienna_1208_6Mar24.RData")



# WRANGLE THE REFERENCE SET DATA ####
data_reference <- vienna_1208 %>%
    pData() %>%
    tibble() %>%
    dplyr::select(
        Center, STUDY_EVALUATION_ID, Felzartamab_presumed, CEL, Group,
        ABMRpm, ggt0, ptcgt0, NKB, DSAST,
        TCMRt, tgt1, igt1, TCB, QCAT
    ) %>%
    dplyr::rename(Trial_Center = Center, Felzartamab = Felzartamab_presumed)



# WRANGLE THE CEL IDS FROM SPSS DATA ####
data_CEL <- data %>%
    dplyr::select(
        Index_CEL, FU1_CEL, FU1bBx_CEL, FU2_CEL,
        Trial_Center,
        STUDY_EVALUATION_ID,
        Felzartamab
    ) %>%
    pivot_longer(
        cols = contains("CEL"),
        names_to = c("Group"),
        names_pattern = "^(\\w+)_.*",
        values_to = "CEL"
    ) %>%
    mutate(
        Group = Group %>%
            str_remove_all("_ABMRpm"),
        Followup = case_when(
            Group == "Index" ~ "Day0",
            Group == "FU0" ~ "Week12",
            Group == "FU1" ~ "Week24",
            Group == "FU2" ~ "Week52",
        )
    ) %>%
    dplyr::filter(Followup %>% str_detect("FU1b", negate = TRUE))


# WRANGLE THE ABMRpm SCORES FROM SPSS DATA ####
data_ABMRpm <- data %>%
    dplyr::select(
        contains("ABMRpm_1208Set"),
        Trial_Center,
        STUDY_EVALUATION_ID,
        Felzartamab
    ) %>%
    mutate_at(vars(contains("ABMRpm_1208Set")), ~ as.numeric(.) %>% suppressWarnings()) %>%
    pivot_longer(
        cols = contains("ABMRpm_1208Set"),
        names_to = c("Group"),
        names_pattern = "^(\\w+)_.*",
        values_to = "ABMRpm"
    ) %>%
    mutate(
        Group = Group %>%
            str_remove_all("_ABMRpm"),
        Followup = case_when(
            Group == "Index" ~ "Day0",
            Group == "FU0" ~ "Week12",
            Group == "FU1" ~ "Week24",
            Group == "FU2" ~ "Week52",
        )
    ) %>%
    dplyr::filter(Followup %>% str_detect("FU1b", negate = TRUE))


# WRANGLE THE ggt0 SCORES FROM SPSS DATA ####
data_ggt0 <- data %>%
    dplyr::select(
        contains("ggt0_1208Set"),
        Trial_Center,
        STUDY_EVALUATION_ID,
        Felzartamab
    ) %>%
    mutate_at(vars(contains("Aggt0_1208Set")), ~ as.numeric(.) %>% suppressWarnings()) %>%
    pivot_longer(
        cols = contains("ggt0_1208Set"),
        names_to = c("Group"),
        names_pattern = "^(\\w+)_.*",
        values_to = "ggt0"
    ) %>%
    mutate(
        Group = Group %>%
            str_remove_all("_ggt0"),
        Followup = case_when(
            Group == "Index" ~ "Day0",
            Group == "FU0" ~ "Week12",
            Group == "FU1" ~ "Week24",
            Group == "FU2" ~ "Week52",
        )
    ) %>%
    dplyr::filter(Followup %>% str_detect("FU1b", negate = TRUE))


# WRANGLE THE ptcgt0 SCORES FROM SPSS DATA ####
data_ptcgt0 <- data %>%
    dplyr::select(
        contains("ptcgt0_1208Set"),
        Trial_Center,
        STUDY_EVALUATION_ID,
        Felzartamab
    ) %>%
    mutate_at(vars(contains("ptcgt0_1208Set")), ~ as.numeric(.) %>% suppressWarnings()) %>%
    pivot_longer(
        cols = contains("ptcgt0_1208Set"),
        names_to = c("Group"),
        names_pattern = "^(\\w+)_.*",
        values_to = "ptcgt0"
    ) %>%
    mutate(
        Group = Group %>%
            str_remove_all("_ptcgt0"),
        Followup = case_when(
            Group == "Index" ~ "Day0",
            Group == "FU0" ~ "Week12",
            Group == "FU1" ~ "Week24",
            Group == "FU2" ~ "Week52",
        )
    ) %>%
    dplyr::filter(Followup %>% str_detect("FU1b", negate = TRUE))


# WRANGLE THE NKB SCORES FROM SPSS DATA ####
data_NKB <- data %>%
    dplyr::select(
        contains("NKB_1208Set"),
        Trial_Center,
        STUDY_EVALUATION_ID,
        Felzartamab
    ) %>%
    mutate_at(vars(contains("NKB_1208Set")), ~ as.numeric(.) %>% suppressWarnings()) %>%
    pivot_longer(
        cols = contains("NKB_1208Set"),
        names_to = c("Group"),
        names_pattern = "^(\\w+)_.*",
        values_to = "NKB"
    ) %>%
    mutate(
        Group = Group %>%
            str_remove_all("_NKB"),
        Followup = case_when(
            Group == "Index" ~ "Day0",
            Group == "FU0" ~ "Week12",
            Group == "FU1" ~ "Week24",
            Group == "FU2" ~ "Week52",
        )
    ) %>%
    dplyr::filter(Followup %>% str_detect("FU1b", negate = TRUE))


# WRANGLE THE DSAST SCORES FROM SPSS DATA ####
data_DSAST <- data %>%
    dplyr::select(
        contains("DSAST_1208Set"),
        Trial_Center,
        STUDY_EVALUATION_ID,
        Felzartamab
    ) %>%
    mutate_at(vars(contains("DSAST_1208Set")), ~ as.numeric(.) %>% suppressWarnings()) %>%
    pivot_longer(
        cols = contains("DSAST_1208Set"),
        names_to = c("Group"),
        names_pattern = "^(\\w+)_.*",
        values_to = "DSAST"
    ) %>%
    mutate(
        Group = Group %>%
            str_remove_all("_DSAST"),
        Followup = case_when(
            Group == "Index" ~ "Day0",
            Group == "FU0" ~ "Week12",
            Group == "FU1" ~ "Week24",
            Group == "FU2" ~ "Week52",
        )
    ) %>%
    dplyr::filter(Followup %>% str_detect("FU1b", negate = TRUE))


# WRANGLE THE TCMRt SCORES FROM SPSS DATA ####
data_TCMRt <- data %>%
    dplyr::select(
        contains("TCMRt_1208Set"),
        Trial_Center,
        STUDY_EVALUATION_ID,
        Felzartamab
    ) %>%
    mutate_at(vars(contains("TCMRt_1208Set")), ~ as.numeric(.) %>% suppressWarnings()) %>%
    pivot_longer(
        cols = contains("TCMRt_1208Set"),
        names_to = c("Group"),
        names_pattern = "^(\\w+)_.*",
        values_to = "TCMRt"
    ) %>%
    mutate(
        Group = Group %>%
            str_remove_all("_TCMRt"),
        Followup = case_when(
            Group == "Index" ~ "Day0",
            Group == "FU0" ~ "Week12",
            Group == "FU1" ~ "Week24",
            Group == "FU2" ~ "Week52",
        )
    ) %>%
    dplyr::filter(Followup %>% str_detect("FU1b", negate = TRUE))


# WRANGLE THE igt1 SCORES FROM SPSS DATA ####
data_igt1 <- data %>%
    dplyr::select(
        contains("igt1_1208Set"),
        Trial_Center,
        STUDY_EVALUATION_ID,
        Felzartamab
    ) %>%
    mutate_at(vars(contains("igt1_1208Set")), ~ as.numeric(.) %>% suppressWarnings()) %>%
    pivot_longer(
        cols = contains("igt1_1208Set"),
        names_to = c("Group"),
        names_pattern = "^(\\w+)_.*",
        values_to = "igt1"
    ) %>%
    mutate(
        Group = Group %>%
            str_remove_all("_igt1"),
        Followup = case_when(
            Group == "Index" ~ "Day0",
            Group == "FU0" ~ "Week12",
            Group == "FU1" ~ "Week24",
            Group == "FU2" ~ "Week52",
        )
    ) %>%
    dplyr::filter(Followup %>% str_detect("FU1b", negate = TRUE))


# WRANGLE THE tgt1 SCORES FROM SPSS DATA ####
data_tgt1 <- data %>%
    dplyr::select(
        contains("tgt1_1208Set"),
        Trial_Center,
        STUDY_EVALUATION_ID,
        Felzartamab
    ) %>%
    mutate_at(vars(contains("tgt1_1208Set")), ~ as.numeric(.) %>% suppressWarnings()) %>%
    pivot_longer(
        cols = contains("tgt1_1208Set"),
        names_to = c("Group"),
        names_pattern = "^(\\w+)_.*",
        values_to = "tgt1"
    ) %>%
    mutate(
        Group = Group %>%
            str_remove_all("_tgt1"),
        Followup = case_when(
            Group == "Index" ~ "Day0",
            Group == "FU0" ~ "Week12",
            Group == "FU1" ~ "Week24",
            Group == "FU2" ~ "Week52",
        )
    ) %>%
    dplyr::filter(Followup %>% str_detect("FU1b", negate = TRUE))


# WRANGLE THE TCB SCORES FROM SPSS DATA ####
data_TCB <- data %>%
    dplyr::select(
        contains("TCB_1208Set"),
        Trial_Center,
        STUDY_EVALUATION_ID,
        Felzartamab
    ) %>%
    mutate_at(vars(contains("TCB_1208Set")), ~ as.numeric(.) %>% suppressWarnings()) %>%
    pivot_longer(
        cols = contains("TCB_1208Set"),
        names_to = c("Group"),
        names_pattern = "^(\\w+)_.*",
        values_to = "TCB"
    ) %>%
    mutate(
        Group = Group %>%
            str_remove_all("_TCB"),
        Followup = case_when(
            Group == "Index" ~ "Day0",
            Group == "FU0" ~ "Week12",
            Group == "FU1" ~ "Week24",
            Group == "FU2" ~ "Week52",
        )
    ) %>%
    dplyr::filter(Followup %>% str_detect("FU1b", negate = TRUE))


# WRANGLE THE QCAT SCORES FROM SPSS DATA ####
data_QCAT <- data %>%
    dplyr::select(
        contains("QCAT_1208Set"),
        Trial_Center,
        STUDY_EVALUATION_ID,
        Felzartamab
    ) %>%
    mutate_at(vars(contains("QCAT_1208Set")), ~ as.numeric(.) %>% suppressWarnings()) %>%
    pivot_longer(
        cols = contains("QCAT_1208Set"),
        names_to = c("Group"),
        names_pattern = "^(\\w+)_.*",
        values_to = "QCAT"
    ) %>%
    mutate(
        Group = Group %>%
            str_remove_all("_QCAT"),
        Followup = case_when(
            Group == "Index" ~ "Day0",
            Group == "FU0" ~ "Week12",
            Group == "FU1" ~ "Week24",
            Group == "FU2" ~ "Week52",
        )
    ) %>%
    dplyr::filter(Followup %>% str_detect("FU1b", negate = TRUE))


# JOIN THE SPSS CELL AND MOLECULAR SCORE DATA ####
data_K1208 <- reduce(
    list(
        data_CEL,
        data_ABMRpm, data_ggt0, data_ptcgt0, data_NKB, data_DSAST,
        data_TCMRt, data_tgt1, data_igt1, data_TCB, data_QCAT
    ),
    left_join,
    by = c("Trial_Center", "STUDY_EVALUATION_ID", "Felzartamab", "Group", "Followup")
) %>%
    mutate(
        Trial_Center = Trial_Center %>% as.character(),
        STUDY_EVALUATION_ID = STUDY_EVALUATION_ID %>% as.character(),
        Felzartamab = Felzartamab %>% as.character()
    ) 



# SAVE THE DATA ####
save(data_K1208, file = "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0Georg Felz CD38 Vienna/G_Rstuff/data/data_K1208.RData")







# # JOIN THE SPSS REFERENCE SET DATA ####
# data_joined <- data_spss %>%
#     left_join(data_reference,
#         by = c("CEL", "Trial_Center", "STUDY_EVALUATION_ID", "Felzartamab", "Group")
#     )


# data_joined %>%
#     print(n="all")
