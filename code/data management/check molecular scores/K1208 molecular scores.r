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
        ABMRpm, ggt0, TCMRt) %>%
    dplyr::rename(Trial_Center = Center,Felzartamab = Felzartamab_presumed)



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
        values_to = "ABMRpm_1208"
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
        values_to = "ggt0_1208"
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
        values_to = "TCMRt_1208"
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




# JOIN THE SPSS CELL AND ABMRpm DATA ####
data_spss <- reduce(
    list(data_CEL, data_ABMRpm, data_ggt0, data_TCMRt),
    left_join,
    by = c("Trial_Center", "STUDY_EVALUATION_ID", "Felzartamab", "Group", "Followup")
)


    

# JOIN THE SPSS REFERENCE SET DATA ####
data_joined <- data_spss %>%
    left_join(data_reference,
        by = c("CEL", "Trial_Center", "STUDY_EVALUATION_ID", "Felzartamab", "Group")
    )


data_joined %>%
    print(n="all")


data_joined %>%
    dplyr::select(contains("ABMR")) %>%
    arrange(ABMRpm) %>%
    data.frame




# SAVE THE DATA ####
save(data_histo, file = "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0Georg Felz CD38 Vienna/G_Rstuff/data/data_histo.RData")
