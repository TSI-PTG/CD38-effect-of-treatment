# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(readr) # install.packages("readr")
library(haven) # install.packages("haven")
# Bioconductor libraries
library(Biobase) # BiocManager::install("Biobase")
# Custom operators, functions, and datasets
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
# load data
data <- read_spss("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0Georg Felz CD38 Vienna/G_Rstuff/data/Generalfile_Felzartamab SPSS.sav")
# load reference set
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0Georg Felz CD38 Vienna/G_Rstuff/data/vienna_1208_6Mar24.RData")





# TODO: find a way to program the pivot_longer across multiple distinct variables at the same time


# DEFINE THE PATIENT VARIABLES ####
vars_patient <- Hmisc::.q(
    Female_gender, Age_at_Tx, Donor_Age, Screening_GFR
)


# DEFINE HISTOLOGY VARIABLES ####
vars_histology <- Hmisc::.q(
    IndexBx_v, FU1Bx_v, FU2Bx_v,
    IndexBx_g, FU1Bx_g, FU2Bx_g,
    IndexBx_ptc, FU1Bx_ptc, FU2Bx_ptc,
    IndexBx_i, FU1Bx_i, FU2Bx_i,
    IndexBx_t, FU1Bx_t, FU2Bx_t,
    IndexBx_cg, FU1Bx_cg, FU2Bx_cg
)


# DEFINE THE MOLECULAR VARIABLES####
vars_molecular <- Hmisc::.q(
    ABMRpm, ggt0, ptcgt0, NKB, DSAST,
    TCMRt, tgt1, igt1, TCB, QCAT
)


# DEFINE THE ABMR ACTIVITY GENES ####
vars_molecular <- Hmisc::.q(
    
)




# WRANGLE THE PATIENT DATA ####
data_patient <- data %>%
    dplyr::select(
        Trial_Center,
        STUDY_EVALUATION_ID,
        Felzartamab,
        all_of(vars_patient)
    )


# WRANGLE THE REFERENCE SET DATA ####
data_reference <- vienna_1208 %>%
    pData() %>%
    tibble() %>%
    dplyr::select(
        CEL, Center, STUDY_EVALUATION_ID, Felzartamab_presumed, Group,
        all_of(vars_molecular)
    ) %>%
    dplyr::rename(
        Trial_Center = Center,
        Felzartamab = Felzartamab_presumed
    )


# WRANGLE THE CEL IDS FROM SPSS DATA ####
data_CEL <- data %>%
    dplyr::select(
        Index_CEL, FU1_CEL, FU2_CEL,
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
    )

data_dateBx <- data %>%
    dplyr::select(
        Trial_Center,
        STUDY_EVALUATION_ID,
        IndexBx_Date, FU1Bx_Date, FU2Bx_Date,
        Felzartamab
    ) %>%
    pivot_longer(
        cols = contains("Date"),
        names_to = c("Group"),
        names_pattern = "^(\\w+)_.*",
        values_to = "dateBx"
    ) %>%
    mutate(
        Group = Group %>%
            str_remove_all("Bx"),
        Followup = case_when(
            Group == "Index" ~ "Day0",
            Group == "FU0" ~ "Week12",
            Group == "FU1" ~ "Week24",
            Group == "FU2" ~ "Week52",
        )
    )


# WRANGLE THE g SCORES FROM SPSS DATA ####
data_g <- data %>%
    dplyr::select(
        IndexBx_g, FU1Bx_g, FU2Bx_g,
        Trial_Center,
        STUDY_EVALUATION_ID,
        Felzartamab
    ) %>%
    pivot_longer(
        cols = c("IndexBx_g", "FU1Bx_g", "FU2Bx_g"),
        names_to = c("Group"),
        names_pattern = "^(\\w+)_.*",
        values_to = "g"
    ) %>%
    mutate(
        Group = Group %>%
            str_remove_all("Bx"),
        Followup = case_when(
            Group == "Index" ~ "Day0",
            Group == "FU0" ~ "Week12",
            Group == "FU1" ~ "Week24",
            Group == "FU2" ~ "Week52",
        )
    )

# WRANGLE THE cg SCORES FROM SPSS DATA ####
data_cg <- data %>%
    dplyr::select(
        IndexBx_cg, FU1Bx_cg, FU2Bx_cg,
        Trial_Center,
        STUDY_EVALUATION_ID,
        Felzartamab
    ) %>%
    pivot_longer(
        cols = c("IndexBx_cg", "FU1Bx_cg", "FU2Bx_cg"),
        names_to = c("Group"),
        names_pattern = "^(\\w+)_.*",
        values_to = "cg"
    ) %>%
    mutate(
        Group = Group %>%
            str_remove_all("Bx"),
        Followup = case_when(
            Group == "Index" ~ "Day0",
            Group == "FU0" ~ "Week12",
            Group == "FU1" ~ "Week24",
            Group == "FU2" ~ "Week52",
        )
    )


# WRANGLE THE ptc SCORES FROM SPSS DATA ####
data_ptc <- data %>%
    dplyr::select(
        IndexBx_ptc, FU1Bx_ptc, FU2Bx_ptc,
        Trial_Center,
        STUDY_EVALUATION_ID,
        Felzartamab
    ) %>%
    # mutate_at(vars(contains("ABMRpm_1208Set")), ~ as.numeric(.) %>% suppressWarnings()) %>%
    pivot_longer(
        cols = c("IndexBx_ptc", "FU1Bx_ptc", "FU2Bx_ptc"),
        names_to = c("Group"),
        names_pattern = "^(\\w+)_.*",
        values_to = "ptc"
    ) %>%
    mutate(
        Group = Group %>%
            str_remove_all("Bx"),
        Followup = case_when(
            Group == "Index" ~ "Day0",
            Group == "FU0" ~ "Week12",
            Group == "FU1" ~ "Week24",
            Group == "FU2" ~ "Week52",
        )
    )


# WRANGLE THE v SCORES FROM SPSS DATA ####
data_v <- data %>%
    dplyr::select(
        IndexBx_v, FU1Bx_v, FU2Bx_v,
        Trial_Center,
        STUDY_EVALUATION_ID,
        Felzartamab
    ) %>%
    # mutate_at(vars(contains("ABMRpm_1208Set")), ~ as.numeric(.) %>% suppressWarnings()) %>%
    pivot_longer(
        cols = c("IndexBx_v", "FU1Bx_v", "FU2Bx_v"),
        names_to = c("Group"),
        names_pattern = "^(\\w+)_.*",
        values_to = "v"
    ) %>%
    mutate(
        Group = Group %>%
            str_remove_all("Bx"),
        Followup = case_when(
            Group == "Index" ~ "Day0",
            Group == "FU0" ~ "Week12",
            Group == "FU1" ~ "Week24",
            Group == "FU2" ~ "Week52",
        )
    )


# WRANGLE THE i SCORES FROM SPSS DATA ####
data_i <- data %>%
    dplyr::select(
        IndexBx_i, FU1Bx_i, FU2Bx_i,
        Trial_Center,
        STUDY_EVALUATION_ID,
        Felzartamab
    ) %>%
    pivot_longer(
        cols = c("IndexBx_i", "FU1Bx_i", "FU2Bx_i"),
        names_to = c("Group"),
        names_pattern = "^(\\w+)_.*",
        values_to = "i"
    ) %>%
    mutate(
        Group = Group %>%
            str_remove_all("Bx"),
        Followup = case_when(
            Group == "Index" ~ "Day0",
            Group == "FU0" ~ "Week12",
            Group == "FU1" ~ "Week24",
            Group == "FU2" ~ "Week52",
        )
    )


# WRANGLE THE t SCORES FROM SPSS DATA ####
data_t <- data %>%
    dplyr::select(
        IndexBx_t, FU1Bx_t, FU2Bx_t,
        Trial_Center,
        STUDY_EVALUATION_ID,
        Felzartamab
    ) %>%
    pivot_longer(
        cols = c("IndexBx_t", "FU1Bx_t", "FU2Bx_t"),
        names_to = c("Group"),
        names_pattern = "^(\\w+)_.*",
        values_to = "t"
    ) %>%
    mutate(
        Group = Group %>%
            str_remove_all("Bx"),
        Followup = case_when(
            Group == "Index" ~ "Day0",
            Group == "FU0" ~ "Week12",
            Group == "FU1" ~ "Week24",
            Group == "FU2" ~ "Week52",
        )
    )


# WRANGLE ABMRActivity_Banff19 FROM SPSS DATA ####
data_ABMRActivity_Banff19 <- data %>%
    dplyr::select(
        contains("ABMRActivity_Banff19"),
        Trial_Center,
        STUDY_EVALUATION_ID,
        Felzartamab
    ) %>%
    mutate_at(vars(contains("ABMRactivity"), ends_with("_g")), ~ as.numeric(.) %>% suppressWarnings()) %>%
    pivot_longer(
        cols = contains("ABMRactivity"),
        names_to = c("Group"),
        values_to = "ABMRActivity_Banff19"
    ) %>%
    mutate(
        Group = Group %>%
            str_remove_all("Bx_Morphologic_ABMRactivity_Banff19") %>%
            str_remove_all("Bx_Morphologic_ABMRActivity_Banff19"),
        Followup = case_when(
            Group == "Index" ~ "Day0",
            Group == "FU0" ~ "Week12",
            Group == "FU1" ~ "Week24",
            Group == "FU2" ~ "Week52",
        )
    ) %>%
    dplyr::filter(Followup %>% str_detect("FU1b", negate = TRUE))


# WRANGLE Active_ABMR_Banff19 FROM SPSS DATA ####
data_Active_ABMR_Banff19 <- data %>%
    dplyr::select(
        contains("Bx_Active_ABMR_Banff19"),
        Trial_Center,
        STUDY_EVALUATION_ID,
        Felzartamab
    ) %>%
    mutate_at(vars(contains("Bx_Active_ABMR_Banff19")), ~ as.numeric(.) %>% suppressWarnings()) %>%
    pivot_longer(
        cols = contains("Bx_Active_ABMR_Banff19"),
        names_to = c("Group"),
        values_to = "Active_ABMR_Banff19"
    ) %>%
    mutate(
        Group = Group %>%
            str_remove_all("Bx_Active_ABMR_Banff19"),
        Followup = case_when(
            Group == "Index" ~ "Day0",
            Group == "FU0" ~ "Week12",
            Group == "FU1" ~ "Week24",
            Group == "FU2" ~ "Week52",
        )
    ) %>%
    dplyr::filter(Followup %>% str_detect("FU1b", negate = TRUE))


# WRANGLE Active_ABMR_Banff19 FROM SPSS DATA ####
data_Chronic_active_ABMR_Banff19 <- data %>%
    dplyr::select(
        contains("Bx_Chronic_active_ABMR_Banff19"),
        Trial_Center,
        STUDY_EVALUATION_ID,
        Felzartamab
    ) %>%
    mutate_at(vars(contains("Bx_Chronic_active_ABMR_Banff19")), ~ as.numeric(.) %>% suppressWarnings()) %>%
    pivot_longer(
        cols = contains("Bx_Chronic_active_ABMR_Banff19"),
        names_to = c("Group"),
        values_to = "Chronic_active_ABMR_Banff19"
    ) %>%
    mutate(
        Group = Group %>%
            str_remove_all("Bx_Chronic_active_ABMR_Banff19"),
        Followup = case_when(
            Group == "Index" ~ "Day0",
            Group == "FU0" ~ "Week12",
            Group == "FU1" ~ "Week24",
            Group == "FU2" ~ "Week52",
        )
    ) %>%
    dplyr::filter(Followup %>% str_detect("FU1b", negate = TRUE))


# WRANGLE Borderline_Banff FROM SPSS DATA ####
data_Borderline_Banff <- data %>%
    dplyr::select(
        contains("Bx_Borderline_Banff"),
        Trial_Center,
        STUDY_EVALUATION_ID,
        Felzartamab
    ) %>%
    mutate_at(vars(contains("Bx_Borderline_Banff")), ~ as.numeric(.) %>% suppressWarnings()) %>%
    pivot_longer(
        cols = contains("Bx_Borderline_Banff"),
        names_to = c("Group"),
        values_to = "Borderline_Banff"
    ) %>%
    mutate(
        Group = Group %>%
            str_remove_all("Bx_Borderline_Banff"),
        Followup = case_when(
            Group == "Index" ~ "Day0",
            Group == "FU0" ~ "Week12",
            Group == "FU1" ~ "Week24",
            Group == "FU2" ~ "Week52",
        )
    ) %>%
    dplyr::filter(Followup %>% str_detect("FU1b", negate = TRUE))


# JOIN THE SPSS AND MOLECULAR SCORE DATA ####
data_K1208 <- reduce(
    list(
        data_CEL,
        data_dateBx,
        data_g, data_cg, data_ptc, data_v, data_i, data_t,
        data_ABMRActivity_Banff19, data_Active_ABMR_Banff19, data_Chronic_active_ABMR_Banff19,
        data_Borderline_Banff
    ),
    left_join,
    by = c("Trial_Center", "STUDY_EVALUATION_ID", "Felzartamab", "Group", "Followup")
) %>%
    left_join(data_patient, by = c("Trial_Center", "STUDY_EVALUATION_ID", "Felzartamab")) %>%
    left_join(data_reference, by = c("CEL", "Trial_Center", "STUDY_EVALUATION_ID", "Felzartamab", "Group")) %>%
    mutate(
        Trial_Center = Trial_Center %>% as.character(),
        STUDY_EVALUATION_ID = STUDY_EVALUATION_ID %>% as.numeric(),
        Felzartamab = Felzartamab %>% as.numeric(),
        Group = Group %>% factor(levels = c("Index", "FU1", "FU2"))
    ) %>%
    relocate(CEL, names(data_patient), .before = 1) %>%
    relocate(Group, dateBx, .after = "Felzartamab") %>%
    dplyr::rename(patient_ID = STUDY_EVALUATION_ID) %>%
    mutate(
        Felzartamab = Felzartamab %>%
            factor(
                labels = c("placebo", "felzartamab")
            )
    ) %>%
    dplyr::filter(patient_ID %nin% c(15, 18))




# SAVE THE DATA ####
save(data_K1208, file = "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0Georg Felz CD38 Vienna/G_Rstuff/data/data_K1208.RData")



# # JOIN THE SPSS REFERENCE SET DATA ####
# data_joined <- data_spss %>%
#     left_join(data_reference,
#         by = c("CEL", "Trial_Center", "STUDY_EVALUATION_ID", "Felzartamab", "Group")
#     )


# data_joined %>%
#     print(n="all")
