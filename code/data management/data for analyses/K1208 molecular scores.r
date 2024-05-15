# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(readr) # install.packages("readr")
library(haven) # install.packages("haven")
# Custom operators, functions, and datasets
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
source("C:/R/CD38-effect-of-treatment/code/functions/complex_pivot.R")
# load data
data <- read_spss("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0Georg Felz CD38 Vienna/G_Rstuff/data/Generalfile_Felzartamab SPSS.sav")
# load reference set
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0Georg Felz CD38 Vienna/G_Rstuff/data/vienna_1208_6Mar24.RData")


# TODO find a way to program the pivot_longer across multiple distinct variables at the same time


# DEFINE PATIENT VARIABLES ####
vars_patient <- Hmisc::.q(
    Female_gender, Age_at_Tx, LD_Tx, Donor_Age, MM_ABDR,
    Age_at_Screening, Years_to_ScreeningVisit, Screening_GFR, Screening_Prot_Krea,
    DSA_HLA_class_I_Only_Screening, DSA_HLA_class_II_Only_Screening, DSA_HLA_class_I_and_II_Screening, HLA_DQ_DSA_Screening,
    # mfi_immunodominant,
    HLA_DQ_DSA_Screening, DSA_Number_Screening
)
# data %>% dplyr::select(contains("mfi"))


# DEFINE VARIABLES ####
vars <- c(
    "_ABMRpm_1208Set",
    "_ggt0_1208Set",
    "_ptcgt0_1208Set",
    "_NKB_1208Set",
    "_DSAST_1208Set",
    "_TCMRt_1208Set",
    "_igt1_1208Set",
    "_tgt1_1208Set",
    "_TCB_1208Set",
    "_QCAT_1208Set",
    "_AMAT1_1208Set",
    "_QCMAT_1208Set",
    "_BAT_1208Set",
    "_cigt1_1208Set",
    "_ctgt1_1208Set",
    "_IGT_1208Set",
    "_MCAT_1208Set",
    "_ENDAT_1208Set",
    "_KT1_1208Set",
    "_KT2_1208Set",
    "_FICOL_1208Set",
    "_IRRAT30_1208Set",
    "_IRITD3_1208Set",
    "_IRITD5_1208Set",
    "_GRIT3_1208Set",
    "_Rej_RAT_1208Set",
    "_RejAA_NR_1208Set",
    "_RejAA_TCMR1_1208Set",
    "_RejAA_Mixed_1208Set",
    "_RejAA_EABMR_1208Set",
    "_RejAA_FABMR_1208Set",
    "_RejAA_LABMR_1208Set",
    "ABMRActivity_Banff19",
    "Bx_Active_ABMR_Banff19",
    "Bx_Chronic_active_ABMR_Banff19",
    "Bx_Borderline_Banff"
)






# WRANGLE THE PATIENT DATA #### %>%
data_patient <- data %>%
    dplyr::select(
        Trial_Center,
        STUDY_EVALUATION_ID,
        Felzartamab,
        all_of(vars_patient)
    )


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
        names_to = "Group",
        names_pattern = "^(\\w+)_.*",
        values_to = "CEL"
    )


# WRANGLE THE SPSS DATA ####
data_scores <- data %>% complex_pivot(
    target = vars,
    vars = c(
        "Trial_Center",
        "STUDY_EVALUATION_ID",
        "Felzartamab"
    )
)

data_scores %>% print(n = "all")



# JOIN THE SPSS AND MOLECULAR SCORE DATA ####
data_K1208 <- data_scores %>%
    left_join(data_patient, by = c("Trial_Center", "STUDY_EVALUATION_ID", "Felzartamab")) %>%
    dplyr::rename(Center = Trial_Center, Patient = STUDY_EVALUATION_ID) %>%
    mutate(
        Patient = Patient %>% factor(),
        Followup = Followup %>%
            factor(levels = c("Day0", "Week12", "Week24", "Week52")),
        Group = case_when(
            Followup == "Day0" ~ "Index",
            Followup == "Week12" ~ "FU0",
            Followup == "Week24" ~ "FU1",
            Followup == "Week52" ~ "FU2",
        ) %>% factor(levels = c("Index", "FU0", "FU1", "FU2")),
        Felzartamab = Felzartamab %>% factor(labels = c("Placebo", "Felzartamab")),
        Group_Felzartamab = paste(Group, Felzartamab, sep = ":") %>%
            factor(levels = c(
                "Index:Placebo", "FU0:Placebo", "FU1:Placebo", "FU2:Placebo",
                "Index:Felzartamab", "FU0:Felzartamab", "FU1:Felzartamab", "FU2:Felzartamab"
            )),
        Group_Followup = paste(Followup, Felzartamab, sep = ":") %>%
            factor(levels = c(
                "Day0:Placebo", "Week12:Placebo", "Week24:Placebo", "Week52:Placebo",
                "Day0:Felzartamab", "Week12:Felzartamab", "Week24:Felzartamab", "Week52:Felzartamab"
            ))
    ) %>%
    relocate(
        Center, Patient, Felzartamab, Group, Followup, Group_Felzartamab, Group_Followup,
        all_of(vars_patient),
        .before = 1
    )
data_K1208$Followup




# # SAVE THE DATA ####
# saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/"
# save(data_K1208, file = paste(saveDir, "data_K1208.RData", sep = ""))
