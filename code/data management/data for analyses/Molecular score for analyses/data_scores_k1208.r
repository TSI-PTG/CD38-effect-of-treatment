# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(haven) # install.packages("haven")
# Custom operators, functions, and datasets
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
source("C:/R/CD38-effect-of-treatment/code/functions/complex_pivot.R")
# load data
data <- read_spss("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0Georg Felz CD38 Vienna/G_Rstuff/data/Generalfile_Felzartamab SPSS.sav")
# load the processed cfDNA data
# load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/data_cfdna_cpml.RData")



# DEFINE PATIENT VARIABLES ####
vars_patient <- Hmisc::.q(
    Date_Tx,
    Female_gender,
    Age_at_Tx, Donor_Age,
    MM_ABDR, LD_Tx,
    Age_at_Screening, Years_to_ScreeningVisit, Screening_GFR, Screening_Prot_Krea,
    DSA_HLA_class_I_Only_Screening, DSA_HLA_class_II_Only_Screening, DSA_HLA_class_I_and_II_Screening, HLA_DQ_DSA_Screening,
    HLA_DQ_DSA_Screening, DSA_Number_Screening
)
data %>% dplyr::select(contains("mvi"))


# DEFINE VARIABLES ####
vars <- c(
    # "GcfDNA_cp_ml",
    # "Bx_date",
    "_RejPC1_1208Set",
    "_RejPC2_1208Set",
    "_RejPC3_1208Set",
    "_InjPC1_5086Set",
    "_InjPC2_5086Set",
    "_InjPC3_5086Set",
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
    "Bx_Borderline_Banff",
    "Bx_MVI_Score"
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
data %>% dplyr::select(contains("cfdna"))



# WRANGLE THE CP/ML cfDNA FROM SPSS DATA ####
data_cfdna_cpml <- data %>%
    dplyr::select(
        Trial_Center,
        STUDY_EVALUATION_ID,
        Felzartamab,
        contains("GcfDNA_cp_ml")
    ) %>%
    mutate_at(vars(contains("cfDNA")), ~ as.numeric(.) %>% suppressWarnings()) %>%
    pivot_longer(
        cols = contains("cfDNA"),
        names_to = "Group",
        values_to = "cfDNA_cpml"
    ) %>%
    mutate(
        Group = Group %>%
            str_remove("GcfDNA_cp_ml_") %>%
            str_replace("Day0", "Index") %>%
            str_replace("Week12", "FU1b") %>%
            str_replace("Week24", "FU1") %>%
            str_replace("Week52", "FU2")
    )


# WRANGLE THE PERCENT cfDNA FROM SPSS DATA ####
data_cfdna_percent <- data %>%
    dplyr::select(
        Trial_Center,
        STUDY_EVALUATION_ID,
        Felzartamab,
        contains("Percent_GcfDNA")
    ) %>%
    mutate_at(vars(contains("cfDNA")), ~ as.numeric(.) %>% suppressWarnings()) %>%
    pivot_longer(
        cols = contains("cfDNA"),
        names_to = "Group",
        values_to = "cfDNA_percent"
    ) %>%
    mutate(
        Group = Group %>%
            str_remove("Percent_GcfDNA_") %>%
            str_replace("Day0", "Index") %>%
            str_replace("Week12", "FU1b") %>%
            str_replace("Week24", "FU1") %>%
            str_replace("Week52", "FU2")
    )




# WRANGLE THE PATIENT DATA ####
data_patient <- data %>%
    dplyr::select(
        Trial_Center,
        STUDY_EVALUATION_ID,
        Felzartamab,
        all_of(vars_patient)
    )
data_patient %>% colnames()


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
    ) %>%
    mutate(Group = Group %>% str_remove("Bx"))


# WRANGLE THE BIOPSY DATES FROM SPSS DATA ####
data_dateBx <- data %>%
    dplyr::select(
        IndexBx_Date, FU1Bx_Date, FU1b_Bx_Date, FU2Bx_Date,
        Trial_Center,
        STUDY_EVALUATION_ID,
        Felzartamab
    ) %>%
    mutate(FU1b_Bx_Date = NA %>% as.Date()) %>%
    pivot_longer(
        cols = contains("Date"),
        names_to = "Group",
        names_pattern = "^(\\w+)_.*",
        values_to = "dateBx"
    ) %>%
    mutate(Group = Group %>% str_remove("Bx") %>% str_remove("_"))



# JOIN THE SPSS AND MOLECULAR SCORE DATA ####
data_k1208 <- data_scores %>%
    left_join(data_patient, by = c("Trial_Center", "STUDY_EVALUATION_ID", "Felzartamab")) %>%
    left_join(data_CEL, by = c("Trial_Center", "STUDY_EVALUATION_ID", "Felzartamab", "Group")) %>%
    left_join(data_dateBx, by = c("Trial_Center", "STUDY_EVALUATION_ID", "Felzartamab", "Group")) %>%
    left_join(data_cfdna_cpml, by = c("Trial_Center", "STUDY_EVALUATION_ID", "Felzartamab", "Group")) %>%
    left_join(data_cfdna_percent, by = c("Trial_Center", "STUDY_EVALUATION_ID", "Felzartamab", "Group")) %>%
    dplyr::rename(Center = Trial_Center, Patient = STUDY_EVALUATION_ID) %>%
    mutate(
        Patient = Patient %>% factor(),
        Group = Group %>%
            factor(levels = c("Index", "FU1b", "FU1", "FU2")),
        Followup = case_when(
            Group == "Index" ~ "Baseline",
            Group == "FU1b" ~ "Week12",
            Group == "FU1" ~ "Week24",
            Group == "FU2" ~ "Week52"
        ) %>% factor(levels = c("Baseline", "Week12", "Week24", "Week52")),
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
    relocate(
        Center, Patient, CEL, dateBx, Felzartamab, Group, Followup, Felzartamab_Group, Felzartamab_Followup,
        all_of(vars_patient),
        .before = 1
    ) %>%
    arrange(Felzartamab, Patient, Group) %>%
    relocate(cfDNA_percent, cfDNA_cpml, .before = "ABMRpm")


# RENAME THE MOLECULAR DATA ####
data_scores_k1208 <- data_k1208
data_scores_k1208  %>% colnames



# SAVE THE DATA ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/"
save(data_scores_k1208, file = paste(saveDir, "data_scores_k1208.RData", sep = ""))
