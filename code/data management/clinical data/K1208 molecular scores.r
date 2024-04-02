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


# TODO find a way to program the pivot_longer across multiple distinct variables at the same time


# DEFINE PATIENT VARIABLES ####
vars_patient <- Hmisc::.q(
    Female_gender, Age_at_Tx, LD_Tx, Donor_Age, MM_ABDR,
    Age_at_Screening, Years_to_ScreeningVisit, Screening_GFR, Screening_Prot_Krea,
    DSA_HLA_class_I_Only_Screening, DSA_HLA_class_II_Only_Screening, DSA_HLA_class_I_and_II_Screening, HLA_DQ_DSA_Screening,
    # mfi_immunodominant,
    HLA_DQ_DSA_Screening, DSA_Number_Screening
)
data %>% dplyr::select(contains("mfi"))


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


# WRANGLE THE AMAT1 SCORES FROM SPSS DATA ####
data_AMAT1 <- data %>%
    dplyr::select(
        contains("AMAT1_1208Set"),
        Trial_Center,
        STUDY_EVALUATION_ID,
        Felzartamab
    ) %>%
    mutate_at(vars(contains("AMAT1_1208Set")), ~ as.numeric(.) %>% suppressWarnings()) %>%
    pivot_longer(
        cols = contains("AMAT1_1208Set"),
        names_to = c("Group"),
        names_pattern = "^(\\w+)_.*",
        values_to = "AMAT1"
    ) %>%
    mutate(
        Group = Group %>%
            str_remove_all("_AMAT1"),
        Followup = case_when(
            Group == "Index" ~ "Day0",
            Group == "FU0" ~ "Week12",
            Group == "FU1" ~ "Week24",
            Group == "FU2" ~ "Week52",
        )
    ) %>%
    dplyr::filter(Followup %>% str_detect("FU1b", negate = TRUE))


# WRANGLE THE QCMAT SCORES FROM SPSS DATA ####
data_QCMAT <- data %>%
    dplyr::select(
        contains("QCMAT_1208Set"),
        Trial_Center,
        STUDY_EVALUATION_ID,
        Felzartamab
    ) %>%
    mutate_at(vars(contains("QCMAT_1208Set")), ~ as.numeric(.) %>% suppressWarnings()) %>%
    pivot_longer(
        cols = contains("QCMAT_1208Set"),
        names_to = c("Group"),
        names_pattern = "^(\\w+)_.*",
        values_to = "QCMAT"
    ) %>%
    mutate(
        Group = Group %>%
            str_remove_all("_QCMAT"),
        Followup = case_when(
            Group == "Index" ~ "Day0",
            Group == "FU0" ~ "Week12",
            Group == "FU1" ~ "Week24",
            Group == "FU2" ~ "Week52",
        )
    ) %>%
    dplyr::filter(Followup %>% str_detect("FU1b", negate = TRUE))


# WRANGLE THE BAT SCORES FROM SPSS DATA ####
data_BAT <- data %>%
    dplyr::select(
        contains("BAT_1208Set"),
        Trial_Center,
        STUDY_EVALUATION_ID,
        Felzartamab
    ) %>%
    mutate_at(vars(contains("BAT_1208Set")), ~ as.numeric(.) %>% suppressWarnings()) %>%
    pivot_longer(
        cols = contains("BAT_1208Set"),
        names_to = c("Group"),
        names_pattern = "^(\\w+)_.*",
        values_to = "BAT"
    ) %>%
    mutate(
        Group = Group %>%
            str_remove_all("_BAT"),
        Followup = case_when(
            Group == "Index" ~ "Day0",
            Group == "FU0" ~ "Week12",
            Group == "FU1" ~ "Week24",
            Group == "FU2" ~ "Week52",
        )
    ) %>%
    dplyr::filter(Followup %>% str_detect("FU1b", negate = TRUE))


# WRANGLE THE cigt1 SCORES FROM SPSS DATA ####
data_cigt1 <- data %>%
    dplyr::select(
        contains("cigt1_1208Set"),
        Trial_Center,
        STUDY_EVALUATION_ID,
        Felzartamab
    ) %>%
    mutate_at(vars(contains("cigt1_1208Set")), ~ as.numeric(.) %>% suppressWarnings()) %>%
    pivot_longer(
        cols = contains("cigt1_1208Set"),
        names_to = c("Group"),
        names_pattern = "^(\\w+)_.*",
        values_to = "cigt1"
    ) %>%
    mutate(
        Group = Group %>%
            str_remove_all("_cigt1"),
        Followup = case_when(
            Group == "Index" ~ "Day0",
            Group == "FU0" ~ "Week12",
            Group == "FU1" ~ "Week24",
            Group == "FU2" ~ "Week52",
        )
    ) %>%
    dplyr::filter(Followup %>% str_detect("FU1b", negate = TRUE))


# WRANGLE THE ctgt1 SCORES FROM SPSS DATA ####
data_ctgt1 <- data %>%
    dplyr::select(
        contains("ctgt1_1208Set"),
        Trial_Center,
        STUDY_EVALUATION_ID,
        Felzartamab
    ) %>%
    mutate_at(vars(contains("ctgt1_1208Set")), ~ as.numeric(.) %>% suppressWarnings()) %>%
    pivot_longer(
        cols = contains("ctgt1_1208Set"),
        names_to = c("Group"),
        names_pattern = "^(\\w+)_.*",
        values_to = "ctgt1"
    ) %>%
    mutate(
        Group = Group %>%
            str_remove_all("_ctgt1"),
        Followup = case_when(
            Group == "Index" ~ "Day0",
            Group == "FU0" ~ "Week12",
            Group == "FU1" ~ "Week24",
            Group == "FU2" ~ "Week52",
        )
    ) %>%
    dplyr::filter(Followup %>% str_detect("FU1b", negate = TRUE))


# WRANGLE THE IGT SCORES FROM SPSS DATA ####
data_IGT <- data %>%
    dplyr::select(
        contains("IGT_1208Set"),
        Trial_Center,
        STUDY_EVALUATION_ID,
        Felzartamab
    ) %>%
    mutate_at(vars(contains("IGT_1208Set")), ~ as.numeric(.) %>% suppressWarnings()) %>%
    pivot_longer(
        cols = contains("IGT_1208Set"),
        names_to = c("Group"),
        names_pattern = "^(\\w+)_.*",
        values_to = "IGT"
    ) %>%
    mutate(
        Group = Group %>%
            str_remove_all("_IGT"),
        Followup = case_when(
            Group == "Index" ~ "Day0",
            Group == "FU0" ~ "Week12",
            Group == "FU1" ~ "Week24",
            Group == "FU2" ~ "Week52",
        )
    ) %>%
    dplyr::filter(Followup %>% str_detect("FU1b", negate = TRUE))


# WRANGLE THE MCAT SCORES FROM SPSS DATA ####
data_MCAT <- data %>%
    dplyr::select(
        contains("MCAT_1208Set"),
        Trial_Center,
        STUDY_EVALUATION_ID,
        Felzartamab
    ) %>%
    mutate_at(vars(contains("MCAT_1208Set")), ~ as.numeric(.) %>% suppressWarnings()) %>%
    pivot_longer(
        cols = contains("MCAT_1208Set"),
        names_to = c("Group"),
        names_pattern = "^(\\w+)_.*",
        values_to = "MCAT"
    ) %>%
    mutate(
        Group = Group %>%
            str_remove_all("_MCAT"),
        Followup = case_when(
            Group == "Index" ~ "Day0",
            Group == "FU0" ~ "Week12",
            Group == "FU1" ~ "Week24",
            Group == "FU2" ~ "Week52",
        )
    ) %>%
    dplyr::filter(Followup %>% str_detect("FU1b", negate = TRUE))



# WRANGLE THE ENDAT SCORES FROM SPSS DATA ####
data_ENDAT <- data %>%
    dplyr::select(
        contains("ENDAT_1208Set"),
        Trial_Center,
        STUDY_EVALUATION_ID,
        Felzartamab
    ) %>%
    mutate_at(vars(contains("ENDAT_1208Set")), ~ as.numeric(.) %>% suppressWarnings()) %>%
    pivot_longer(
        cols = contains("ENDAT_1208Set"),
        names_to = c("Group"),
        names_pattern = "^(\\w+)_.*",
        values_to = "ENDAT"
    ) %>%
    mutate(
        Group = Group %>%
            str_remove_all("_ENDAT"),
        Followup = case_when(
            Group == "Index" ~ "Day0",
            Group == "FU0" ~ "Week12",
            Group == "FU1" ~ "Week24",
            Group == "FU2" ~ "Week52",
        )
    ) %>%
    dplyr::filter(Followup %>% str_detect("FU1b", negate = TRUE))

# WRANGLE THE KT1 SCORES FROM SPSS DATA ####
data_KT1 <- data %>%
    dplyr::select(
        contains("KT1_1208Set"),
        Trial_Center,
        STUDY_EVALUATION_ID,
        Felzartamab
    ) %>%
    mutate_at(vars(contains("KT1_1208Set")), ~ as.numeric(.) %>% suppressWarnings()) %>%
    pivot_longer(
        cols = contains("KT1_1208Set"),
        names_to = c("Group"),
        names_pattern = "^(\\w+)_.*",
        values_to = "KT1"
    ) %>%
    mutate(
        Group = Group %>%
            str_remove_all("_KT1"),
        Followup = case_when(
            Group == "Index" ~ "Day0",
            Group == "FU0" ~ "Week12",
            Group == "FU1" ~ "Week24",
            Group == "FU2" ~ "Week52",
        )
    ) %>%
    dplyr::filter(Followup %>% str_detect("FU1b", negate = TRUE))

# WRANGLE THE KT2 SCORES FROM SPSS DATA ####
data_KT2 <- data %>%
    dplyr::select(
        contains("KT2_1208Set"),
        Trial_Center,
        STUDY_EVALUATION_ID,
        Felzartamab
    ) %>%
    mutate_at(vars(contains("KT2_1208Set")), ~ as.numeric(.) %>% suppressWarnings()) %>%
    pivot_longer(
        cols = contains("KT2_1208Set"),
        names_to = c("Group"),
        names_pattern = "^(\\w+)_.*",
        values_to = "KT2"
    ) %>%
    mutate(
        Group = Group %>%
            str_remove_all("_KT2"),
        Followup = case_when(
            Group == "Index" ~ "Day0",
            Group == "FU0" ~ "Week12",
            Group == "FU1" ~ "Week24",
            Group == "FU2" ~ "Week52",
        )
    ) %>%
    dplyr::filter(Followup %>% str_detect("FU1b", negate = TRUE))





# WRANGLE THE FICOL SCORES FROM SPSS DATA ####
data_FICOL <- data %>%
    dplyr::select(
        contains("FICOL_1208Set"),
        Trial_Center,
        STUDY_EVALUATION_ID,
        Felzartamab
    ) %>%
    mutate_at(vars(contains("FICOL_1208Set")), ~ as.numeric(.) %>% suppressWarnings()) %>%
    pivot_longer(
        cols = contains("FICOL_1208Set"),
        names_to = c("Group"),
        names_pattern = "^(\\w+)_.*",
        values_to = "FICOL"
    ) %>%
    mutate(
        Group = Group %>%
            str_remove_all("_FICOL"),
        Followup = case_when(
            Group == "Index" ~ "Day0",
            Group == "FU0" ~ "Week12",
            Group == "FU1" ~ "Week24",
            Group == "FU2" ~ "Week52",
        )
    ) %>%
    dplyr::filter(Followup %>% str_detect("FU1b", negate = TRUE))


# WRANGLE THE IRRAT30 SCORES FROM SPSS DATA ####
data_IRRAT30 <- data %>%
    dplyr::select(
        contains("IRRAT30_1208Set"),
        Trial_Center,
        STUDY_EVALUATION_ID,
        Felzartamab
    ) %>%
    mutate_at(vars(contains("IRRAT30_1208Set")), ~ as.numeric(.) %>% suppressWarnings()) %>%
    pivot_longer(
        cols = contains("IRRAT30_1208Set"),
        names_to = c("Group"),
        names_pattern = "^(\\w+)_.*",
        values_to = "IRRAT30"
    ) %>%
    mutate(
        Group = Group %>%
            str_remove_all("_IRRAT30"),
        Followup = case_when(
            Group == "Index" ~ "Day0",
            Group == "FU0" ~ "Week12",
            Group == "FU1" ~ "Week24",
            Group == "FU2" ~ "Week52",
        )
    ) %>%
    dplyr::filter(Followup %>% str_detect("FU1b", negate = TRUE))

# WRANGLE THE IRITD3 SCORES FROM SPSS DATA ####
data_IRITD3 <- data %>%
    dplyr::select(
        contains("IRITD3_1208Set"),
        Trial_Center,
        STUDY_EVALUATION_ID,
        Felzartamab
    ) %>%
    mutate_at(vars(contains("IRITD3_1208Set")), ~ as.numeric(.) %>% suppressWarnings()) %>%
    pivot_longer(
        cols = contains("IRITD3_1208Set"),
        names_to = c("Group"),
        names_pattern = "^(\\w+)_.*",
        values_to = "IRITD3"
    ) %>%
    mutate(
        Group = Group %>%
            str_remove_all("_IRITD3"),
        Followup = case_when(
            Group == "Index" ~ "Day0",
            Group == "FU0" ~ "Week12",
            Group == "FU1" ~ "Week24",
            Group == "FU2" ~ "Week52",
        )
    ) %>%
    dplyr::filter(Followup %>% str_detect("FU1b", negate = TRUE))


# WRANGLE THE IRITD5 SCORES FROM SPSS DATA ####
data_IRITD5 <- data %>%
    dplyr::select(
        contains("IRITD5_1208Set"),
        Trial_Center,
        STUDY_EVALUATION_ID,
        Felzartamab
    ) %>%
    mutate_at(vars(contains("IRITD5_1208Set")), ~ as.numeric(.) %>% suppressWarnings()) %>%
    pivot_longer(
        cols = contains("IRITD5_1208Set"),
        names_to = c("Group"),
        names_pattern = "^(\\w+)_.*",
        values_to = "IRITD5"
    ) %>%
    mutate(
        Group = Group %>%
            str_remove_all("_IRITD5"),
        Followup = case_when(
            Group == "Index" ~ "Day0",
            Group == "FU0" ~ "Week12",
            Group == "FU1" ~ "Week24",
            Group == "FU2" ~ "Week52",
        )
    ) %>%
    dplyr::filter(Followup %>% str_detect("FU1b", negate = TRUE))


# WRANGLE THE GRIT3 SCORES FROM SPSS DATA ####
data_GRIT3 <- data %>%
    dplyr::select(
        contains("GRIT3_1208Set"),
        Trial_Center,
        STUDY_EVALUATION_ID,
        Felzartamab
    ) %>%
    mutate_at(vars(contains("GRIT3_1208Set")), ~ as.numeric(.) %>% suppressWarnings()) %>%
    pivot_longer(
        cols = contains("GRIT3_1208Set"),
        names_to = c("Group"),
        names_pattern = "^(\\w+)_.*",
        values_to = "GRIT3"
    ) %>%
    mutate(
        Group = Group %>%
            str_remove_all("_GRIT3"),
        Followup = case_when(
            Group == "Index" ~ "Day0",
            Group == "FU0" ~ "Week12",
            Group == "FU1" ~ "Week24",
            Group == "FU2" ~ "Week52",
        )
    ) %>%
    dplyr::filter(Followup %>% str_detect("FU1b", negate = TRUE))


# WRANGLE THE Rej_RAT SCORES FROM SPSS DATA ####
data_Rej_RAT <- data %>%
    dplyr::select(
        contains("Rej_RAT_1208Set"),
        Trial_Center,
        STUDY_EVALUATION_ID,
        Felzartamab
    ) %>%
    mutate_at(vars(contains("Rej_RAT_1208Set")), ~ as.numeric(.) %>% suppressWarnings()) %>%
    pivot_longer(
        cols = contains("Rej_RAT_1208Set"),
        names_to = c("Group"),
        names_pattern = "^(\\w+)_.*",
        values_to = "Rej_RAT"
    ) %>%
    mutate(
        Group = Group %>%
            str_remove_all("_Rej_RAT"),
        Followup = case_when(
            Group == "Index" ~ "Day0",
            Group == "FU0" ~ "Week12",
            Group == "FU1" ~ "Week24",
            Group == "FU2" ~ "Week52",
        )
    ) %>%
    dplyr::filter(Followup %>% str_detect("FU1b", negate = TRUE))

# WRANGLE THE RejAA_NR SCORES FROM SPSS DATA ####
data_RejAA_NR <- data %>%
    dplyr::select(
        contains("RejAA_NR_1208Set"),
        Trial_Center,
        STUDY_EVALUATION_ID,
        Felzartamab
    ) %>%
    mutate_at(vars(contains("RejAA_NR_1208Set")), ~ as.numeric(.) %>% suppressWarnings()) %>%
    pivot_longer(
        cols = contains("RejAA_NR_1208Set"),
        names_to = c("Group"),
        names_pattern = "^(\\w+)_.*",
        values_to = "RejAA_NR"
    ) %>%
    mutate(
        Group = Group %>%
            str_remove_all("_RejAA_NR"),
        Followup = case_when(
            Group == "Index" ~ "Day0",
            Group == "FU0" ~ "Week12",
            Group == "FU1" ~ "Week24",
            Group == "FU2" ~ "Week52",
        )
    ) %>%
    dplyr::filter(Followup %>% str_detect("FU1b", negate = TRUE))



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
        data_ABMRpm, data_ggt0, data_ptcgt0, data_NKB, data_DSAST,
        data_TCMRt, data_tgt1, data_igt1, data_TCB, data_QCAT,
        data_AMAT1, data_QCMAT,
        data_BAT, data_cigt1, data_ctgt1, data_IGT, data_MCAT,
        data_ENDAT, data_KT1, data_KT2, 
        data_FICOL, data_IRRAT30, data_IRITD3, data_IRITD5, data_GRIT3, 
        data_Rej_RAT, data_RejAA_NR,
        data_ABMRActivity_Banff19, data_Active_ABMR_Banff19, data_Chronic_active_ABMR_Banff19,
        data_Borderline_Banff
    ),
    left_join,
    by = c("Trial_Center", "STUDY_EVALUATION_ID", "Felzartamab", "Group", "Followup")
) %>%
    left_join(data_patient, by = c("Trial_Center", "STUDY_EVALUATION_ID", "Felzartamab")) %>%
    mutate(
        Trial_Center = Trial_Center %>% as.character(),
        STUDY_EVALUATION_ID = STUDY_EVALUATION_ID %>% as.numeric(),
        Felzartamab = Felzartamab %>% as.numeric(),
        Group = Group  %>% factor(levels = c("Index", "FU1", "FU2"))
    ) %>%
    relocate(names(data_patient), .before = 1)



# SAVE THE DATA ####
save(data_K1208, file = "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0Georg Felz CD38 Vienna/G_Rstuff/data/data_K1208.RData")



# # JOIN THE SPSS REFERENCE SET DATA ####
# data_joined <- data_spss %>%
#     left_join(data_reference,
#         by = c("CEL", "Trial_Center", "STUDY_EVALUATION_ID", "Felzartamab", "Group")
#     )


# data_joined %>%
#     print(n="all")
