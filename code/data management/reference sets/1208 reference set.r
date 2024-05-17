# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(haven) # install.packages("haven")
# Bioconductor libraries
library(Biobase) # BiocManager::install("Biobase")
# Custom operators, functions, and datasets
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
source("C:/R/CD38-effect-of-treatment/code/functions/complex_pivot.R")
# Load new reference set
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0Georg Felz CD38 Vienna/G_Rstuff/data/OLD_Felz_6March24.RData")
# load SPSS data from NEJM paper
data <- read_spss("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0Georg Felz CD38 Vienna/G_Rstuff/data/Generalfile_Felzartamab SPSS.sav")
# load reference data
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/data_scores_k1208.RData")
data_scores_k1208 %>% colnames()



# DEFINE VARIABLES OF INTEREST ####
vars <- Hmisc::.q(
    STUDY_EVALUATION_ID,
    Felzartamab_presumed,
    Center,
    Date_birth,
    Date_Tx,
    Date_Bx
)

# DEFINE SET ####
set <- NEW
set$CEL <- set %>% sampleNames()
# NEW %>% pData  %>% colnames


# WRANGLE THE CEL IDS FROM SPSS DATA ####
cels <- data %>%
    dplyr::select(
        Index_CEL, FU1_CEL, FU1bBx_CEL, FU2_CEL
    ) %>%
    pivot_longer(
        cols = contains("CEL"),
        names_to = "Group",
        names_pattern = "^(\\w+)_.*",
        values_to = "CEL"
    ) %>%
    mutate(
        Group = Group %>% str_remove("Bx")
    ) %>%
    mutate_at("CEL", na_if, "") %>%
    drop_na(CEL) %>%
    pull(CEL)

set$CEL[set$CEL %nin% cels]


# TRIM SET TO PHENOTYPE DATA ####
set_felzartamab <- set[, which(set$CEL %in% c(cels))]
# set_felzartamab  %>% pData



# WRANGLE THE PHENOTYPE DATA ####
set_felzartamab$STUDY_EVALUATION_ID <- NA
set_felzartamab$Felzartamab_presumed <- NA
set_felzartamab$Center <- NA
set_felzartamab$Date_birth <- as.Date(NA)
set_felzartamab$Date_Tx <- as.Date(NA)
set_felzartamab$Date_Bx <- as.Date(NA)
set_felzartamab$cg <- NA
set_felzartamab$MVI_Score <- NA
set_felzartamab$ID_MMDx <- NA
set_felzartamab$Group <- NA

# set_felzartamab$DateBx <- NULL # REDCap entries
# set_felzartamab$DateTx <- NULL

sel <- which(IDs$IndexCEL %in% sampleNames(set_felzartamab))
sel2 <- match(IDs$IndexCEL[sel], sampleNames(set_felzartamab))
set_felzartamab$STUDY_EVALUATION_ID[sel2] <- IDs$STUDY_EVALUATION_ID[sel]
set_felzartamab$Felzartamab_presumed[sel2] <- IDs$Felzartamab_presumed[sel]
set_felzartamab$Center[sel2] <- IDs$Center[sel]
set_felzartamab$Date_birth[sel2] <- IDs$Date_birth[sel]
set_felzartamab$Date_Tx[sel2] <- IDs$Date_Tx[sel]
set_felzartamab$Date_Bx[sel2] <- IDs$IndexBx_Date[sel]
set_felzartamab$cg[sel2] <- IDs$IndexBx_cg[sel]
set_felzartamab$MVI_Score[sel2] <- IDs$IndexBx_MVI_Score[sel]
set_felzartamab$ID_MMDx[sel2] <- IDs$IndexBx_ID_MMDx[sel]
set_felzartamab$Group[sel2] <- "Index"

sel <- which(IDs$FU1CEL %in% sampleNames(set_felzartamab))
sel2 <- match(IDs$FU1CEL[sel], sampleNames(set_felzartamab))
set_felzartamab$STUDY_EVALUATION_ID[sel2] <- IDs$STUDY_EVALUATION_ID[sel]
set_felzartamab$Felzartamab_presumed[sel2] <- IDs$Felzartamab_presumed[sel]
set_felzartamab$Center[sel2] <- IDs$Center[sel]
set_felzartamab$Date_birth[sel2] <- IDs$Date_birth[sel]
set_felzartamab$Date_Tx[sel2] <- IDs$Date_Tx[sel]
set_felzartamab$Date_Bx[sel2] <- IDs$FU1Bx_Date[sel]
set_felzartamab$cg[sel2] <- IDs$FU1Bx_cg[sel]
set_felzartamab$MVI_Score[sel2] <- IDs$FU1Bx_MVI_Score[sel]
set_felzartamab$ID_MMDx[sel2] <- IDs$FU1Bx_ID_MMDx[sel]
set_felzartamab$Group[sel2] <- "FU1"

sel <- which(IDs$FU2CEL %in% sampleNames(set_felzartamab))
sel2 <- match(IDs$FU2CEL[sel], sampleNames(set_felzartamab))
set_felzartamab$STUDY_EVALUATION_ID[sel2] <- IDs$STUDY_EVALUATION_ID[sel]
set_felzartamab$Felzartamab_presumed[sel2] <- IDs$Felzartamab_presumed[sel]
set_felzartamab$Center[sel2] <- IDs$Center[sel]
set_felzartamab$Date_birth[sel2] <- IDs$Date_birth[sel]
set_felzartamab$Date_Tx[sel2] <- IDs$Date_Tx[sel]
set_felzartamab$Date_Bx[sel2] <- IDs$FU2Bx_Date[sel]
# set_felzartamab$cg[sel2]<-IDs$FU2Bx_cg[sel]
# set_felzartamab$MVI_Score[sel2]<-IDs$FU2Bx_MVI_Score[sel]
set_felzartamab$ID_MMDx[sel2] <- IDs$FU2Bx_ID_MMDx[sel]
set_felzartamab$Group[sel2] <- "FU2"


set_felzartamab$TxBx <- set_felzartamab$Date_Bx - set_felzartamab$Date_Tx
set_felzartamab$RejAA_NR <- set_felzartamab$AARej6[, 1]
set_felzartamab$RejAA_TCMR <- set_felzartamab$AARej6[, 2]
set_felzartamab$RejAA_Mixed <- set_felzartamab$AARej6[, 3]
set_felzartamab$RejAA_EABMR <- set_felzartamab$AARej6[, 4]
set_felzartamab$RejAA_FABMR <- set_felzartamab$AARej6[, 5]
set_felzartamab$RejAA_LABMR <- set_felzartamab$AARej6[, 6]
set_felzartamab$RejAA_Clust <- ifelse(set_felzartamab$AARejClust == 1, "NR",
    ifelse(set_felzartamab$AARejClust == 2, "TCMR",
        ifelse(set_felzartamab$AARejClust == 3, "Mixed",
            ifelse(set_felzartamab$AARejClust == 4, "EABMR",
                ifelse(set_felzartamab$AARejClust == 5, "FABMR",
                    ifelse(set_felzartamab$AARejClust == 6, "LABMR", "Missing")
                )
            )
        )
    )
)

# set_felzartamab$InjAA_MildCKD <- set_felzartamab$Inj5Score[, 1]
# set_felzartamab$InjAA_CKDAKI <- set_felzartamab$Inj5Score[, 2]
# set_felzartamab$InjAA_AKI1 <- set_felzartamab$Inj5Score[, 3]
# set_felzartamab$InjAA_AKI2 <- set_felzartamab$Inj5Score[, 4]
# set_felzartamab$InjAA_Normal <- set_felzartamab$Inj5Score[, 5]

# set_felzartamab$InjAA_Clust <- ifelse(set_felzartamab$Inj5Clust == 1, "MildCKD",
#     ifelse(set_felzartamab$Inj5Clust == 2, "CKDAKI",
#         ifelse(set_felzartamab$Inj5Clust == 3, "AKI1",
#             ifelse(set_felzartamab$Inj5Clust == 4, "AKI2",
#                 ifelse(set_felzartamab$Inj5Clust == 5, "Normal", "Missing")
#             )
#         )
#     )
# )


# SAVE NEW DATALOCK ####
vienna_1208 <- set_felzartamab
save(vienna_1208, file = "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0Georg Felz CD38 Vienna/G_Rstuff/data/vienna_1208_6Mar24.RData")


# CREATE AN EXCEL FILE WITH DATA ####
df_vienna_1208 <- tibble(
    CEL = sampleNames(vienna_1208), ID = vienna_1208$STUDY_EVALUATION_ID, Group = vienna_1208$Group, Center = vienna_1208$Center,
    ID_MMDx = vienna_1208$ID_MMDx,
    Felz_presumed = vienna_1208$Felzartamab_presumed, Date_birth = vienna_1208$Date_birth, Date_Tx = vienna_1208$Date_Tx,
    Date_Bx = vienna_1208$Date_Bx, TxBx = vienna_1208$TxBx, cg = vienna_1208$cg, MVI_Score = vienna_1208$MVI_Score,
    RejAA_Clust = vienna_1208$RejAA_Clust,
    RejAA_NR = vienna_1208$RejAA_NR, RejAA_TCMR = vienna_1208$RejAA_TCMR,
    RejAA_Mixed = vienna_1208$RejAA_Mixed, RejAA_EABMR = vienna_1208$RejAA_EABMR,
    RejAA_FABMR = vienna_1208$RejAA_FABMR, RejAA_LABMR = vienna_1208$RejAA_LABMR, # RejAA_MinorABMR = vienna_1208$RejAA_MinorABMR,
    # InjAA_Clust = vienna_1208$InjAA_Clust, InjAA_MildCKD = vienna_1208$InjAA_MildCKD, InjAA_CKDAKI = vienna_1208$InjAA_CKDAKI,
    # InjAA_AKI1 = vienna_1208$InjAA_AKI1, InjAA_AKI2 = vienna_1208$InjAA_AKI2, InjAA_Normal = vienna_1208$InjAA_Normal,
    ABMRpm = vienna_1208$ABMRpm, TCMRt = vienna_1208$TCMRt, Rejr = vienna_1208$Rejr, tgt1 = vienna_1208$tgt1, igt1 = vienna_1208$igt1,
    ggt0 = vienna_1208$ggt0, cggt0 = vienna_1208$cggt0, ptcgt0 = vienna_1208$ptcgt0, cigt1 = vienna_1208$cigt1, ctgt1 = vienna_1208$ctgt1,
    mmgt1 = vienna_1208$mmgt1, ahgt0 = vienna_1208$ahgt0, Protprob = vienna_1208$Protprob, lowGFRprob = vienna_1208$Protprob,
    DSAprob = vienna_1208$DSAprob, Surv3yrprob = vienna_1208$S3Combmin, RFTCMRprob = vienna_1208$RFTCMR, RFABMRprob = vienna_1208$RFABMR,
    InjPC1 = vienna_1208$InjPC1, InjPC2 = vienna_1208$InjPC2, InjPC3 = vienna_1208$InjPC3, RejPC1 = vienna_1208$RejPC1, RejPC2 = vienna_1208$RejPC2,
    RejPC3 = vienna_1208$RejPC3, GlobalDist = vienna_1208$GD, Cortexprob = vienna_1208$Cort,
    ABMR_RAT = vienna_1208$`ABMR-RAT`, AMAT1 = vienna_1208$AMAT1, BAT = vienna_1208$BAT, DAMP = vienna_1208$DAMP, DSAST = vienna_1208$DSAST, eDSAST = vienna_1208$eDSAST, ENDAT = vienna_1208$ENDAT,
    GRIT1 = vienna_1208$GRIT1, GRIT2 = vienna_1208$GRIT2, GRIT3 = vienna_1208$GRIT3, IGT = vienna_1208$IGT, IRITD3 = vienna_1208$IRITD3, IRITD5 = vienna_1208$IRITD5,
    IRRAT30 = vienna_1208$IRRAT30, KT1 = vienna_1208$KT1, KT2 = vienna_1208$KT2,
    MCAT = vienna_1208$MCAT, NKB = vienna_1208$NKB, QCAT = vienna_1208$QCAT, QCMAT = vienna_1208$QCMAT,
    Rej_RAT = vienna_1208$`Rej-RAT`, TCB = vienna_1208$TCB, TCMR_RAT = vienna_1208$`TCMR-RAT`, FICOL = vienna_1208$FICOL
) %>%
    mutate_if(is.numeric, ~ round(., 3))


# SAVE THE EXCEL FILE ####
write.xlsx(df_vienna_1208, file = "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0Georg Felz CD38 Vienna/G_Rstuff/output/dfvienna_1208_6Mar24.xlsx")
