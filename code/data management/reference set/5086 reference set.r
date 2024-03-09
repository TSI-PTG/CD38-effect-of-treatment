# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(openxlsx)
# Bioconductor libraries
library(Biobase) # BiocManager::install("Biobase")
# Custom operators, functions, and datasets
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
# Load new reference set
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0Georg Felz CD38 Vienna/G_Rstuff/data/NEW_Felz_6March24.RData") # N=52
# Load previous reference set
# load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0Georg Felz CD38 Vienna/G_Rstuff/data/Vienna.RData") # N=52
# Load additional phenotype data
IDs <- read.xlsx("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0Georg Felz CD38 Vienna/G_Rstuff/data/Felzartamab excel sheet Updated 5March2024.xlsx", detectDates = T)


# DEFINE SET ####
set <- NEW
# NEW %>% pData  %>% colnames


# DEFINE INDEX AND FOLLOW-UP CELS ####
cels_index <- IDs %>%
    dplyr::filter(IndexCEL %>% str_detect(".CEL")) %>%
    pull(IndexCEL)

cels_fu1 <- IDs %>%
    dplyr::filter(FU1CEL %>% str_detect(".CEL")) %>%
    pull(FU1CEL)

cels_fu2 <- IDs %>%
    dplyr::filter(FU2CEL %>% str_detect(".CEL")) %>%
    pull(FU2CEL)


# TRIM SET TO PHENOTYPE DATA ####
set$CEL <- set %>% sampleNames()
set_vienna <- set[, which(set$CEL %in% c(cels_index, cels_fu1, cels_fu2))]
# set_vienna  %>% pData


# WRANGLE THE PHENOTYPE DATA ####
set_vienna$STUDY_EVALUATION_ID <- NA
set_vienna$Felzartamab_presumed <- NA
set_vienna$Center <- NA
set_vienna$Date_birth <- as.Date(NA)
set_vienna$Date_Tx <- as.Date(NA)
set_vienna$Date_Bx <- as.Date(NA)
set_vienna$cg <- NA
set_vienna$MVI_Score <- NA
set_vienna$ID_MMDx <- NA
set_vienna$Group <- NA

set_vienna$DateBx <- NULL # REDCap entries
set_vienna$DateTx <- NULL

sel <- which(IDs$IndexCEL %in% sampleNames(set_vienna))
sel2 <- match(IDs$IndexCEL[sel], sampleNames(set_vienna))
set_vienna$STUDY_EVALUATION_ID[sel2] <- IDs$STUDY_EVALUATION_ID[sel]
set_vienna$Felzartamab_presumed[sel2] <- IDs$Felzartamab_presumed[sel]
set_vienna$Center[sel2] <- IDs$Center[sel]
set_vienna$Date_birth[sel2] <- IDs$Date_birth[sel]
set_vienna$Date_Tx[sel2] <- IDs$Date_Tx[sel]
set_vienna$Date_Bx[sel2] <- IDs$IndexBx_Date[sel]
set_vienna$cg[sel2] <- IDs$IndexBx_cg[sel]
set_vienna$MVI_Score[sel2] <- IDs$IndexBx_MVI_Score[sel]
set_vienna$ID_MMDx[sel2] <- IDs$IndexBx_ID_MMDx[sel]
set_vienna$Group[sel2] <- "Index"

sel <- which(IDs$FU1CEL %in% sampleNames(set_vienna))
sel2 <- match(IDs$FU1CEL[sel], sampleNames(set_vienna))
set_vienna$STUDY_EVALUATION_ID[sel2] <- IDs$STUDY_EVALUATION_ID[sel]
set_vienna$Felzartamab_presumed[sel2] <- IDs$Felzartamab_presumed[sel]
set_vienna$Center[sel2] <- IDs$Center[sel]
set_vienna$Date_birth[sel2] <- IDs$Date_birth[sel]
set_vienna$Date_Tx[sel2] <- IDs$Date_Tx[sel]
set_vienna$Date_Bx[sel2] <- IDs$FU1Bx_Date[sel]
set_vienna$cg[sel2] <- IDs$FU1Bx_cg[sel]
set_vienna$MVI_Score[sel2] <- IDs$FU1Bx_MVI_Score[sel]
set_vienna$ID_MMDx[sel2] <- IDs$FU1Bx_ID_MMDx[sel]
set_vienna$Group[sel2] <- "FU1"

sel <- which(IDs$FU2CEL %in% sampleNames(set_vienna))
sel2 <- match(IDs$FU2CEL[sel], sampleNames(set_vienna))
set_vienna$STUDY_EVALUATION_ID[sel2] <- IDs$STUDY_EVALUATION_ID[sel]
set_vienna$Felzartamab_presumed[sel2] <- IDs$Felzartamab_presumed[sel]
set_vienna$Center[sel2] <- IDs$Center[sel]
set_vienna$Date_birth[sel2] <- IDs$Date_birth[sel]
set_vienna$Date_Tx[sel2] <- IDs$Date_Tx[sel]
set_vienna$Date_Bx[sel2] <- IDs$FU2Bx_Date[sel]
# set_vienna$cg[sel2]<-IDs$FU2Bx_cg[sel]
# set_vienna$MVI_Score[sel2]<-IDs$FU2Bx_MVI_Score[sel]
set_vienna$ID_MMDx[sel2] <- IDs$FU2Bx_ID_MMDx[sel]
set_vienna$Group[sel2] <- "FU2"


set_vienna$TxBx <- set_vienna$Date_Bx - set_vienna$Date_Tx
set_vienna$RejAA_NR <- set_vienna$Rej7Score[, 1]
set_vienna$RejAA_TCMR1 <- set_vienna$Rej7Score[, 2]
set_vienna$RejAA_TCMR2 <- set_vienna$Rej7Score[, 3]
set_vienna$RejAA_EABMR <- set_vienna$Rej7Score[, 4]
set_vienna$RejAA_FABMR <- set_vienna$Rej7Score[, 5]
set_vienna$RejAA_LABMR <- set_vienna$Rej7Score[, 6]
set_vienna$RejAA_MinorABMR <- set_vienna$Rej7Score[, 7]
set_vienna$RejAA_Clust <- ifelse(set_vienna$Rej7Clust == 1, "NR",
    ifelse(set_vienna$Rej7Clust == 2, "TCMR1",
        ifelse(set_vienna$Rej7Clust == 3, "TCMR2",
            ifelse(set_vienna$Rej7Clust == 4, "EABMR",
                ifelse(set_vienna$Rej7Clust == 5, "FABMR",
                    ifelse(set_vienna$Rej7Clust == 6, "LABMR",
                        ifelse(set_vienna$Rej7Clust == 7, "MinorABMR", "Missing")
                    )
                )
            )
        )
    )
)

# set_vienna$InjAA_MildCKD <- set_vienna$Inj5Score[, 1]
# set_vienna$InjAA_CKDAKI <- set_vienna$Inj5Score[, 2]
# set_vienna$InjAA_AKI1 <- set_vienna$Inj5Score[, 3]
# set_vienna$InjAA_AKI2 <- set_vienna$Inj5Score[, 4]
# set_vienna$InjAA_Normal <- set_vienna$Inj5Score[, 5]

# set_vienna$InjAA_Clust <- ifelse(set_vienna$Inj5Clust == 1, "MildCKD",
#     ifelse(set_vienna$Inj5Clust == 2, "CKDAKI",
#         ifelse(set_vienna$Inj5Clust == 3, "AKI1",
#             ifelse(set_vienna$Inj5Clust == 4, "AKI2",
#                 ifelse(set_vienna$Inj5Clust == 5, "Normal", "Missing")
#             )
#         )
#     )
# )


# SAVE NEW DATALOCK ####
vienna_5086 <- set_vienna
save(vienna_5086, file = "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0Georg Felz CD38 Vienna/G_Rstuff/data/vienna_5086_6Mar24.RData")


# CREATE AN EXCEL FILE WITH DATA ####
df_vienna_5086 <- tibble(
    CEL = sampleNames(vienna_5086), ID = vienna_5086$STUDY_EVALUATION_ID, Group = vienna_5086$Group, Center = vienna_5086$Center,
    ID_MMDx = vienna_5086$ID_MMDx,
    Felz_presumed = vienna_5086$Felzartamab_presumed, Date_birth = vienna_5086$Date_birth, Date_Tx = vienna_5086$Date_Tx,
    Date_Bx = vienna_5086$Date_Bx, TxBx = vienna_5086$TxBx, cg = vienna_5086$cg, MVI_Score = vienna_5086$MVI_Score,
    RejAA_Clust = vienna_5086$RejAA_Clust,
    RejAA_NR = vienna_5086$RejAA_NR, 
    RejAA_TCMR1 = vienna_5086$RejAA_TCMR1, RejAA_TCMR2 = vienna_5086$RejAA_TCMR2, 
    RejAA_EABMR = vienna_5086$RejAA_EABMR, RejAA_FABMR = vienna_5086$RejAA_FABMR, RejAA_LABMR = vienna_5086$RejAA_LABMR, 
    RejAA_MinorABMR = vienna_5086$RejAA_MinorABMR,
    # InjAA_Clust = vienna_5086$InjAA_Clust, InjAA_MildCKD = vienna_5086$InjAA_MildCKD, InjAA_CKDAKI = vienna_5086$InjAA_CKDAKI,
    # InjAA_AKI1 = vienna_5086$InjAA_AKI1, InjAA_AKI2 = vienna_5086$InjAA_AKI2, InjAA_Normal = vienna_5086$InjAA_Normal,
    ABMRpm = vienna_5086$ABMRpm, TCMRt = vienna_5086$TCMRt, Rejr = vienna_5086$Rejr, tgt1 = vienna_5086$tgt1, igt1 = vienna_5086$igt1,
    ggt0 = vienna_5086$ggt0, cggt0 = vienna_5086$cggt0, ptcgt0 = vienna_5086$ptcgt0, cigt1 = vienna_5086$cigt1, ctgt1 = vienna_5086$ctgt1,
    mmgt1 = vienna_5086$mmgt1, ahgt0 = vienna_5086$ahgt0, Protprob = vienna_5086$Protprob, lowGFRprob = vienna_5086$Protprob,
    DSAprob = vienna_5086$DSAprob, Surv3yrprob = vienna_5086$S3Combmin, RFTCMRprob = vienna_5086$RFTCMR, RFABMRprob = vienna_5086$RFABMR,
    InjPC1 = vienna_5086$InjPC1, InjPC2 = vienna_5086$InjPC2, InjPC3 = vienna_5086$InjPC3, RejPC1 = vienna_5086$RejPC1, RejPC2 = vienna_5086$RejPC2,
    RejPC3 = vienna_5086$RejPC3, GlobalDist = vienna_5086$GD, Cortexprob = vienna_5086$Cort,
    ABMR_RAT = vienna_5086$`ABMR-RAT`, AMAT1 = vienna_5086$AMAT1, BAT = vienna_5086$BAT, DAMP = vienna_5086$DAMP, DSAST = vienna_5086$DSAST, eDSAST = vienna_5086$eDSAST, ENDAT = vienna_5086$ENDAT,
    GRIT1 = vienna_5086$GRIT1, GRIT2 = vienna_5086$GRIT2, GRIT3 = vienna_5086$GRIT3, IGT = vienna_5086$IGT, IRITD3 = vienna_5086$IRITD3, IRITD5 = vienna_5086$IRITD5,
    IRRAT30 = vienna_5086$IRRAT30, KT1 = vienna_5086$KT1, KT2 = vienna_5086$KT2,
    MCAT = vienna_5086$MCAT, NKB = vienna_5086$NKB, QCAT = vienna_5086$QCAT, QCMAT = vienna_5086$QCMAT,
    Rej_RAT = vienna_5086$`Rej-RAT`, TCB = vienna_5086$TCB, TCMR_RAT = vienna_5086$`TCMR-RAT`, FICOL = vienna_5086$FICOL
) %>%
    mutate_if(is.numeric, ~ round(., 3))


# SAVE THE EXCEL FILE ####
write.xlsx(df_vienna_5086, file = "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0Georg Felz CD38 Vienna/G_Rstuff/output/dfvienna_5086_6Mar24.xlsx")
