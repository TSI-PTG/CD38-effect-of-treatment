# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(readxl) # install.packages("readxl")

# Bioconductor libraries
library(Biobase) # BiocManager::install("Biobase")
# Custom operators, functions, and datasets
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
# Load new reference set
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0Georg Felz CD38 Vienna/G_Rstuff/data/Vienna_new_raw_18Oct23.RData") # N=52
# Load additional phenotype data
IDs <- read_excel("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0Georg Felz CD38 Vienna/G_Rstuff/data/Felzartamab_MMDx_File_4Feb 24.xlsx")



dir_cel <- "Z:/CELs/219/KIDNEY INTERCOMEX"
CEL_kidney <- dir_cel  %>% dir

CEL_Vienna_new <- Vienna  %>% sampleNames()


IDs$CEL %in% CEL_Vienna_new



# DEFINE SET ####
set <- NEW


# DEFINE INDEX AND FOLLOW-UP CELS ####
cels_index <- IDs %>%
    dplyr::filter(IndexCEL  %>% str_detect(".CEL"))  %>% 
    pull(IndexCEL) 

cels_fu1 <- IDs %>%
    dplyr::filter(FU1CEL %>% str_detect(".CEL")) %>%
    pull(FU1CEL)


# TRIM SET TO PHENOTYPE DATA ####
set$CEL <- set  %>% sampleNames
set_vienna <- set[, which(set$CEL %in% c(cels_index, cels_fu1))]


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

set_vienna$InjAA_MildCKD <- set_vienna$Inj5Score[, 1]
set_vienna$InjAA_CKDAKI <- set_vienna$Inj5Score[, 2]
set_vienna$InjAA_AKI1 <- set_vienna$Inj5Score[, 3]
set_vienna$InjAA_AKI2 <- set_vienna$Inj5Score[, 4]
set_vienna$InjAA_Normal <- set_vienna$Inj5Score[, 5]

set_vienna$InjAA_Clust <- ifelse(set_vienna$Inj5Clust == 1, "MildCKD",
    ifelse(set_vienna$Inj5Clust == 2, "CKDAKI",
        ifelse(set_vienna$Inj5Clust == 3, "AKI1",
            ifelse(set_vienna$Inj5Clust == 4, "AKI2",
                ifelse(set_vienna$Inj5Clust == 5, "Normal", "Missing")
            )
        )
    )
)


# SAVE NEW DATALOCK ####
Vienna44 <- set_vienna
save(Vienna44, file = "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0Georg Felz CD38 Vienna/G_Rstuff/data/Vienna44_18Oct23.RData")


# CREATE AN EXCEL FILE WITH DATA ####
df_Vienna44 <- tibble(
    CEL = sampleNames(Vienna44), ID = Vienna44$STUDY_EVALUATION_ID, Group = Vienna44$Group, Center = Vienna44$Center,
    ID_MMDx = Vienna44$ID_MMDx,
    Felz_presumed = Vienna44$Felzartamab_presumed, Date_birth = Vienna44$Date_birth, Date_Tx = Vienna44$Date_Tx,
    Date_Bx = Vienna44$Date_Bx, TxBx = Vienna44$TxBx, cg = Vienna44$cg, MVI_Score = Vienna44$MVI_Score,
    RejAA_Clust = Vienna44$RejAA_Clust, RejAA_NR = Vienna44$RejAA_NR, RejAA_TCMR1 = Vienna44$RejAA_TCMR1,
    RejAA_TCMR2 = Vienna44$RejAA_TCMR2, RejAA_EABMR = Vienna44$RejAA_EABMR,
    RejAA_FABMR = Vienna44$RejAA_FABMR, RejAA_LABMR = Vienna44$RejAA_LABMR, RejAA_MinorABMR = Vienna44$RejAA_MinorABMR,
    InjAA_Clust = Vienna44$InjAA_Clust, InjAA_MildCKD = Vienna44$InjAA_MildCKD, InjAA_CKDAKI = Vienna44$InjAA_CKDAKI,
    InjAA_AKI1 = Vienna44$InjAA_AKI1, InjAA_AKI2 = Vienna44$InjAA_AKI2, InjAA_Normal = Vienna44$InjAA_Normal,
    ABMRpm = Vienna44$ABMRpm, TCMRt = Vienna44$TCMRt, Rejr = Vienna44$Rejr, tgt1 = Vienna44$tgt1, igt1 = Vienna44$igt1,
    ggt0 = Vienna44$ggt0, cggt0 = Vienna44$cggt0, ptcgt0 = Vienna44$ptcgt0, cigt1 = Vienna44$cigt1, ctgt1 = Vienna44$ctgt1,
    mmgt1 = Vienna44$mmgt1, ahgt0 = Vienna44$ahgt0, Protprob = Vienna44$Protprob, lowGFRprob = Vienna44$Protprob,
    DSAprob = Vienna44$DSAprob, Surv3yrprob = Vienna44$S3Combmin, RFTCMRprob = Vienna44$RFTCMR, RFABMRprob = Vienna44$RFABMR,
    InjPC1 = Vienna44$InjPC1, InjPC2 = Vienna44$InjPC2, InjPC3 = Vienna44$InjPC3, RejPC1 = Vienna44$RejPC1, RejPC2 = Vienna44$RejPC2,
    RejPC3 = Vienna44$RejPC3, GlobalDist = Vienna44$GD, Cortexprob = Vienna44$Cort,
    ABMR_RAT = Vienna44$`ABMR-RAT`, AMAT1 = Vienna44$AMAT1, BAT = Vienna44$BAT, DAMP = Vienna44$DAMP, DSAST = Vienna44$DSAST, eDSAST = Vienna44$eDSAST, ENDAT = Vienna44$ENDAT,
    GRIT1 = Vienna44$GRIT1, GRIT2 = Vienna44$GRIT2, GRIT3 = Vienna44$GRIT3, IGT = Vienna44$IGT, IRITD3 = Vienna44$IRITD3, IRITD5 = Vienna44$IRITD5,
    IRRAT30 = Vienna44$IRRAT30, KT1 = Vienna44$KT1, KT2 = Vienna44$KT2,
    MCAT = Vienna44$MCAT, NKB = Vienna44$NKB, QCAT = Vienna44$QCAT, QCMAT = Vienna44$QCMAT,
    Rej_RAT = Vienna44$`Rej-RAT`, TCB = Vienna44$TCB, TCMR_RAT = Vienna44$`TCMR-RAT`, FICOL = Vienna44$FICOL
)%>%
    mutate_if(is.numeric, ~ round(., 3))


# SAVE THE EXCEL FILE ####
write.xlsx(df_Vienna44, file = "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0Georg Felz CD38 Vienna/G_Rstuff/output/dfVienna44_18Oct23.xlsx")
