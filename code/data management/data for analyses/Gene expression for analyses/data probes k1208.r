# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(openxlsx) # install.packages("openxlsx")
# Bioconductor libraries
library(Biobase) # BiocManager::install("Biobase")
# Custom operators, functions, and datasets
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
# load the reference set
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/data_cfdna.RData")




# DEFINE SET ####
set <- NEW
# NEW %>% pData  %>% colnames


# DEFINE INDEX AND FOLLOW-UP CELS ####
cels_index <- IDs %>%
    dplyr::filter(IndexCEL  %>% str_detect(".CEL"))  %>% 
    pull(IndexCEL) 

cels_fu1 <- IDs %>%
    dplyr::filter(FU1CEL %>% str_detect(".CEL")) %>%
    pull(FU1CEL)

cels_fu2 <- IDs %>%
    dplyr::filter(FU2CEL %>% str_detect(".CEL")) %>%
    pull(FU2CEL)



# TRIM SET TO PHENOTYPE DATA ####
set$CEL <- set  %>% sampleNames
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
set_vienna$RejAA_NR <- set_vienna$AARej6[, 1]
set_vienna$RejAA_TCMR <- set_vienna$AARej6[, 2]
set_vienna$RejAA_Mixed <- set_vienna$AARej6[, 3]
set_vienna$RejAA_EABMR <- set_vienna$AARej6[, 4]
set_vienna$RejAA_FABMR <- set_vienna$AARej6[, 5]
set_vienna$RejAA_LABMR <- set_vienna$AARej6[, 6]
set_vienna$RejAA_Clust <- ifelse(set_vienna$AARejClust == 1, "NR",
    ifelse(set_vienna$AARejClust == 2, "TCMR",
        ifelse(set_vienna$AARejClust == 3, "Mixed",
            ifelse(set_vienna$AARejClust == 4, "EABMR",
                ifelse(set_vienna$AARejClust == 5, "FABMR",
                    ifelse(set_vienna$AARejClust == 6, "LABMR", "Missing")
                    
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
vienna_1208 <- set_vienna
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
    RejAA_FABMR = vienna_1208$RejAA_FABMR, RejAA_LABMR = vienna_1208$RejAA_LABMR, #RejAA_MinorABMR = vienna_1208$RejAA_MinorABMR,
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
)%>%
    mutate_if(is.numeric, ~ round(., 3))


# SAVE THE EXCEL FILE ####
write.xlsx(df_vienna_1208, file = "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0Georg Felz CD38 Vienna/G_Rstuff/output/dfvienna_1208_6Mar24.xlsx")
