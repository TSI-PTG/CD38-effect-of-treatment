# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
# Custom operators, functions, and datasets
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
# load affymap
load("Z:/DATA/Datalocks/Other data/affymap219_21Oct2019_1306_JR.RData")
# load new PBT lists
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/PBTlist219_ABMRactivity_26Jul2024_PTG.RData")
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/PBTlist219_AKIinduced_26Jul2024_PTG.RData")


# ADD INDIVIDUAL PBTS FROM PBTLISTs TO ENVIRONMENT ####
PBTlist219_AKIinduced %>% list2env(envir = .GlobalEnv)
PBTlist219_ABMRactivity %>% list2env(envir = .GlobalEnv)


# ADD PBT ANNOTATION TO affymap219 ####
PBTlist219_AKIinduced %>% names()
PBTlist219_ABMRactivity %>% names()

affymap219_new <- affymap219 %>%
    mutate(
        PBT = ifelse(AffyID %in% AAG, paste(PBT, "AAG", sep = ","), PBT),
        PBT = ifelse(AffyID %in% IIAAG, paste(PBT, "IIAAG", sep = ","), PBT),
        PBT = ifelse(AffyID %in% NKAAG, paste(PBT, "NKAAG", sep = ","), PBT),
        PBT = ifelse(AffyID %in% AEG, paste(PBT, "AEG", sep = ","), PBT),
        PBT = ifelse(AffyID %in% PT_New1, paste(PBT, "PT_New1", sep = ","), PBT),
        PBT = ifelse(AffyID %in% PT_New2, paste(PBT, "PT_New2", sep = ","), PBT),
        PBT = ifelse(AffyID %in% PT_New3, paste(PBT, "PT_New3", sep = ","), PBT),
        PBT = ifelse(AffyID %in% PT_New4, paste(PBT, "PT_New4", sep = ","), PBT),
        PBT = ifelse(AffyID %in% TL_New1, paste(PBT, "TL_New1", sep = ","), PBT),
        PBT = ifelse(AffyID %in% TAL_New1, paste(PBT, "TAL_New1", sep = ","), PBT),
        PBT = ifelse(AffyID %in% TAL_New2, paste(PBT, "TAL_New2", sep = ","), PBT),
        PBT = ifelse(AffyID %in% TAL_New3, paste(PBT, "TAL_New3", sep = ","), PBT),
        PBT = ifelse(AffyID %in% TAL_New4, paste(PBT, "TAL_New4", sep = ","), PBT),
        PBT = ifelse(AffyID %in% DCT_New1, paste(PBT, "DCT_New1", sep = ","), PBT),
        PBT = ifelse(AffyID %in% DCT_New2, paste(PBT, "DCT_New2", sep = ","), PBT),
        PBT = ifelse(AffyID %in% DCT_New3, paste(PBT, "DCT_New3", sep = ","), PBT),
        PBT = ifelse(AffyID %in% DCT_New4, paste(PBT, "DCT_New4", sep = ","), PBT),
        PBT = ifelse(AffyID %in% CNT_New1, paste(PBT, "CNT_New1", sep = ","), PBT),
        PBT = ifelse(AffyID %in% CNT_New2, paste(PBT, "CNT_New2", sep = ","), PBT),
        PBT = ifelse(AffyID %in% CNT_New3, paste(PBT, "CNT_New3", sep = ","), PBT),
        PBT = ifelse(AffyID %in% CD_IC_New1, paste(PBT, "CD_IC_New1", sep = ","), PBT),
        PBT = ifelse(AffyID %in% CD_IC_New2, paste(PBT, "CD_IC_New2", sep = ","), PBT),
        PBT = ifelse(AffyID %in% CD_PC_New1, paste(PBT, "CD_PC_New1", sep = ","), PBT),
        PBT = ifelse(AffyID %in% CD_PC_New2, paste(PBT, "CD_PC_New2", sep = ","), PBT),
        PBT = ifelse(AffyID %in% EC_New1, paste(PBT, "EC_New1", sep = ","), PBT),
        PBT = PBT %>% str_remove("^,")
    )

# affymap219_new %>%
#     tibble() %>%
#     dplyr::filter(AffyID %in% EC_New1) %>%
#     pull(PBT)



# SAVE AN UPDATED AFFYMAP ####
affymap219 <- affymap219_new
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/"
save(affymap219, file = paste(saveDir, "affymap219_26Jul2024_PTG.RData", sep = ""))
saveDir2 <- "Z:/DATA/Datalocks/Other data/"
save(affymap219, file = paste(saveDir2, "affymap219_26Jul2024_PTG.RData", sep = ""))
