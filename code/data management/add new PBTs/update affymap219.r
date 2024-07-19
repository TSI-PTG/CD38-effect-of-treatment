# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
# Custom operators, functions, and datasets
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
# load affymap
load("Z:/DATA/Datalocks/Other data/affymap219_21Oct2019_1306_JR.RData")
# load gene lists
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/ABMR_NK_genes.RData")
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/ABMR_endothelial_genes.RData")
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/ABMR_activity_genes.RData")
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/ABMR_IFNG_genes.RData")


# DEFINE ABMR ACTIVITY PBTLIST ####
ABMRact <- genes_ABMR_activity$AffyID
ABMRifng <- genes_ABMR_IFNG$AffyID
ABMRnk <- genes_ABMR_NK$AffyID
ABMRendo <- genes_ABMR_endothelial$AffyID


# ADD PBT ANNOTATION TO affymap219 ####
affymap219 <- affymap219 %>%
    mutate(
        PBT = case_when(
            AffyID %in% ABMRact ~ paste(PBT, "AAG", sep = ","),
            AffyID %in% ABMRifng ~ paste(PBT, "IIAAG", sep = ","),
            AffyID %in% ABMRnk ~ paste(PBT, "NKAAG", sep = ","),
            AffyID %in% ABMRendo ~ paste(PBT, "AEG", sep = ","),
            TRUE ~ PBT
        )
    )

# affymap219 %>%
#     tibble() %>%
#     dplyr::filter(AffyID %in% ABMRact) %>%
#     pull(PBT)




# SAVE AN UPDATED AFFYMAP SET ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/"
save(affymap219, file = paste(saveDir, "affymap219_19Jul2024_PTG.RData", sep = ""))
