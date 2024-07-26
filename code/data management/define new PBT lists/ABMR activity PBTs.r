# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
# Custom operators, functions, and datasets
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
# load gene lists
# load gene lists
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/genes_NK_GEP.RData")
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/ABMR_NK_genes.RData")
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/ABMR_endothelial_genes.RData")
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/ABMR_activity_genes.RData")
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/ABMR_IFNG_genes.RData")
# load pbt lists
load("Z:/DATA/Datalocks/Other data/PBTlist219_14Dec2020_JR.RData")


# DEFINE ABMR ACTIVITY PBTLIST ####
PBTlist219_ABMRactivity <- list(
    AAG = genes_ABMR_activity$AffyID,
    IIAAG = genes_ABMR_IFNG$AffyID,
    NKAAG = genes_ABMR_NK$AffyID,
    AEG = genes_ABMR_endothelial$AffyID
)


# SAVE THE PBT LIST ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/"
save(PBTlist219_ABMRactivity, file = paste(saveDir, "PBTlist219_ABMRactivity_26Jul2024_PTG.RData", sep = ""))