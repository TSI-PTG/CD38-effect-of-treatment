# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
# Custom operators, functions, and datasets
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
pbt.generate <- function(data_CEL, reference_set, pbtlist) {
    for (i in 1:length(names(pbtlist))) {
        Biobase::pData(data_CEL)[names(pbtlist)[i]] <- NA
        conavg <- apply(Biobase::exprs(reference_set)[pbtlist[[i]], ], 1, mean)
        smallset <- data_CEL[pbtlist[[i]], ]
        dog <- sweep(Biobase::exprs(smallset), 1, conavg)
        dog <- apply(dog, 2, mean)
        sel <- which(names(Biobase::pData(data_CEL)) == names(pbtlist)[i])
        Biobase::pData(data_CEL)[, sel] <- dog
    }
    return(data_CEL)
}
# load reference set
load("Z:/DATA/Datalocks/Kidney/K5086Mar13_2024_JR.RData")
# load nephrectomies
load("Z:/Genome-Archive/RefData/KidneyReports2021REDCap/data/Conset.RData") # 4 Controls used for PBT calculations
# load gene lists
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/genes_NK_GEP.RData")
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/ABMR_NK_genes.RData")
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/ABMR_endothelial_genes.RData")
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/ABMR_activity_genes.RData")
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/ABMR_IFNG_genes.RData")
# load pbt lists
load("Z:/DATA/Datalocks/Other data/PBTlist219_14Dec2020_JR.RData")


# DEFINE ABMR ACTIVITY PBTLIST ####
ABMRactivity_genesets <- list(
    ABMRact = genes_ABMR_activity$AffyID,
    ABMRifng = genes_ABMR_IFNG$AffyID, 
    ABMRnk = genes_ABMR_NK$AffyID, 
    ABMRendo = genes_ABMR_endothelial$AffyID
)


# DEFINE SET ####
set <- K5086

# GENERATE PBT SCORES ####
set <- pbt.generate(set, Conset, ABMRactivity_genesets)


# REDFINE THE SET ###
K5086 <- set


# SAVE THE SET ####
save(K5086, file = "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/K5086Jul12_2024_PTG.RData")
