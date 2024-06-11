# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
# Bioconductor libraries
library(Biobase) # BiocManager::install("Biobase")
library(genefilter) # BiocManager::install("genefilter")

# Custom operators, functions, and datasets
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
# load affymap
load("Z:/DATA/Datalocks/Other data/affymap219_21Oct2019_1306_JR.RData")
# load reference set
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0Georg Felz CD38 Vienna/G_Rstuff/data/vienna_1208_6Mar24.RData")



# DEFINE THE ABMR ACTIVITY GENES ####
genes_endothelial <- Hmisc::.q(
    LYPD5,
    RASIP1,
    ICAM2,
    CDH5,
    HSPA12B,
    GNG11,
    GJD3,
    ADGRL4,
    ECSCR,
    RAPGEF5,
    FGR,
    MMRN2,
    ACSL5,
    ROBO4,
    ERG,
    NOS3,
    CDH13,
    PECAM1
)


# DEFINE THE ABMR ACTIVITY GENES BY MEAN EXPRESSION ####
genes_ABMR_endothelial <- vienna_1208 %>%
    exprs() %>%
    as.data.frame() %>%
    mutate(means = rowMeans(.), .before = 1) %>%
    as_tibble(rownames = "AffyID") %>%
    right_join(affymap219 %>% dplyr::select(AffyID, Symb) %>% tibble(), ., by = "AffyID") %>%
    dplyr::filter(Symb %in% genes_endothelial) %>%
    arrange(Symb, means) %>%
    dplyr::select(AffyID, Symb, means)


# SAVE THE DATA ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/"
save(genes_ABMR_endothelial, file = paste(saveDir, "ABMR_endothelial_genes.RData", sep = ""))


# EXPLORE THE DATA ####
genes_ABMR_endothelial  %>% slice_max(means, by = "Symb")
