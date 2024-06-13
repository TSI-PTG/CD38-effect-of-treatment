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
genes_NK_L765 <- Hmisc::.q(
    KLRF1,
    PRF1,
    TBX21,
    GNLY,
    CD247,
    KLRD1,
    NKG7,
    GZMB,
    CX3CR1,
    SH2D1B,
    IL2RB,
    FCRL6,
    S1PR5,
    TRDC,
    FGFBP2,
    SAMD3,
    GZMA,
    NCR1,
    CD160,
    PYHIN1,
    MYBL1
)


# DEFINE THE ABMR ACTIVITY GENES BY MEAN EXPRESSION ####
genes_NK_L765 <- vienna_1208 %>%
    exprs() %>%
    as.data.frame() %>%
    mutate(means = rowMeans(.), .before = 1) %>%
    as_tibble(rownames = "AffyID") %>%
    right_join(affymap219 %>% dplyr::select(AffyID, Symb) %>% tibble(), ., by = "AffyID") %>%
    dplyr::filter(Symb %in% genes_NK_L765) %>%
    arrange(Symb, means) %>%
    dplyr::select(AffyID, Symb, means)


# SAVE THE DATA ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/"
save(genes_NK_L765, file = paste(saveDir, "NK_genes_L765.RData", sep = ""))


# EXPLORE THE DATA ####
genes_NK_L765 %>% slice_max(means, by = "Symb")
