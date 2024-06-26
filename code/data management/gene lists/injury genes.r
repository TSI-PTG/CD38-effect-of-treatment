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


# DEFINE GENERIC INJURY MARKER GENES ####
genes_markers <- Hmisc::.q(
    LCN2,
    HAVCR1,
    IGFBP7,
    TIMP2,
    IL18,
    CHI3L1,
    SPP1,
    S100A9
)


# DEFINE PROXIMAL TUBULE INJURY GENES ####
genes_PT <- Hmisc::.q(
    SPARC,
    VCAM1,
    IFITM2,
    IFITM3,
    VIM,
    "HLA-A",
    BTG1,
    SERPINA1,
    RPL5,
    RPL4,
    MYO5B,
    AKAP12,
    ERO1A,
    KLF6,
    NQ01,
    MT1F,
    MT1HL1,
    FTL,
    FTH1,
    GPX3,
    SCL7A13,
    SCL22A24,
    SLC5A11,
    SLC22A7,
    SLC34A1,
    SLC7A7,
    SLC5A2
)


# DEFINE THICK ASCENDING LIMB INJURY GENES ####
genes_TAL <- Hmisc::.q(
    MET,
    ITGA2,
    ITGB1,
    LAMC2,
    CD44,
    TMP1,
    NNMT,
    SERPINA1,
    "HLA-A",
    IFITM3,
    IFITM2,
    RPL30,
    RPL5,
    RPL4,
    ERO1A,
    SLC2A1,
    VEGFA,
    ALDOB,
    FTH1,
    MT1E,
    FTL,
    MT1G,
    GPX3,
    CLDN10,
    SLC12A1,
    FGF9,
    CLDN16,
    CASR,
    ROBO2,
    PAPP2,
    NOS1
)


# DEFINE THE INJURY GENES BY MEAN EXPRESSION ####
genes_injury_PT <- vienna_1208 %>%
    exprs() %>%
    as.data.frame() %>%
    mutate(means = rowMeans(.), .before = 1) %>%
    as_tibble(rownames = "AffyID") %>%
    right_join(affymap219 %>% dplyr::select(AffyID, Symb) %>% tibble(), ., by = "AffyID") %>%
    dplyr::filter(Symb %in% genes_PT) %>%
    arrange(Symb, means) %>%
    dplyr::select(AffyID, Symb, means) %>%
    mutate(celltype = "PT", .before = 1)

genes_injury_TAL <- vienna_1208 %>%
    exprs() %>%
    as.data.frame() %>%
    mutate(means = rowMeans(.), .before = 1) %>%
    as_tibble(rownames = "AffyID") %>%
    right_join(affymap219 %>% dplyr::select(AffyID, Symb) %>% tibble(), ., by = "AffyID") %>%
    dplyr::filter(Symb %in% genes_TAL) %>%
    arrange(Symb, means) %>%
    dplyr::select(AffyID, Symb, means) %>%
    mutate(celltype = "TAL", .before = 1)

genes_injury <- bind_rows(genes_injury_PT, genes_injury_TAL)


# SAVE THE DATA ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/"
save(genes_injury, file = paste(saveDir, "injury_genes.RData", sep = ""))


# EXPLORE THE DATA ####
genes_injury %>% slice_max(means, by = c("celltype", "Symb"))
