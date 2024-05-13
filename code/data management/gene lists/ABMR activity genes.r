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


# # IQR FILTER THE DATA ####
# f1 <- function(x) (IQR(x) > 0.5)
# ff <- filterfun(f1)
# if (!exists("selected")) {
#     selected <- genefilter(vienna_1208, ff)
# }
# set00 <- vienna_1208[selected, vienna_1208$STUDY_EVALUATION_ID %nin% c(15, 18)]



# DEFINE THE MOLECULAR VARIABLES####
vars_molecular <- Hmisc::.q(
    ABMRpm, ggt0, ptcgt0, NKB, DSAST,
    TCMRt, tgt1, igt1, TCB, QCAT
)

# DEFINE THE ABMR ACTIVITY GENES ####
vars_abmr <- Hmisc::.q(
    CCL4,
    KLRD1,
    GNLY,
    CXCL11,
    CXCL9,
    PRF1,
    WARS,
    NKG7,
    GBP4,
    XCL1,
    IRF1,
    PSMB9,
    IDO,
    GZMB,
    FGFBP2,
    SH2D1B
)


# DEFINE THE ABMR ACTIVITY GENES BY MEAN EXPRESSION ####
genes_ABMR_activity <- vienna_1208 %>%
    exprs() %>%
    as.data.frame() %>%
    mutate(means = rowMeans(.), .before = 1) %>%
    as_tibble(rownames = "AffyID") %>%
    right_join(affymap219 %>% dplyr::select(AffyID, Symb) %>% tibble(), ., by = "AffyID") %>%
    dplyr::filter(Symb %in% vars_abmr) %>%
    arrange(Symb, means) %>%
    dplyr::select(AffyID, Symb, means)


# SAVE THE DATA ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/"
save(genes_ABMR_activity, file = paste(saveDir, "ABMR_activity_genes.RData", sep = ""))


# EXPLORE THE DATA ####
genes_ABMR_activity  %>% slice_max(means, by = "Symb")
