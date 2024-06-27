# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(readxl) # install.packages("readxl")
# Custom operators, functions, and datasets
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
# load affymap
load("Z:/DATA/Datalocks/Other data/affymap219_21Oct2019_1306_JR.RData")
# load SCC data
simplefile_path <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/0000 simple XL files/Kidney 5086/MASTER COPY K5086 SimpleCorrAAInjRej 5AAInjNR 7AARej.xlsx"
simplefile <- read_excel(path = simplefile_path)
# load mean expression in K5086
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/mean_expression_by_probe_5086.RData")


# WRANGLE THE MEAN EXPRESSION DATA ####
means_K5086 <- mean_exprs_5086 %>%
    slice_max(mean_expression, by = "Symb")


# WRANGLE THE SIMPLE FILE DATA ####
K5086 <- simplefile %>%
    dplyr::select(
        Affy, SYMB, Name, PBT,
        `corrRej7AA4-EABMR`, `pvalRej7AA4-EABMR`,
        `corrRej7AA5-FABMR`, `pvalRej7AA5-FABMR`
        ) %>%
    dplyr::rename(AffyID = Affy, Symb = SYMB, Gene = Name) %>%
    dplyr::filter(AffyID %in% means_K5086$AffyID)



# DEFINE THE ABMR ACTIVITY GENES ####
# vars_abmr <- Hmisc::.q(
#     CCL4,
#     KLRD1,
#     GNLY,
#     CXCL11,
#     CXCL9,
#     PRF1,
#     WARS,
#     NKG7,
#     GBP4,
#     XCL1,
#     IRF1,
#     PSMB9,
#     IDO,
#     GZMB,
#     FGFBP2,
#     SH2D1B
# )


# DEFINE THE ABMR ACTIVITY GENES BY MEAN EXPRESSION ####
genes_ABMR_activity <- K5086 %>%
    slice_min(`pvalRej7AA4-EABMR`, n = 20) %>%
    dplyr::select(AffyID, Symb)


# SAVE THE DATA ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/"
save(genes_ABMR_activity, file = paste(saveDir, "ABMR_activity_genes.RData", sep = ""))


