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
# load in house cell panel
atagc <- read_excel("Z:/MISC/Patrick Gauthier/R/affymap219-CELL-PANEL/backup/UPDATED 2017 ANNOTATIONS - MASTERFILE - U133 HUMAN CELL PANEL - ALL PROBESETS (nonIQR) pfhptg.xlsx")



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


# WRANGLE THE CELL PANEL DATA ####
cell_panel <- atagc %>%
    dplyr::select(`Affy Probeset ID`, `Index`, `NK cell`, CD4, CD8, `Unstim HUVEC`) %>%
    dplyr::rename(
        AffyID133 = `Affy Probeset ID`,
        Symb = Index,
        NK = `NK cell`,
        `HUVEC (unstimulated)` = `Unstim HUVEC`
    ) %>%
    slice_max(`HUVEC (unstimulated)`, by = "Symb", with_ties = FALSE)



# DEFINE THE ABMR ACTIVITY GENES ####
# genes_endothelial <- Hmisc::.q(
#     LYPD5,
#     RASIP1,
#     ICAM2,
#     CDH5,
#     HSPA12B,
#     GNG11,
#     GJD3,
#     ADGRL4,
#     ECSCR,
#     RAPGEF5,
#     FGR,
#     MMRN2,
#     ACSL5,
#     ROBO4,
#     ERG,
#     NOS3,
#     CDH13,
#     PECAM1
# )


# JOIN K5086 AND CELL PANEL DATA ####
data <- K5086 %>%
    left_join(cell_panel, by = "Symb")


# DEFINE THE ABMR ACTIVITY GENES BY MEAN EXPRESSION ####
genes_FABMR <- data %>%
    slice_min(`pvalRej7AA5-FABMR`, n = 100)


genes_ABMR_endothelial <- genes_FABMR %>%
    dplyr::filter(`HUVEC (unstimulated)` > NK)  %>% 
    slice_max(`HUVEC (unstimulated)`, n = 20)



# SAVE THE DATA ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/"
# save(genes_ABMR_endothelial, file = paste(saveDir, "ABMR_endothelial_genes.RData", sep = ""))
