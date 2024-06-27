# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(readxl) # install.packages("readxl")
library(flextable) # install.packages("flextable")
library(officer) # install.packages("officer")
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
    dplyr::select(`Affy Probeset ID`, `Index`, `NK cell`, CD4, CD8, `Unstim HUVEC`, `HUVEC + IFNg`) %>%
    dplyr::rename(
        AffyID_U133 = `Affy Probeset ID`,
        Symb = Index,
        NK = `NK cell`,
        `HUVEC (unstimulated)` = `Unstim HUVEC`,
        `HUVEC (IFNg stimulated)` = `HUVEC + IFNg`
    ) %>%
    slice_max(NK, by = "Symb", with_ties = FALSE)



# JOIN K5086 AND CELL PANEL DATA ####
data <- K5086 %>%
    left_join(cell_panel, by = "Symb")


# DEFINE THE ABMR ACTIVITY GENES BY MEAN EXPRESSION ####
genes_EABMR <- data %>%
    slice_min(`pvalRej7AA4-EABMR`, n = 100)

genes_ABMR_NK <- genes_EABMR %>%
    dplyr::filter(
        `HUVEC (unstimulated)` < quantile(`HUVEC (unstimulated)`, 0.5, na.rm = TRUE),
        `HUVEC (IFNg stimulated)` < quantile(`HUVEC (IFNg stimulated)`, 0.5, na.rm = TRUE),
    ) %>%
    slice_max(NK, n = 20) %>%
    relocate(AffyID_U133, .after = AffyID) %>%
    relocate(`corrRej7AA4-EABMR`, `pvalRej7AA4-EABMR`, `corrRej7AA5-FABMR`, `pvalRej7AA5-FABMR`,
        .after = last_col()
    ) %>%
    mutate_at(vars(contains("corrRej")), ~ round(., 2)) %>%
    mutate_at(vars(contains("pvalRej")), ~ formatC(., digits = 0, format = "e")) %>%
    mutate(
        Gene = Gene %>% stringr::str_remove("///.*"),
    )


# SAVE THE DATA ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/"
save(genes_ABMR_NK, file = paste(saveDir, "ABMR_NK_genes.RData", sep = ""))


# UNIVERSAL VARIABLES FOR FLEXTABLE ####
header1 <- c(
    "Gene\nsymbol", "Gene", "PBT",
    rep("Cell panel expression", 5),
    rep("Correlation with ABMR in K5086", 4)
)

header2 <- c(
    "Gene\nsymbol", "Gene", "PBT",
    "NK", "CD4", "CD8", "HUVEC\n(unstimulated)", "HUVEC\n(IFNg stimulated)",
    "SCC\nEABMR\nK5086", "p\nEABMR\nK5086",
    "SCC\nFABMR\nK5086", "p\nFABMR\nK5086"
)

title <- c("Table Si. Definition of NK cell-expressed ABMR activity genes")
cellWidths <- c(1.5, 10, 4.6, 1, 1, 1, 2, 2, rep(1.3, 4))


# MAKE FLEXTABLE ####
flextable <- genes_ABMR_NK %>%
    dplyr::select(!contains("AffyID")) %>%
    mutate(
        PBT = PBT %>%
            stringr::str_remove("RAT") %>%
            stringr::str_remove("Rej-RAT") %>%
            stringr::str_replace(",,", ",")
    ) %>%
    flextable::flextable() %>%
    flextable::delete_part("header") %>%
    flextable::add_header_row(top = TRUE, values = header2) %>%
    flextable::add_header_row(top = TRUE, values = header1) %>%
    flextable::add_header_row(top = TRUE, values = rep(title, ncol_keys(.))) %>%
    flextable::merge_v(part = "header") %>%
    flextable::merge_h(part = "header") %>%
    flextable::border_remove() %>%
    flextable::border(part = "header", border = fp_border()) %>%
    flextable::border(part = "body", border = fp_border()) %>%
    flextable::align(align = "center") %>%
    flextable::align(align = "center", part = "header") %>%
    flextable::font(fontname = "Arial", part = "all") %>%
    flextable::fontsize(size = 8, part = "all") %>%
    flextable::bold(part = "header") %>%
    flextable::bg(bg = "white", part = "all") %>%
    flextable::padding(padding = 0, part = "all") %>%
    flextable::width(width = cellWidths, unit = "cm")
#   %>%
# flextable::width(., width = dim(.)$widths * 33 / (flextable::flextable_dim(.)$widths), unit = "cm")

flextable %>% print(preview = "pptx")
