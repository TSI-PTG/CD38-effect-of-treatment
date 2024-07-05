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
# load K5086 rejection simple file
simplefile_path <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/0000 simple XL files/Kidney 5086/MASTER COPY K5086 SimpleCorrAAInjRej 5AAInjNR 7AARej.xlsx"
simplefile_rejection <- readxl::read_excel(path = simplefile_path)
# load K4502 injury simple file
simplefile_path <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/0000 simple XL files/Kidney 4502/MASTER COPY K4502 Injset SimpleCorrAAInjRej 5AAInj 7AARej.xlsx"
simplefile_injury <- read_excel(path = simplefile_path, sheet = "simpleCorrAAInjRejInjset")
# load N604 cfdna simple file
simplefile_path <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/0000 simple XL files/Kidney N604/MASTER COPY SimpleCorr 6AARej N604.xlsx"
simplefile_cfdna <- read_excel(path = simplefile_path, sheet = "simpleCorrAARej")
# load mean expression in K5086
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/mean_expression_by_probe_5086.RData")
# load in house cell panel
atagc <- read_excel("Z:/MISC/Patrick Gauthier/R/affymap219-CELL-PANEL/backup/UPDATED 2017 ANNOTATIONS - MASTERFILE - U133 HUMAN CELL PANEL - ALL PROBESETS (nonIQR) pfhptg.xlsx")



# WRANGLE THE MEAN EXPRESSION DATA ####
means_K5086 <- mean_exprs_5086 %>%
    dplyr::slice_max(mean_expression, by = "Symb")


# WRANGLE THE REJECTION SIMPLE FILE DATA ####
scc_rejection <- simplefile_rejection %>%
    dplyr::select(
        Affy, SYMB, Name, PBT,
        `corrRej7AA4-EABMR`, `pvalRej7AA4-EABMR`,
        `corrRej7AA5-FABMR`, `pvalRej7AA5-FABMR`
    ) %>%
    dplyr::rename(AffyID = Affy, Symb = SYMB, Gene = Name) %>%
    dplyr::filter(AffyID %in% means_K5086$AffyID)


# WRANGLE THE INJURY SIMPLE FILE DATA ####
scc_injury <- simplefile_injury %>%
    dplyr::select(
        Affy, SYMB, Name, PBT,
        `corrRej7AA4-EABMR`, `pvalRej7AA4-EABMR`,
        `corrRej7AA5-FABMR`, `pvalRej7AA5-FABMR`
    ) %>%
    dplyr::rename(AffyID = Affy, Symb = SYMB, Gene = Name) %>%
    dplyr::filter(AffyID %in% means_K5086$AffyID)



# WRANGLE THE CFDNA  SIMPLE FILE DATA ####
scc_cfdna <- simplefile_cfdna %>%
    dplyr::select(
        Affy, SYMB, Name, PBT,
        "corrcfDNA", "pvalcfDNA",
        "corrQuant", "pvalQuant"
    ) %>%
    dplyr::rename(AffyID = Affy, Symb = SYMB, Gene = Name) %>%
    dplyr::filter(AffyID %in% means_K5086$AffyID)


# JOIN THE SCC DATA ####
scc <- left_join(scc_rejection, scc_cfdna, by = c("AffyID", "Symb", "Gene", "PBT"))


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
    dplyr::slice_max(`HUVEC (IFNg stimulated)`, by = "Symb", with_ties = FALSE)


# JOIN SCC AND CELL PANEL DATA ####
data <- scc %>% dplyr::left_join(cell_panel, by = "Symb")


# DEFINE THE ABMR ACTIVITY GENES BY MEAN EXPRESSION ####
genes_EABMR <- data %>%
    dplyr::filter(`corrRej7AA4-EABMR` >= 0.2) %>%
    dplyr::slice_min(`pvalRej7AA4-EABMR`, n = 100)

genes_ABMR_IFNG <- genes_EABMR %>%
    dplyr::filter(
        # `HUVEC (unstimulated)` < quantile(`HUVEC (unstimulated)`, 0.5, na.rm = TRUE),
        # NK < quantile(NK, 0.5, na.rm = TRUE),
        `HUVEC (IFNg stimulated)` > (NK * 5),
        # `HUVEC (IFNg stimulated)` > quantile(`HUVEC (unstimulated)`, 0.75, na.rm = TRUE),
        `HUVEC (IFNg stimulated)` > (`HUVEC (unstimulated)` * 5),
    ) %>%
    dplyr::slice_max(`HUVEC (IFNg stimulated)`, n = 20) %>%
    dplyr::relocate(AffyID_U133, .after = AffyID) %>%
    dplyr::relocate(
        `corrRej7AA4-EABMR`, `pvalRej7AA4-EABMR`, `corrRej7AA5-FABMR`, `pvalRej7AA5-FABMR`,
        corrcfDNA, pvalcfDNA, corrQuant, pvalQuant,
        .after = dplyr::last_col()
    ) %>%
    dplyr::mutate_at(dplyr::vars(dplyr::contains("corr")), ~ round(., 2)) %>%
    dplyr::mutate_at(dplyr::vars(dplyr::contains("pval")), ~ formatC(., digits = 0, format = "e")) %>%
    dplyr::mutate(Gene = Gene %>% stringr::str_remove("///.*"),)


# SAVE THE DATA ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/"
save(genes_ABMR_IFNG, file = paste(saveDir, "ABMR_IFNG_genes.RData", sep = ""))



# UNIVERSAL VARIABLES FOR FLEXTABLE ####
header1 <- c(
    "Gene\nsymbol", "Gene", "PBT",
    rep("Cell panel expression", 5),
    rep("Correlation with ABMR in K5086", 4),
    rep("Correlation with dd-cfDNA in N604", 4)
)

header2 <- c(
    "Gene\nsymbol", "Gene", "PBT",
    "NK", "CD4", "CD8", "HUVEC\n(unstimulated)", "HUVEC\n(IFNg stimulated)",
    "SCC\nEABMR\nK5086", "p\nEABMR\nK5086",
    "SCC\nFABMR\nK5086", "p\nFABMR\nK5086",
    "SCC\ndd-cfDNA (%)\nN604", "p\ndd-cfDNA (%)\nN604",
    "SCC\ndd-cfDNA (cp/mL)\nN604", "p\ndd-cfDNA (cp/mL)\nN604"
)

title <- c("Table Si. Definition of IFNG-inducible ABMR activity genes")
cellWidths <- c(1.5, 9, 5.4, 1, 1, 1, 2, 2, rep(1.3, 8))


# MAKE FLEXTABLE ####
flextable <- genes_ABMR_IFNG %>%
    dplyr::select(!dplyr::contains("AffyID")) %>%
    dplyr::mutate(
        PBT = PBT %>%
            stringr::str_remove(",RAT") %>%
            stringr::str_remove(",Rej-RAT") %>%
            stringr::str_replace(",,", ",")
    ) %>%
    flextable::flextable() %>%
    flextable::delete_part("header") %>%
    flextable::add_header_row(top = TRUE, values = header2) %>%
    flextable::add_header_row(top = TRUE, values = header1) %>%
    flextable::add_header_row(top = TRUE, values = rep(title, flextable::ncol_keys(.))) %>%
    flextable::merge_v(part = "header") %>%
    flextable::merge_h(part = "header") %>%
    flextable::border_remove() %>%
    flextable::border(part = "header", border = officer::fp_border()) %>%
    flextable::border(part = "body", border = officer::fp_border()) %>%
    flextable::align(align = "center") %>%
    flextable::align(align = "center", part = "header") %>%
    flextable::font(fontname = "Arial", part = "all") %>%
    flextable::fontsize(size = 8, part = "all") %>%
    flextable::bold(part = "header") %>%
    flextable::bg(bg = "white", part = "all") %>%
    flextable::padding(padding = 0, part = "all") %>%
    flextable::width(width = cellWidths, unit = "cm")

# flextable %>% print()
flextable %>% print(preview = "pptx")
