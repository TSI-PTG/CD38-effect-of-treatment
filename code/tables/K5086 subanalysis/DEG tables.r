# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(flextable) # install.packages("flextable")
library(officer) # install.packages("officer")
library(openxlsx) # install.packages("openxlsx")
library(readxl) # install.packages("readxl")
# Custom operators, functions, and datasets
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
# load affymap
load("Z:/DATA/Datalocks/Other data/affymap219_21Oct2019_1306_JR.RData")
# load limma results
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/limma_K5086.RData")
# load gene lists
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/Hinze_injury_markers.RData")
# load simple files for PBT correlations
data_scc_all_import <- read_excel(
    "Z:/MISC/Phil/AA All papers in progress/A GC papers/0000 simple XL files/Kidney 5086/MASTER COPY K5086 corr with ABMR activity genesets.xlsx",
    sheet = 1
)
data_scc_FABMR_import <- read_excel(
    "Z:/MISC/Phil/AA All papers in progress/A GC papers/0000 simple XL files/Kidney 5086/MASTER COPY K5086 corr with ABMR activity genesets Rej7AA.xlsx",
    sheet = "FABMR (redundant)"
)
data_scc_EABMR_import <- read_excel(
    "Z:/MISC/Phil/AA All papers in progress/A GC papers/0000 simple XL files/Kidney 5086/MASTER COPY K5086 corr with ABMR activity genesets Rej7AA.xlsx",
    sheet = "EABMR (redundant)"
)


# WRANGLE THE SCC DATA ####
data_scc_all <- data_scc_all_import %>%
    dplyr::select(-SYMB, -Name, -PBT) %>%
    dplyr::rename(AffyID = Affy) %>%
    dplyr::mutate_at(dplyr::vars(dplyr::contains("corr")), ~ round(., 2)) %>%
    dplyr::mutate_at(dplyr::vars(dplyr::contains("pval")), ~ formatC(., digits = 0, format = "e"))

data_scc_FABMR <- data_scc_FABMR_import %>%
    dplyr::select(-SYMB, -Name, -PBT) %>%
    dplyr::rename(AffyID = Affy) %>%
    dplyr::mutate_at(dplyr::vars(dplyr::contains("corr")), ~ round(., 2)) %>%
    dplyr::mutate_at(dplyr::vars(dplyr::contains("pval")), ~ formatC(., digits = 0, format = "e"))

data_scc_EABMR <- data_scc_EABMR_import %>%
    dplyr::select(-SYMB, -Name, -PBT) %>%
    dplyr::rename(AffyID = Affy) %>%
    dplyr::mutate_at(dplyr::vars(dplyr::contains("corr")), ~ round(., 2)) %>%
    dplyr::mutate_at(dplyr::vars(dplyr::contains("pval")), ~ formatC(., digits = 0, format = "e"))


# WRANGLE THE INJURY MARKER DATA ####
injury_markers <- genes_injury_markers %>%
    unnest(data) %>%
    drop_na(AffyID) %>%
    dplyr::filter(
        celltypename %>% str_detect(c("leukocytes"), negate = TRUE),
        cluster %>% str_detect(c("New"), negate = FALSE),
        # abs(log2FC) > 1
    ) %>%
    dplyr::select(celltypename:Symb) %>%
    nest(.by = AffyID) %>%
    mutate(
        "cellular expression" = map_chr(
            data,
            function(data) {
                data %>%
                    pull(cluster) %>%
                    paste(collapse = ", ")
            }
        )
    ) %>%
    unnest(data) %>%
    distinct(Symb, .keep_all = TRUE) %>%
    dplyr::select(AffyID, `cellular expression`)


# DEFINE NUMBER OF ROWS YOU WANT IN TABLE ####
ngene <- 40


# FORMAT TABLES TO MAKE FLEXTABLES ####
limma_tables_formatted <- limma_tables %>%
    expand_grid(direction = c("all", "increased", "decreased")) %>%
    relocate(direction, .after = design) %>%
    mutate(
        gene_tables = pmap(
            list(direction, table),
            function(direction, table) {
                if (direction == "increased") {
                    table <- table %>%
                        dplyr::filter(logFC > 0) %>%
                        arrange(p, logFC  %>% desc)
                } else if (direction == "decreased") {
                    table <- table %>%
                        dplyr::filter(logFC < 0) %>%
                        arrange(p, logFC)
                }
                table <- table %>%
                    left_join(injury_markers, by = "AffyID") %>%
                    left_join(data_scc_all, by = "AffyID") %>%
                    dplyr::select(-t, -contains("AffyID"), -contains("MMDx"), -contains("macrophage")) %>%
                    distinct(Symb, .keep_all = TRUE) %>% 
                    dplyr::slice(1:ngene) %>%
                    mutate(
                        Gene = Gene %>% str_remove("///.*"),
                        logFC = logFC %>% round(2),
                        p = case_when(
                            p < 0.0001 ~ p %>% formatC(digits = 0, format = "e"),
                            TRUE ~ p %>% formatC(digits = 4, format = "f")
                        ),
                        FDR = case_when(
                            FDR < 0.0001 ~ FDR %>% formatC(digits = 0, format = "e"),
                            TRUE ~ FDR %>% formatC(digits = 4, format = "f")
                        )
                    ) %>%
                    relocate(
                        c(
                            # `HUVEC (unstimulated)`, `HUVEC (IFNg stimulated)`,
                            #  `RPTEC (unstimulated)`,`RPTEC (IFNg stimulated)`,
                            `cellular expression`, "logFC", "p", "FDR"
                        ),
                        .after = PBT
                    ) %>%
                    dplyr::select(-`cellular expression`)
            }
        )
    )
limma_tables_formatted$gene_tables[[2]] %>% colnames()


# GLOBAL PARAMETERS FOR FLEXTABLES ####
header1 <- c(
    # "AffyID",
    "Gene\nsymbol", "Gene", "PBT",
    rep("Differential expression", 3),
    rep("Mean expression by group", 7),
    rep("Cell panel expression", 4),
    rep("Correlation with ABMR activity genesets in all K5086", 8)
)
header2 <- c(
    # "AffyID",
    "Gene\nsymbol", "Gene", "PBT",
    "logFC", "P", "FDR",
    "NR", "EABMR", "FABMR", "LABMR", "TCMR1", "TCMR2", "Minor",
    "HUVEC", "HUVEC\n(+IFNg)", " RPTEC", "RPTEC\n(+IFNg)",
    rep(c("AAG", "IIAAG", "NKAAG", "AEG"), each = 2)
)
header3 <- c(
    # "AffyID",
    "Gene\nsymbol", "Gene", "PBT",
    "logFC", "P", "FDR",
    "NR", "EABMR", "FABMR", "LABMR", "TCMR1", "TCMR2", "Minor",
    "HUVEC", "HUVEC\n(+IFNg)", " RPTEC", "RPTEC\n(+IFNg)",
    rep(c("SCC", "p"), 4)
)

cellWidths <- c(1.5, 5, 3, 1, 1, 1, rep(1, 7), rep(1, 4), rep(1, 8)) # for individual tables up or down
cellWidths %>% length()


# MAKE FORMATTED FLEXTABLES ####
flextables <- limma_tables_formatted %>%
    mutate(
        flextable = pmap(
            list(design, gene_tables, direction),
            function(design, gene_tables, direction) {
                if (direction == "all") {
                    title <- paste("Table i. Top ", ngene, " differentially expressed genes (by P-value) in ", design, sep = "")
                } else if (direction == "increased") {
                    title <- paste("Table i. Top ", ngene, " increased genes (by P-value) in ", design, sep = "")
                } else if (direction == "decreased") {
                    title <- paste("Table i. Top ", ngene, " decreased genes (by P-value) in ", design, sep = "")
                }
                gene_tables %>%
                    flextable::flextable() %>%
                    flextable::delete_part("header") %>%
                    flextable::add_header_row(top = TRUE, values = header3) %>%
                    flextable::add_header_row(top = TRUE, values = header2) %>%
                    flextable::add_header_row(top = TRUE, values = header1) %>%
                    flextable::add_header_row(top = TRUE, values = rep(title, ncol_keys(.))) %>%
                    flextable::merge_v(part = "header") %>%
                    flextable::merge_h(part = "header") %>%
                    flextable::border_remove() %>%
                    flextable::border(part = "header", border = officer::fp_border()) %>%
                    flextable::border(part = "body", border = officer::fp_border()) %>%
                    flextable::align(align = "center") %>%
                    flextable::align(align = "center", part = "header") %>%
                    flextable::font(fontname = "Arial", part = "all") %>%
                    flextable::fontsize(size = 7, part = "all") %>%
                    flextable::fontsize(size = 7, part = "footer") %>%
                    flextable::fontsize(i = 1, size = 12, part = "header") %>%
                    flextable::bold(part = "header") %>%
                    flextable::bg(bg = "white", part = "all") %>%
                    flextable::padding(padding = 0, part = "all") %>%
                    flextable::width(width = cellWidths, unit = "cm") %>%
                    flextable::width(., width = dim(.)$widths * 33 / (flextable_dim(.)$widths), unit = "cm")
            }
        )
    )
#
flextables$flextable[[3]]


# PRINT THE DATA TO POWERPOINT ####
flextables %>%
    dplyr::filter(design == "FABMRvNR", direction == "decreased") %>%
    pull(flextable) %>%
    pluck(1) %>%
    print(preview = "pptx")

flextables %>%
    dplyr::filter(design == "FABMRvNONABMR", direction == "decreased") %>%
    pull(flextable) %>%
    pluck(1) %>%
    print(preview = "pptx")

flextables %>%
    dplyr::filter(design == "EABMRvNR", direction == "decreased") %>%
    pull(flextable) %>%
    pluck(1) %>%
    print(preview = "pptx")

flextables %>%
    dplyr::filter(design == "EABMRvNONABMR", direction == "decreased") %>%
    pull(flextable) %>%
    pluck(1) %>%
    print(preview = "pptx")

flextables %>%
    dplyr::filter(design == "ABMRvNR", direction == "decreased") %>%
    pull(flextable) %>%
    pluck(1) %>%
    print(preview = "pptx")

flextables %>%
    dplyr::filter(design == "ABMRvNONABMR", direction == "decreased") %>%
    pull(flextable) %>%
    pluck(1) %>%
    print(preview = "pptx")
