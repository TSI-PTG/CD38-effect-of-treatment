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
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/all_probes_limma_1208.RData")
# load gene lists
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/genes_NK_GEP.RData")
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/ABMR_NK_genes.RData")
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/ABMR_endothelial_genes.RData")
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/ABMR_activity_genes.RData")
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/ABMR_IFNG_genes.RData")
# load in house cell panel
atagc <- read_excel("Z:/MISC/Patrick Gauthier/R/affymap219-CELL-PANEL/backup/UPDATED 2017 ANNOTATIONS - MASTERFILE - U133 HUMAN CELL PANEL - ALL PROBESETS (nonIQR) pfhptg.xlsx")
# load SCC data
simplefile_path <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/0000 simple XL files/Kidney 5086/MASTER COPY K5086 SimpleCorrAAInjRej 5AAInjNR 7AARej.xlsx"
simplefile <- read_excel(path = simplefile_path)
# load mean expression in K1208
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/mean_expression_by_probe_1208.RData")
# load mean expression in K5086
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/mean_expression_by_probe_5086.RData")


# WRANGLE THE MEAN EXPRESSION DATA ####
means_1208 <- mean_exprs_1208 %>%
    slice_max(mean_expression, by = "Symb")

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
    dplyr::select(`Index`, `NK cell`, CD4, CD8, `Unstim HUVEC`, `HUVEC + IFNg`) %>%
    dplyr::rename(
        Symb = Index,
        NK = `NK cell`,
        `HUVEC (unstimulated)` = `Unstim HUVEC`
    )


# JOIN THE CELL PANELA AND SIMPLE FILE ####
# pheno <- cell_panel %>%
#     right_join(K5086, by = "Symb") %>%
#     relocate(Gene, PBT, .after = Symb)



# FILTER NK SELECTIVITY BY GEP ####
genes_NK_ATAGC_U133 <- genes_NK_GEP %>%
    dplyr::filter(panel == "ATAGC_U133") %>%
    pull(data) %>%
    pluck(1) %>%
    dplyr::select(AffyID, Symb)

genes_NK_KTB18_RNAseq <- genes_NK_GEP %>%
    dplyr::filter(panel == "KTB18_RNAseq") %>%
    pull(data) %>%
    pluck(1) %>%
    dplyr::select(AffyID, Symb)

genes_NK_LM22_U133 <- genes_NK_GEP %>%
    dplyr::filter(panel == "LM22_U133") %>%
    pull(data) %>%
    pluck(1) %>%
    dplyr::select(AffyID, Symb)


# FILTER THE GENE TABLES ####
gene_tables <- limma_tables %>%
    expand_grid(
        geneset = c(
            "ABMR_activity",
            "NK_ATAGC_U133", "NK_KTB18_RNAseq", "NK_LM22_U133", "ABMR_NK",
            "ABMR_endothelial", "ABMR_IFNG"
        )
    ) %>%
    mutate(
        genes = case_when(
            geneset == "ABMR_activity" ~ genes_ABMR_activity$AffyID %>% list(),
            geneset == "NK_ATAGC_U133" ~ genes_NK_ATAGC_U133$AffyID %>% list(),
            geneset == "NK_KTB18_RNAseq" ~ genes_NK_KTB18_RNAseq$AffyID %>% list(),
            geneset == "NK_LM22_U133" ~ genes_NK_LM22_U133$AffyID %>% list(),
            geneset == "ABMR_NK" ~ genes_ABMR_NK$AffyID %>% list(),
            geneset == "ABMR_endothelial" ~ genes_ABMR_endothelial$AffyID %>% list(),
            geneset == "ABMR_IFNG" ~ genes_ABMR_IFNG$AffyID %>% list()
        )
    ) %>%
    mutate(
        gene_tables = pmap(
            list(genes, table),
            function(genes, table) {
                table %>%
                    dplyr::filter(AffyID %in% genes)
            }
        )
    )
gene_tables$gene_tables


# FORMAT TABLES FOR MAKING FLEXTABLES ####
gene_flextables00 <- gene_tables %>%
    mutate(
        gene_tables = map(gene_tables, function(gene_tables) {
            gene_tables %>%
                mutate(
                    Gene = Gene %>% stringr::str_remove("///.*"),
                    `<U+0394><U+0394> FC` = `<U+0394><U+0394> FC` %>% round(2),
                    `<U+0394><U+0394> logFC` = `<U+0394><U+0394> logFC` %>% round(2),
                    `<U+0394><U+0394> p` = case_when(
                        `<U+0394><U+0394> p` < 0.0001 ~ `<U+0394><U+0394> p` %>% formatC(digits = 0, format = "e"),
                        TRUE ~ `<U+0394><U+0394> p` %>% formatC(digits = 4, format = "f")
                    ),
                    `<U+0394><U+0394> FDR` = case_when(
                        `<U+0394><U+0394> FDR` < 0.001 ~ `<U+0394><U+0394> FDR` %>% formatC(digits = 0, format = "e"),
                        TRUE ~ `<U+0394><U+0394> FDR` %>% formatC(digits = 3, format = "f")
                    )
                )
        })
    ) %>%
    dplyr::select(design, geneset, gene_tables) %>%
    nest(.by = geneset) %>%
    mutate(
        data = map(
            data,
            function(data) {
                data %>%
                    pivot_wider(names_from = design, values_from = gene_tables) %>%
                    unnest(everything(), names_repair = tidyr_legacy) %>%
                    dplyr::select(
                        -Symb1, -Symb2,
                        -Gene1, -Gene2,
                        -contains("FDR"),
                        -contains("AffyID"),
                        # -contains("Gene"),
                        -contains("PBT"),
                        -contains("MMDx"),
                        -contains("Baseline_"),
                        -contains("Week24_"),
                        -contains("Week52_"),
                    )
            }
        )
    )

# gene_flextables$data[[1]] %>%
#     pivot_wider(names_from = design, values_from = gene_tables) %>%
#     unnest(everything(), names_repair = tidyr_legacy) %>%
#     dplyr::select(
#         -Symb1, -Symb2,
#         -contains("AffyID"),
#         -contains("Gene"),
#         -contains("PBT"),
#         -contains("MMDx"),
#         -contains("Baseline_"),
#         -contains("Week24_"),
#         -contains("Week52_"),
#     )


# UNIVERSAL VARIABLES FOR FLEXTABLE ####
header2 <- c(
    "Gene\nsymbol", "Gene",
    rep("Week24 - Baseline", 5),
    rep("Week52 - Week24", 5),
    rep("Week52 - Baseline", 5)
)

header3 <- c(
    "Gene\nsymbol", "Gene",
    rep(c(
        "\u394 FC\nPlacebo\n(N=10)", "\u394 FC\nFelzartamab\n(N=10)",
        "\u394\u394 logFC\n(N=10)", "\u394\u394 FC\n(N=10)",
        "\u394\u394 p"
        # , "\u394\u394 FDR"
    ), 3)
)

cellWidths <- c(4, 16, rep(c(4, 4, 3, 3, 3), 3))



# MAKE FORMATTED FLEXTABLES ####
gene_flextables <- gene_flextables00 %>%
    mutate(
        gene_flextables = pmap(
            list(geneset, data),
            function(geneset, data) {
                colnames(data) <- LETTERS[1:ncol(data)]
                title <- paste("Table i. Fold change expression in", geneset, "genes in biopsies from felzartamab vs placebo patients")
                data %>%
                    flextable::flextable() %>%
                    flextable::delete_part("header") %>%
                    flextable::add_header_row(top = TRUE, values = header3) %>%
                    flextable::add_header_row(top = TRUE, values = header2) %>%
                    flextable::add_header_row(top = TRUE, values = rep(title, flextable::ncol_keys(.))) %>%
                    # # flextable::add_footer_row(values = footnoteText[[1]], colwidths = ncol_keys(.)) %>%
                    # flextable::merge_v(j = 1:2) %>%
                    flextable::merge_v(part = "header") %>%
                    flextable::merge_h(part = "header") %>%
                    flextable::border_remove() %>%
                    flextable::border(part = "header", border = officer::fp_border()) %>%
                    flextable::border(part = "body", border = officer::fp_border()) %>%
                    flextable::border(part = "footer", border.left = officer::fp_border(), border.right = officer::fp_border()) %>%
                    flextable::border(i = 1, part = "footer", border.bottom = officer::fp_border()) %>%
                    flextable::align(align = "center") %>%
                    flextable::align(align = "center", part = "header") %>%
                    flextable::valign(i = 3, j = c(-1, -2, -7, -12, -17), valign = "bottom", part = "header") %>%
                    flextable::font(fontname = "Arial", part = "all") %>%
                    flextable::fontsize(size = 8, part = "all") %>%
                    flextable::fontsize(size = 8, part = "footer") %>%
                    flextable::fontsize(i = 1, size = 12, part = "header") %>%
                    flextable::bold(part = "header") %>%
                    # flextable::bold(j = 1, part = "body") %>%
                    flextable::bg(bg = "white", part = "all") %>%
                    # flextable::bg(i = ~ as.numeric(`<U+0394><U+0394> p`) < 0.05, j = 2:ncol_keys(.), bg = "grey90", part = "body") %>%
                    flextable::padding(padding = 0, part = "all") %>%
                    flextable::width(width = cellWidths, unit = "cm") %>%
                    flextable::width(., width = dim(.)$widths * 33 / (flextable::flextable_dim(.)$widths), unit = "cm")
            }
        )
    )

gene_flextables$gene_flextables[[1]]


# EXPORT THE DATA AS .RData FILE ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/"
save(gene_tables, file = paste(saveDir, "gene_tables_limma_1208.RData", sep = ""))


# PRINT THE DATA TO POWERPOINT ####
gene_flextables %>%
    dplyr::filter(geneset == "ABMR_activity") %>%
    pull(gene_flextables) %>%
    pluck(1)  %>% print(preview ="pptx")

gene_flextables %>%
    dplyr::filter(geneset == "ABMR_IFNG") %>%
    pull(gene_flextables) %>%
    pluck(1) %>%
    print(preview = "pptx")

gene_flextables %>%
    dplyr::filter(geneset == "ABMR_NK") %>%
    pull(gene_flextables) %>%
    pluck(1) %>%
    print(preview = "pptx")

gene_flextables %>%
    dplyr::filter(geneset == "ABMR_endothelial") %>%
    pull(gene_flextables) %>%
    pluck(1) %>%
    print(preview = "pptx")




# SPLIT GENE TABLES BY GENESET FOR EXPORTING TABLES ####
gene_tables_ABMR_activity <- gene_tables %>% dplyr::filter(geneset == "ABMR_activity")
gene_tables_NK_ATAGC_U133 <- gene_tables %>% dplyr::filter(geneset == "NK_ATAGC_U133")
gene_tables_NK_KTB18_RNAseq <- gene_tables %>% dplyr::filter(geneset == "NK_KTB18_RNAseq")
gene_tables_NK_LM22_U133 <- gene_tables %>% dplyr::filter(geneset == "NK_LM22_U133")
gene_tables_ABMR_NK <- gene_tables %>% dplyr::filter(geneset == "ABMR_NK")
gene_tables_ABMR_Endothelial <- gene_tables %>% dplyr::filter(geneset == "ABMR_endothelial")
gene_tables_ABMR_IFNG <- gene_tables %>% dplyr::filter(geneset == "ABMR_IFNG")


# EXPORT THE DATA AS AN EXCEL SHEET ####
saveDir1 <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/output/"
# openxlsx::write.xlsx(gene_tables_ABMR_activity$gene_tables,
#     asTable = TRUE,
#     file = paste(saveDir1, "ABMR_activity_genes_limma_1208_12Jun24",
#         # Sys.Date(),
#         # format(Sys.time(), "_%I%M%p"),
#         ".xlsx",
#         sep = ""
#     )
# )
# openxlsx::write.xlsx(gene_tables_NK_ATAGC_U133$gene_tables,
#     asTable = TRUE,
#     file = paste(saveDir1, "NK_ATAGC_U133_genes_limma_1208_12Jun24",
#         # Sys.Date(),
#         # format(Sys.time(), "_%I%M%p"),
#         ".xlsx",
#         sep = ""
#     )
# )
# openxlsx::write.xlsx(gene_tables_NK_KTB18_RNAseq$gene_tables,
#     asTable = TRUE,
#     file = paste(saveDir1, "NK_KTB18_RNAseq_genes_limma_1208_12Jun24",
#         # Sys.Date(),
#         # format(Sys.time(), "_%I%M%p"),
#         ".xlsx",
#         sep = ""
#     )
# )
# openxlsx::write.xlsx(gene_tables_NK_LM22_U133$gene_tables,
#     asTable = TRUE,
#     file = paste(saveDir1, "NK_LM22_U133_genes_limma_1208_12Jun24",
#         # Sys.Date(),
#         # format(Sys.time(), "_%I%M%p"),
#         ".xlsx",
#         sep = ""
#     )
# )



# openxlsx::write.xlsx(gene_tables_NK_L765$gene_tables,
#     asTable = TRUE,
#     file = paste(saveDir1, "NK_L765_genes_limma_1208_12Jun24",
#         # Sys.Date(),
#         # format(Sys.time(), "_%I%M%p"),
#         ".xlsx",
#         sep = ""
#     )
# )
# openxlsx::write.xlsx(gene_tables_Endothelial$gene_tables,
#     asTable = TRUE,
#     file = paste(saveDir1, "Endothelial_genes_limma_1208_12Jun24",
#         # Sys.Date(),
#         # format(Sys.time(), "_%I%M%p"),
#         ".xlsx",
#         sep = ""
#     )
# )
# openxlsx::write.xlsx(gene_tables_IFNG$gene_tables,
#     asTable = TRUE,
#     file = paste(saveDir1, "IFNG_genes_limma_1208_12Jun24",
#         # Sys.Date(),
#         # format(Sys.time(), "_%I%M%p"),
#         ".xlsx",
#         sep = ""
#     )
# )