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
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/All_probes_cortex_corrected_limma_1208.RData")
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

limma_tables$table[[1]] %>%
    dplyr::filter(AffyID == "11756632_a_at")


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
            geneset == "ABMR_activity" ~ genes_ABMR_activity$Symb %>% list(),
            geneset == "NK_ATAGC_U133" ~ genes_NK_ATAGC_U133$Symb %>% list(),
            geneset == "NK_KTB18_RNAseq" ~ genes_NK_KTB18_RNAseq$Symb %>% list(),
            geneset == "NK_LM22_U133" ~ genes_NK_LM22_U133$Symb %>% list(),
            geneset == "ABMR_NK" ~ genes_ABMR_NK$Symb %>% list(),
            geneset == "ABMR_endothelial" ~ genes_ABMR_endothelial$Symb %>% list(),
            geneset == "ABMR_IFNG" ~ genes_ABMR_IFNG$Symb %>% list()
        )
    ) %>%
    mutate(
        gene_tables = pmap(
            list(genes, table),
            function(genes, table) {
                table %>%
                    dplyr::filter(Symb %in% genes)
            }
        )
    )


genes_ABMR_activity %>%
    dplyr::filter(Symb %nin% gene_tables$gene_tables[[1]]$Symb)
gene_tables$gene_tables[[1]]$Symb


# FORMAT TABLES FOR MAKING FLEXTABLES ####
gene_flextables00 <- gene_tables %>%
    mutate(
        gene_tables = map(gene_tables, function(gene_tables) {
            gene_tables %>%
                mutate(
                    Gene = Gene %>% stringr::str_remove("///.*"),
                    plogFC = plogFC %>% round(2),
                    flogFC = flogFC %>% round(2),
                    logFC = logFC %>% round(2),
                    p = case_when(
                        p < 0.0001 ~ p %>% formatC(digits = 0, format = "e"),
                        TRUE ~ p %>% formatC(digits = 4, format = "f")
                    ),
                    FDR = case_when(
                        FDR < 0.001 ~ FDR %>% formatC(digits = 0, format = "e"),
                        TRUE ~ FDR %>% formatC(digits = 3, format = "f")
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
                        Symb, Gene,
                        plogFC, flogFC, logFC, p,
                        plogFC1, flogFC1, logFC1, p1,
                        plogFC2, flogFC2, logFC2, p2,
                    )
            }
        )
    )
gene_flextables00$data[[6]]


# UNIVERSAL VARIABLES FOR FLEXTABLE ####
header2 <- c(
    "Gene\nsymbol", "Gene",
    rep("Baseline - Week24", 4),
    rep("Week24 - Week52", 4),
    rep("Baseline - Week52", 4)
)
header2 %>% length()

header3 <- c(
    "Gene\nsymbol", "Gene",
    rep(c(
        "\u394\nlogFC\nPlacebo\n(N=10)",
        "\u394\nlogFC\nFelzartamab\n(N=10)",
        "\u394\u394\nlogFC\n(N=20)",
        "\u394\u394\np"
    ), 3)
)
header3 %>% length()

cellWidths <- c(4, 16, rep(c(4, 4, 3, 3), 3))
cellWidths %>% length()




# MAKE FORMATTED FLEXTABLES ####
gene_flextables <- gene_flextables00 %>%
    mutate(
        gene_flextables = pmap(
            list(geneset, data),
            function(geneset, data) {
                title <- paste("Table i. Fold change expression in", geneset, "genes in biopsies from felzartamab vs placebo patients")
                colmeans <- c(
                    Symb = "column means",
                    Gene = "column means",
                    data %>%
                        dplyr::select(-Symb, -Gene) %>%
                        mutate_all(as.numeric) %>%
                        colMeans() %>%
                        round(2)
                ) %>%
                    bind_rows() %>%
                    mutate_at(vars(contains("logFC")), ~ as.numeric(.)) %>%
                    mutate_at(vars(p, p1, p2), ~NA)
                data_table <- data %>%
                    bind_rows(colmeans) %>%
                    mutate_at(vars(contains("logFC")), ~ round(., 2))
                data_table %>%
                    flextable::flextable() %>%
                    flextable::delete_part("header") %>%
                    flextable::add_header_row(top = TRUE, values = header3) %>%
                    flextable::add_header_row(top = TRUE, values = header2) %>%
                    flextable::add_header_row(top = TRUE, values = rep(title, flextable::ncol_keys(.))) %>%
                    # # flextable::add_footer_row(values = footnoteText[[1]], colwidths = ncol_keys(.)) %>%
                    # flextable::merge_v(j = 1:2) %>%
                    flextable::merge_v(part = "header") %>%
                    flextable::merge_h(part = "header") %>%
                    flextable::merge_at(j = 1:2, i = flextable::nrow_part(., "body"), part = "body") %>%
                    flextable::border_remove() %>%
                    flextable::border(part = "header", border = officer::fp_border()) %>%
                    flextable::border(part = "body", border = officer::fp_border()) %>%
                    flextable::border(part = "footer", border.left = officer::fp_border(), border.right = officer::fp_border()) %>%
                    flextable::border(i = 1, part = "footer", border.bottom = officer::fp_border()) %>%
                    flextable::align(align = "center") %>%
                    flextable::align(align = "center", part = "header") %>%
                    flextable::align(j = 1:2, i = flextable::nrow_part(., "body"), align = "right", part = "body") %>%
                    # flextable::valign(i = 3, j = c(-1, -2, -7, -12, -17), valign = "bottom", part = "header") %>%
                    flextable::font(fontname = "Arial", part = "all") %>%
                    flextable::fontsize(size = 8, part = "all") %>%
                    flextable::fontsize(size = 8, part = "footer") %>%
                    flextable::fontsize(i = 1, size = 12, part = "header") %>%
                    flextable::bold(part = "header") %>%
                    flextable::bold(j = 1:2, i = flextable::nrow_part(., "body"), part = "body") %>%
                    flextable::italic(j = 1:2, i = flextable::nrow_part(., "body"), part = "body") %>%
                    flextable::bg(bg = "white", part = "all") %>%
                    # flextable::bg(i = ~ as.numeric(`<U+0394><U+0394> p`) < 0.05, j = 2:ncol_keys(.), bg = "grey90", part = "body") %>%
                    flextable::padding(padding = 0, part = "all") %>%
                    flextable::padding(j = 1:2, i = flextable::nrow_part(., "body"), padding.right = 5, part = "body") %>%
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
    pluck(1) %>%
    print(preview = "pptx")

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



gene_flextables00$data[[7]]


fuckyou <- c(
    Symb = "column means",
    Gene = "column means",
    gene_flextables00 %>%
        dplyr::filter(geneset == "ABMR_activity") %>%
        pull(data) %>%
        pluck(1) %>%
        dplyr::select(-Symb, -Gene) %>%
        mutate_all(as.numeric) %>%
        colMeans()
)


cunt <- tibble(!!!fuckyou)

cunt %>% mutate_at(vars(contains("logFC")), ~ as.numeric(.))


# SPLIT GENE TABLES BY GENESET FOR EXPORTING TABLES ####
# gene_tables_ABMR_activity <- gene_tables %>% dplyr::filter(geneset == "ABMR_activity")
# gene_tables_NK_ATAGC_U133 <- gene_tables %>% dplyr::filter(geneset == "NK_ATAGC_U133")
# gene_tables_NK_KTB18_RNAseq <- gene_tables %>% dplyr::filter(geneset == "NK_KTB18_RNAseq")
# gene_tables_NK_LM22_U133 <- gene_tables %>% dplyr::filter(geneset == "NK_LM22_U133")
# gene_tables_ABMR_NK <- gene_tables %>% dplyr::filter(geneset == "ABMR_NK")
# gene_tables_ABMR_Endothelial <- gene_tables %>% dplyr::filter(geneset == "ABMR_endothelial")
# gene_tables_ABMR_IFNG <- gene_tables %>% dplyr::filter(geneset == "ABMR_IFNG")


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
