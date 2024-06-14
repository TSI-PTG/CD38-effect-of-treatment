# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(flextable) # install.packages("flextable")
library(officer) # install.packages("officer")
library(openxlsx) # install.packages("openxlsx")
# Custom operators, functions, and datasets
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
# load affymap
load("Z:/DATA/Datalocks/Other data/affymap219_21Oct2019_1306_JR.RData")
# load limma results
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/IQR_filtered_probes_unique_genes_baseline_corrected_cortex_corrected_limma_1208.RData")



# FORMAT TABLES TO MAKE FLEXTABLES ####
limma_tables <- limma_tables %>%
    mutate(
        gene_tables = map(
            table,
            function(table) {
                table %>% dplyr::select(
                    -AffyID, -cortex, 
                    -contains("placebo FC"), -contains("felz FC"),
                    -contains("MMDx")
                )
            }
        )
    )
# limma_tables$gene_tables[[1]]  %>% colnames


# GLOBAL PARAMETERS FOR FLEXTABLES ####
header1 <- c(
    # "AffyID",
    "Gene\nsymbol", "Gene", "PBT", "\u394\u394FC", "\u394\u394P", "\u394\u394FDR",
    rep("Mean expression by group", 4)
)
header2 <- c(
    # "AffyID",
    "Gene\nsymbol", "Gene", "PBT", "\u394\u394FC", "\u394\u394P", "\u394\u394FDR",
    rep("Placebo", 2), rep("Felzartamab", 2)
)

cellWidths <- c(1.5, 5, 3, 1, 1, 1, 1, rep(1, 4)) # for individual tables up or down
cellWidths  %>% length

limma_tables$gene_tables[[1]] %>%
    flextable::flextable() %>%
    flextable::delete_part("header")


# MAKE FORMATTED FLEXTABLES ####
flextables <- limma_tables %>%
    mutate(
        flextables = pmap(
            list(design, gene_tables),
            function(design, gene_tables) {
                colnames(gene_tables) <- LETTERS[1:ncol(gene_tables)]
                if (design == "Baseline_vs_Week24") {
                    title <- paste("Table i. Top 20 differentially expressed genes between baseline and week24 in biopsies from placebo and Felzartamab treated patients (by P-value)", sep = "")
                    header3 <- c(
                        "Symb", "Gene", "PBT", "FC", "P", "FDR",
                        "Baseline\n(N=10)", "Week24\n(N=10)", "Baseline\n(N=10)", "Week24\n(N=10)"
                    )
                } else if (design == "Week24_vs_Week52") {
                    title <- paste("Table i. Top 20 differentially expressed genes between baseline and week24 in biopsies from placebo and Felzartamab treated patients (by P-value)", sep = "")

                    header3 <- c(
                        "Symb", "Gene", "PBT", "FC", "P", "FDR",
                        "Week24\n(N=10)", "Week52\n(N=10)", "Week24\n(N=10)", "Week52\n(N=10)"
                    )
                } else if (design == "Baseline_vs_Week52") {
                    title <- paste("Table i. Top 20 differentially expressed genes between baseline and week24 in biopsies from placebo and Felzartamab treated patients (by P-value)", sep = "")
                    header3 <- c(
                        "Symb", "Gene", "PBT", "FC", "P", "FDR",
                        "Baseline\n(N=10)", "Week52\n(N=10)", "Baseline\n(N=10)", "Week52\n(N=10)"
                    )
                }
                gene_tables %>%
                    flextable::flextable() %>%
                    flextable::delete_part("header") %>%
                    # flextable::add_header_row(top = TRUE, values = header3) %>%
                    # flextable::add_header_row(top = TRUE, values = header2) %>%
                    # flextable::add_header_row(top = TRUE, values = header1) %>%
                    # flextable::add_header_row(top = TRUE, values = rep(title, ncol_keys(.))) %>%
                    # # flextable::add_footer_row(values = footnoteText[[1]], colwidths = ncol_keys(.)) %>%
                    # flextable::merge_v(j = 1:2) %>%
                    flextable::merge_v(part = "header") %>%
                    flextable::merge_h(part = "header") %>%
                    flextable::border_remove() %>%
                    flextable::border(part = "header", border = fp_border()) %>%
                    flextable::border(part = "body", border = fp_border()) %>%
                    # flextable::border(part = "footer", border.left = fp_border(), border.right = fp_border()) %>%
                    # flextable::border(i = 1, part = "footer", border.bottom = fp_border()) %>%
                    flextable::align(align = "center") %>%
                    flextable::align(align = "center", part = "header") %>%
                    # flextable::valign(i = 3, j = c(-1, -2, -7, -12, -17), valign = "bottom", part = "header") %>%
                    flextable::font(fontname = "Arial", part = "all") %>%
                    flextable::fontsize(size = 8, part = "all") %>%
                    flextable::fontsize(size = 8, part = "footer") %>%
                    flextable::fontsize(i = 1, size = 12, part = "header") %>%
                    flextable::bold(part = "header") %>%
                    # flextable::bold(j = 1, part = "body") %>%
                    flextable::bg(bg = "white", part = "all") %>%
                    # flextable::bg(i = ~ as.numeric(`<U+0394><U+0394> p`) < 0.05, j = 2:ncol_keys(.), bg = "grey90", part = "body") %>%
                    flextable::padding(padding = 0, part = "all") %>%
                    # flextable::width(width = cellWidths, unit = "cm") %>%
                    flextable::width(., width = dim(.)$widths * 33 / (flextable_dim(.)$widths), unit = "cm")
            }
        )
    )

flextables$flextables[[1]]



# PRINT THE DATA TO POWERPOINT ####
# gene_flextables$gene_flextables %>% print(preview = "pptx")



# EXPORT THE DATA AS .RData FILE ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/"
save(gene_tables, file = paste(saveDir, "gene_tables_limma_1208.RData", sep = ""))



# SPLIT GENE TABLES BY GENESET FOR EXPORTING TABLES ####
gene_tables_ABMR_activity <- gene_tables %>% dplyr::filter(geneset == "ABMR_activity")
gene_tables_NK_ATAGC_U133 <- gene_tables %>% dplyr::filter(geneset == "NK_ATAGC_U133")
gene_tables_NK_KTB18_RNAseq <- gene_tables %>% dplyr::filter(geneset == "NK_KTB18_RNAseq")
gene_tables_NK_LM22_U133 <- gene_tables %>% dplyr::filter(geneset == "NK_LM22_U133")
gene_tables_NK_L765 <- gene_tables %>% dplyr::filter(geneset == "NK_L765")
gene_tables_Endothelial <- gene_tables %>% dplyr::filter(geneset == "Endothelial")
gene_tables_IFNG <- gene_tables %>% dplyr::filter(geneset == "IFNG")


# EXPORT THE DATA AS AN EXCEL SHEET ####
saveDir1 <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/output/"
openxlsx::write.xlsx(gene_tables_ABMR_activity$gene_tables,
    asTable = TRUE,
    file = paste(saveDir1, "ABMR_activity_genes_limma_1208_12Jun24",
        # Sys.Date(),
        # format(Sys.time(), "_%I%M%p"),
        ".xlsx",
        sep = ""
    )
)
openxlsx::write.xlsx(gene_tables_NK_ATAGC_U133$gene_tables,
    asTable = TRUE,
    file = paste(saveDir1, "NK_ATAGC_U133_genes_limma_1208_12Jun24",
        # Sys.Date(),
        # format(Sys.time(), "_%I%M%p"),
        ".xlsx",
        sep = ""
    )
)
openxlsx::write.xlsx(gene_tables_NK_KTB18_RNAseq$gene_tables,
    asTable = TRUE,
    file = paste(saveDir1, "NK_KTB18_RNAseq_genes_limma_1208_12Jun24",
        # Sys.Date(),
        # format(Sys.time(), "_%I%M%p"),
        ".xlsx",
        sep = ""
    )
)
openxlsx::write.xlsx(gene_tables_NK_LM22_U133$gene_tables,
    asTable = TRUE,
    file = paste(saveDir1, "NK_LM22_U133_genes_limma_1208_12Jun24",
        # Sys.Date(),
        # format(Sys.time(), "_%I%M%p"),
        ".xlsx",
        sep = ""
    )
)













openxlsx::write.xlsx(gene_tables_NK_L765$gene_tables,
    asTable = TRUE,
    file = paste(saveDir1, "NK_L765_genes_limma_1208_12Jun24",
        # Sys.Date(),
        # format(Sys.time(), "_%I%M%p"),
        ".xlsx",
        sep = ""
    )
)
openxlsx::write.xlsx(gene_tables_Endothelial$gene_tables,
    asTable = TRUE,
    file = paste(saveDir1, "Endothelial_genes_limma_1208_12Jun24",
        # Sys.Date(),
        # format(Sys.time(), "_%I%M%p"),
        ".xlsx",
        sep = ""
    )
)
openxlsx::write.xlsx(gene_tables_IFNG$gene_tables,
    asTable = TRUE,
    file = paste(saveDir1, "IFNG_genes_limma_1208_12Jun24",
        # Sys.Date(),
        # format(Sys.time(), "_%I%M%p"),
        ".xlsx",
        sep = ""
    )
)
