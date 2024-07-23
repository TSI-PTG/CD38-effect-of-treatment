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
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/ONLY_FELZ_IQR_filtered_probes_unique_genes_baseline_corrected_cortex_corrected_limma_1208.RData")
# load gene lists
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/Hinze_injury_markers.RData")
# load in house cell panel
atagc <- readxl::read_excel("Z:/MISC/Patrick Gauthier/R/affymap219-CELL-PANEL/backup/UPDATED 2017 ANNOTATIONS - MASTERFILE - U133 HUMAN CELL PANEL - ALL PROBESETS (nonIQR) pfhptg.xlsx")



# WRANGLE THE INJURY MARKER DATA ####
# genes_injury_markers %>%
#     unnest(data) %>%
#     nest(.by = c(celltypename, cluster)) %>%
#     print(n = "all")

injury_markers <- genes_injury_markers %>%
    unnest(data) %>%
    drop_na(AffyID) %>%
    dplyr::filter(
        # celltypename %>% str_detect("leukocytes", negate = TRUE),
        cluster %>% str_detect("New"),
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


# WRANGLE THE CELL PANEL DATA ####
cell_panel <- atagc %>%
    dplyr::slice_max(`Control Nephr`, by = "Index", with_ties = FALSE) %>%
    dplyr::rename(
        AffyID_U133 = `Affy Probeset ID`,
        Symb = Index,
        NK = `NK cell`,
        `HUVEC (unstimulated)` = `Unstim HUVEC`,
        `HUVEC (IFNg stimulated)` = `HUVEC + IFNg`
    ) %>%
    dplyr::select(Symb, NK, CD4, CD8, `HUVEC (unstimulated)`, `HUVEC (IFNg stimulated)`)



# FORMAT TABLES TO MAKE FLEXTABLES ####
limma_tables <- limma_tables %>%
    mutate(
        gene_tables = map(
            table,
            function(table) {
                colnames(table) <- table %>%
                    colnames() %>%
                    str_remove_all("\u394 |\u394")
                table %>%
                    left_join(injury_markers, by = "AffyID") %>%
                    dplyr::select(
                        -contains("AffyID"),
                        -contains("placebo FC"), -contains("felzartamab FC"),
                        -contains("MMDx")
                    ) %>%
                    left_join(cell_panel, by = "Symb") %>%
                    # dplyr::slice(1:100) %>%
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
                        c(`cellular expression`, "logFC", "p", "FDR"),
                        .after = PBT
                    )
            }
        )
    )
limma_tables$gene_tables[[3]] %>% colnames()


# GLOBAL PARAMETERS FOR FLEXTABLES ####
header1 <- c(
    # "AffyID",
    "Gene\nsymbol", "Gene",
    # rep("annotation", 2),
    "PBT", "Marker gene in\ninjured kidney cells",
    # "\u394 felzartamab\nlogFC", "\u394 felzartamab\nP-value", "\u394 felzartamab\nFDR",
    rep("Temporal effect in felzartamab patients", 3),
    rep("Mean expression in felzartamab patients", 3),
    rep("Mean expression by cell type", 5)
)
header2 <- c(
    # "AffyID",
    "Gene\nsymbol", "Gene",
    # rep("annotation", 2),
    "PBT", "Marker gene in\ninjured kidney cells",
    # "\u394 felzartamab\nlogFC", "\u394 felzartamab\nP-value", "\u394 felzartamab\nFDR",
    rep("Temporal effect in felzartamab patients", 3),
    rep("Mean expression in felzartamab patients", 3),
    rep("Mean expression by cell type", 5)
)

cellWidths <- c(2.2, 7, 4, 3.5, rep(1.8, 3), rep(1.3, 3), rep(1.2, 3), rep(2, 2)) # for individual tables up or down
cellWidths %>% length()

# limma_tables$gene_tables[[1]] %>%
#     flextable::flextable() %>%
#     flextable::delete_part("header")


# MAKE FORMATTED FLEXTABLES ####
flextables <- limma_tables %>%
    mutate(
        flextables = pmap(
            list(design, gene_tables),
            function(design, gene_tables) {
                # colnames(gene_tables) <- LETTERS[1:ncol(gene_tables)]
                if (design == "Baseline_vs_Week24") {
                    title <- paste("Table i. Top 30 differentially expressed genes between baseline and week24 in felzartamab treated patients (by P-value)", sep = "")
                    header3 <- c(
                        "Gene\nsymbol", "Gene", "PBT", "Marker gene in\ninjured kidney cells",
                        "\u394 felzartamab\nlogFC", "\u394 felzartamab\nP-value", "\u394 felzartamab\nFDR",
                        "Baseline\n(N=10)", "Week24\n(N=10)","Week52\n(N=10)",
                        "NK", "CD4", "CD8", "HUVEC\n(unstimulated)", "HUVEC\n(IFNg stimulated)"
                    )
                } else if (design == "Week24_vs_Week52") {
                    title <- paste("Table i. Top 30 differentially expressed genes between week24 and week52 in biopsies from placebo and felzartamab treated patients (by P-value)", sep = "")

                    header3 <- c(
                        "Gene\nsymbol", "Gene", "PBT", "Marker gene in\ninjured kidney cells",
                        "\u394 felzartamab\nlogFC", "\u394 felzartamab\nP-value", "\u394 felzartamab\nFDR",
                        "Baseline\n(N=10)", "Week24\n(N=10)", "Week52\n(N=10)",
                        "NK", "CD4", "CD8", "HUVEC\n(unstimulated)", "HUVEC\n(IFNg stimulated)"
                    )
                } else if (design == "Baseline_vs_Week52") {
                    title <- paste("Table i. Top 30 differentially expressed genes between baseline and week52 in biopsies from placebo and felzartamab treated patients (by P-value)", sep = "")
                    header3 <- c(
                        "Gene\nsymbol", "Gene", "PBT", "Marker gene in\ninjured kidney cells",
                        "\u394 felzartamab\nlogFC", "\u394 felzartamab\nP-value", "\u394 felzartamab\nFDR",
                        "Baseline\n(N=10)",  "Week24\n(N=10)", "Week52\n(N=10)",
                        "NK", "CD4", "CD8", "HUVEC\n(unstimulated)", "HUVEC\n(IFNg stimulated)"
                    )
                }
                gene_tables %>%
                    dplyr::mutate(
                        PBT = PBT %>%
                            stringr::str_remove(",RAT") %>%
                            stringr::str_remove("Rej-RAT") %>%
                            stringr::str_replace(",,", ",")
                    ) %>%
                    slice_min(p, n = 30, with_ties = FALSE) %>%
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
                    flextable::fontsize(size = 8, part = "all") %>%
                    flextable::fontsize(size = 8, part = "footer") %>%
                    flextable::fontsize(i = 1, size = 12, part = "header") %>%
                    flextable::bold(part = "header") %>%
                    flextable::bg(bg = "white", part = "all") %>%
                    flextable::padding(padding = 0, part = "all") %>%
                    flextable::width(width = cellWidths, unit = "cm")
                # %>%
                # flextable::width(., width = dim(.)$widths * 33 / (flextable_dim(.)$widths), unit = "cm")
            }
        )
    )

flextables$flextables[[3]]



# PRINT THE DATA TO POWERPOINT ####
flextables$flextables[[3]] %>% print(preview = "pptx")
