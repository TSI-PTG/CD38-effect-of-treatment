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
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/IQR_filtered_probes_unique_genes_baseline_corrected_cortex_corrected_limma_1208.RData")
# load gene lists
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/Hinze_injury_markers.RData")
# load K4502 injury simple file
simplefile_path <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/0000 simple XL files/Kidney 4502/MASTER COPY K4502 Injset SimpleCorrAAInjRej 5AAInj 7AARej.xlsx"
simplefile <- read_excel(path = simplefile_path, sheet = "simpleCorrAAInjRejInjset")
# load in house cell panel
atagc <- readxl::read_excel("Z:/MISC/Patrick Gauthier/R/affymap219-CELL-PANEL/backup/UPDATED 2017 ANNOTATIONS - MASTERFILE - U133 HUMAN CELL PANEL - ALL PROBESETS (nonIQR) pfhptg.xlsx")


# DEFINE CELL STATES OF INTEREST ####
cell_states <- c(
    "PT-New1",
    "PT-New2",
    "PT-New3",
    "PT-New4",
    "TAL-New1",
    "TAL-New2",
    "TAL-New3",
    "TAL-New4",
    "DCT-New1",
    "DCT-New2",
    "DCT-New3",
    "DCT-New4",
    "CNT-New1",
    "CNT-New2",
    "CNT-New3",
    "CD-PC-New1",
    "CD-PC-New2",
    "CD-IC-New1",
    "CD-IC-New2"
)


# WRANGLE THE INJURY MARKER DATA ####
injury_markers <- genes_injury_markers %>%
    unnest(data) %>%
    nest(.by = c(celltypename, cluster), injury_genes = c(-"celltypename", -"cluster")) %>%
    mutate(
        injury_genes = map(
            injury_genes,
            function(injury_genes) {
                injury_genes %>%
                    drop_na(AffyID) %>%
                    arrange(p, log2FC %>% desc()) %>%
                    # dplyr::slice(1:20) %>%
                    dplyr::select(AffyID)
            }
        )
    )
injury_markers$injury_genes[[1]]



# WRANGLE THE SIMPLE FILE DATA ####
K4502 <- simplefile %>%
    dplyr::select(Affy, corrInjPCA1, corrInjPCA2, corrInjPCA3) %>%
    dplyr::rename(AffyID = Affy)


# JOIN DE AND INJURY MARKER DATA ####
data <- injury_markers %>%
    expand_grid(design = limma_tables$design) %>%
    left_join(limma_tables, by = "design") %>%
    dplyr::select(-toptable) %>%
    expand_grid(direction = c("all", "increased", "decreased"))


# FORMAT TABLES TO MAKE FLEXTABLES ####
tables <- data %>%
    mutate(
        gene_tables = pmap(
            list(injury_genes, table, direction),
            function(injury_genes, table, direction) {
                genes <- injury_genes %>% pull(AffyID)
                colnames(table) <- table %>%
                    colnames() %>%
                    str_remove_all("\u394 |\u394")
                if (direction == "inallreased") {
                    table <- table
                } else if (direction == "increased") {
                    table <- table %>% dplyr::filter(logFC > 0)
                } else if (direction == "decreased") {
                    table <- table %>% dplyr::filter(logFC < 0)
                }
                df <- table %>%
                    dplyr::filter(AffyID %in% genes) %>%
                    # left_join(injury_genes, by = "AffyID") %>%
                    arrange(p) %>%
                    distinct(Symb, .keep_all = TRUE) %>%
                    # left_join(K4502, by = "AffyID")  %>%
                    # dplyr::slice_min(p, n = 20) %>%
                    dplyr::select(
                        -AffyID, -cortex,
                        -contains("placebo FC"), -contains("felz FC"),
                        -contains("MMDx")
                    ) %>%
                    mutate(
                        Gene = Gene %>% str_remove("///.*"),
                        logFC = logFC %>% round(2),
                        FC = FC %>% round(2),
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
                        c("logFC", "FC", "p", "FDR"),
                        .after = PBT
                    )
            }
        )
    )
tables$gene_tables


# GLOBAL PARAMETERS FOR FLEXTABLES ####
header1 <- c(
    # "AffyID",
    "Gene\nsymbol", "Gene", "PBT",
    "\u394\u394 logFC", "\u394\u394 FC", "\u394\u394 P", "\u394\u394 FDR",
    rep("Mean expression by group", 4)
)
header2 <- c(
    # "AffyID",
    "Gene\nsymbol", "Gene", "PBT",
    "\u394\u394 logFC", "\u394\u394 FC", "\u394\u394 P", "\u394\u394 FDR",
    rep("Placebo", 2), rep("Felzartamab", 2)
)

cellWidths <- c(1.5, 5, 3, 1, 1, 1, 1, rep(1, 4)) # for individual tables up or down
cellWidths %>% length()

# limma_tables$gene_tables[[1]] %>%
#     flextable::flextable() %>%
#     flextable::delete_part("header")


# MAKE FORMATTED FLEXTABLES ####
flextables <- tables %>%
    mutate(
        flextables = pmap(
            list(design, gene_tables, direction),
            function(design, gene_tables, direction) {
                # colnames(gene_tables) <- LETTERS[1:ncol(gene_tables)]
                if (design == "Baseline_vs_Week24") {
                    title <- paste(
                        "Table i. Top 20 ",
                        direction,
                        " injury genes between baseline and week24 in biopsies from placebo and Felzartamab treated patients (by P-value)",
                        sep = ""
                    )
                    header3 <- c(
                        "Gene\nsymbol", "Gene", "PBT", "\u394\u394 logFC", "\u394\u394 FC", "\u394\u394 P", "\u394\u394 FDR",
                        "Baseline\n(N=10)", "Week24\n(N=10)", "Baseline\n(N=10)", "Week24\n(N=10)"
                    )
                } else if (design == "Week24_vs_Week52") {
                    title <- paste("Table i. Top 20 ",
                        direction,
                        " injury genes between week24 and week52 in biopsies from placebo and Felzartamab treated patients (by P-value)",
                        sep = ""
                    )

                    header3 <- c(
                        "Gene\nsymbol", "Gene", "PBT", "\u394\u394 logFC", "\u394\u394 FC", "\u394\u394 P", "\u394\u394 FDR",
                        "Week24\n(N=10)", "Week52\n(N=10)", "Week24\n(N=10)", "Week52\n(N=10)"
                    )
                } else if (design == "Baseline_vs_Week52") {
                    title <- paste("Table i. Top 20 ",
                        direction,
                        " injury genes between baseline and week52 in biopsies from placebo and Felzartamab treated patients (by P-value)",
                        sep = ""
                    )
                    header3 <- c(
                        "Gene\nsymbol", "Gene", "PBT", "\u394\u394 logFC", "\u394\u394 FC", "\u394\u394 P", "\u394\u394 FDR",
                        "Baseline\n(N=10)", "Week52\n(N=10)", "Baseline\n(N=10)", "Week52\n(N=10)"
                    )
                }
                gene_tables %>%
                    dplyr::slice_min(p, n = 20) %>%
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
                    flextable::width(width = cellWidths, unit = "cm") %>%
                    flextable::width(., width = dim(.)$widths * 33 / (flextable_dim(.)$widths), unit = "cm")
            }
        )
    )

flextables %>%
    dplyr::filter(
        direction == "decreased",
        design == "Baseline_vs_Week52",
        cluster == "TAL-New4"
    ) %>%
    pull(flextables)



# EXPORT THE DATA AS .RData FILE ####
tables$gene_tables[[1]]

gene_tables <- tables
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/"
save(gene_tables, file = paste(saveDir, "injury_gene_tables_limma_1208.RData", sep = ""))



# PRINT THE DATA TO POWERPOINT ####
flextables %>%
    dplyr::filter(design == "Baseline_vs_Week52", direction == "decreased") %>%
    pull(flextables) %>%
    pluck(1)
#  %>%
# print(preview = "pptx")


flextables %>%
    dplyr::filter(design == "Baseline_vs_Week24", direction == "decreased") %>%
    pull(flextables) %>%
    pluck(1)
