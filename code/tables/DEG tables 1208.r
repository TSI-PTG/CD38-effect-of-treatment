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
# load gene lists
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/Hinze_injury_markers.RData")



# WRANGLE THE INJURY MARKER DATA ####
# genes_injury_markers %>%
#     unnest(data) %>%
#     nest(.by = c(celltypename, cluster)) %>%
#     print(n = "all")

injury_markers <- genes_injury_markers %>%
    unnest(data) %>%
    drop_na(AffyID) %>%
    dplyr::filter(
        celltypename %>% str_detect(c("leukocytes"), negate = TRUE),
        cluster %>% str_detect(c("New"), negate = TRUE),
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
                    left_join(injury_markers, by = "AffyID")%>%
                    dplyr::select(
                        -contains("AffyID"), 
                        -cortex, -FC,
                        -contains("placebo FC"), -contains("felz FC"),
                        -contains("MMDx")
                    )  %>% 
                    dplyr::slice(1:20) %>%
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
    "Gene\nsymbol", "Gene", "PBT", "Cellular expression\nin AKI",
    "\u394\u394 logFC", "\u394\u394 P", "\u394\u394 FDR",
    rep("Mean expression by group", 4)
)
header2 <- c(
    # "AffyID",
    "Gene\nsymbol", "Gene", "PBT", "Cellular expression\nin AKI",
    "\u394\u394 logFC", "\u394\u394 P", "\u394\u394 FDR",
    rep("Placebo", 2), rep("Felzartamab", 2)
)

cellWidths <- c(1.5, 5, 3, 1, 1, 1, 1, rep(1, 4)) # for individual tables up or down
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
                    title <- paste("Table i. Top 20 differentially expressed genes between baseline and week24 in biopsies from placebo and Felzartamab treated patients (by P-value)", sep = "")
                    header3 <- c(
                        "Gene\nsymbol", "Gene", "PBT", "Cellular expression\nin AKI", "\u394\u394 logFC", "\u394\u394 P", "\u394\u394 FDR",
                        "Baseline\n(N=10)", "Week24\n(N=10)", "Baseline\n(N=10)", "Week24\n(N=10)"
                    )
                } else if (design == "Week24_vs_Week52") {
                    title <- paste("Table i. Top 20 differentially expressed genes between week24 and week52 in biopsies from placebo and Felzartamab treated patients (by P-value)", sep = "")

                    header3 <- c(
                        "Gene\nsymbol", "Gene", "PBT","Cellular expression\nin AKI", "\u394\u394 logFC",  "\u394\u394 P", "\u394\u394 FDR",
                        "Week24\n(N=10)", "Week52\n(N=10)", "Week24\n(N=10)", "Week52\n(N=10)"
                    )
                } else if (design == "Baseline_vs_Week52") {
                    title <- paste("Table i. Top 20 differentially expressed genes between baseline and week52 in biopsies from placebo and Felzartamab treated patients (by P-value)", sep = "")
                    header3 <- c(
                        "Gene\nsymbol", "Gene", "PBT", "Cellular expression\nin AKI", "\u394\u394 logFC", "\u394\u394 P", "\u394\u394 FDR",
                        "Baseline\n(N=10)", "Week52\n(N=10)", "Baseline\n(N=10)", "Week52\n(N=10)"
                    )
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
                    flextable::border(part = "header", border = fp_border()) %>%
                    flextable::border(part = "body", border = fp_border()) %>%
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

flextables$flextables[[3]]



# PRINT THE DATA TO POWERPOINT ####
flextables$flextables[[1]] %>% print(preview = "pptx")
