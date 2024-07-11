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


# FORMAT TABLES TO MAKE FLEXTABLES ####
limma_tables$table[[1]]
limma_tables <- limma_tables %>%
    mutate(
        gene_tables = map(
            table,
            function(table) {
                # colnames(table) <- table %>%
                #     colnames() %>%
                #     str_remove_all("\u394 |\u394")
                table %>%
                    left_join(injury_markers, by = "AffyID") %>%
                    dplyr::select(-t, -contains("AffyID"), -contains("MMDx")) %>%
                    dplyr::slice(1:20) %>%
                    mutate(
                        Gene = Gene %>% str_remove("///.*"),
                        plogFC = plogFC %>% round(2),
                        flogFC = flogFC %>% round(2),
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
                            `cellular expression`, "plogFC", "flogFC", "logFC", "p", "FDR"
                        ),
                        .after = PBT
                    )
            }
        )
    )
limma_tables$gene_tables[[1]] %>% colnames()


# GLOBAL PARAMETERS FOR FLEXTABLES ####
header1 <- c(
    # "AffyID",
    "Gene\nsymbol", "Gene", "PBT",
    "AKI-induced cell states",
    rep("Differential expression", 5),
    rep("Mean expression by group", 6),
    rep("Cell panel expression", 4)
)
header2 <- c(
    # "AffyID",
    "Gene\nsymbol", "Gene", "PBT",
    "AKI-induced cell states",
    "\u394\nplacebo\nlogFC", "\u394\nfelzartamab\nlogFC",
    "\u394\u394\nlogFC", "\u394\u394\nP", "\u394\u394\nFDR",
    rep("Placebo", 3), rep("Felzartamab", 3),
    "HUVEC", "HUVEC\n(+IFNg)", " RPTEC", "RPTEC\n(+IFNg)"
)
header3 <- c(
    "Gene\nsymbol", "Gene", "PBT",
    "AKI-induced cell states",
    "\u394\nplacebo\nlogFC", "\u394\nfelzartamab\nlogFC",
    "\u394\u394\nlogFC", "\u394\u394\nP", "\u394\u394\nFDR",
    "Baseline\n(N=10)", "Week24\n(N=10)", "Week52\n(N=10)", "Baseline\n(N=10)", "Week24\n(N=10)", "Week52\n(N=10)",
    "HUVEC", "HUVEC\n(+IFNg)", " RPTEC", "RPTEC\n(+IFNg)"
)

cellWidths <- c(1.5, 5, 3, 1.5, 1, 1, 1, 1, 1, rep(1, 6), rep(1, 4)) # for individual tables up or down
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
                if (design == "Baseline_vs_Week24") {
                    title <- paste("Table i. Top 20 differentially expressed genes between baseline and week24 in biopsies from placebo and Felzartamab treated patients (by P-value)", sep = "")
                } else if (design == "Week24_vs_Week52") {
                    title <- paste("Table i. Top 20 differentially expressed genes between week24 and week52 in biopsies from placebo and Felzartamab treated patients (by P-value)", sep = "")
                } else if (design == "Baseline_vs_Week52") {
                    title <- paste("Table i. Top 20 differentially expressed genes between baseline and week52 in biopsies from placebo and Felzartamab treated patients (by P-value)", sep = "")
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
flextables$flextables[[3]] %>% print(preview = "pptx")
