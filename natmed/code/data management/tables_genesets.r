# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(flextable) # install.packages("flextable")
library(officer) # install.packages("officer")
library(openxlsx) # install.packages("openxlsx")
library(readxl) # install.packages("readxl")
# Custom operators and functions
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
# load affymap
load("natmed/data/affymap219.RData")
# load DE results
load("natmed/results/deg.RData")
# load gene lists
load("natmed/data/genesets_abmr.RData")


# WRANGLE THE GENE TABLES ####
tables_genesets_abmr <- deg %>%
    tidyr::crossing(gene_tables) %>%
    mutate(
        table_deg = pmap(
            list(table, genes),
            function(table, genes) {
                table %>%
                    dplyr::filter(Symb %in% genes)
            }
        )
    ) %>%
    dplyr::select(geneset, design, genes, table_genes, table_deg)


# SAVE THE TABLES ####
saveDir <- "natmed/data/"
save(tables_genesets_abmr, file = paste0(saveDir, "tables_genesets_abmr.RData"))


# FORMAT TABLES FOR MAKING FLEXTABLES ####
gene_flextables00 <- tables_genesets_abmr %>%
    mutate(
        gene_tables = map(table_deg, function(table_deg) {
            table_deg %>%
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
gene_flextables00$data[[1]]


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
