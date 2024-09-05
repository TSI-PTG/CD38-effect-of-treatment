# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(flextable) # install.packages("flextable") 
library(officer) # install.packages("officer")
# Custom operators and functions
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
# load gsea results
loadDir <- "natmed/results/"
load(paste(loadDir, "gsea_dose.RData", sep = ""))
load(paste(loadDir, "gsea_kegg.RData", sep = ""))
load(paste(loadDir, "gsea_go.RData", sep = ""))
load(paste(loadDir, "gsea_wiki.RData", sep = ""))
load(paste(loadDir, "gsea_msigdb.RData", sep = ""))
load(paste(loadDir, "gsea_reactome.RData", sep = ""))
load(paste(loadDir, "gsea_aki.RData", sep = ""))


# JOIN THE GSEA RESULTS ####
gsea <- purrr::reduce(
    list(
        gsea_go,
        gsea_msigdb,
        gsea_wiki,
        gsea_dose,
        gsea_kegg,
        gsea_reactome,
        gsea_aki
    ), bind_rows
) %>%
    mutate(keep = map_dbl(gsea_tables, nrow)) %>%
    dplyr::filter(keep > 0) %>%
    dplyr::select(db, design, genes_gsea, gsea_tables) %>%
    mutate(gsea_tables = map(gsea_tables, mutate, Description = Description %>% as.character()))


# DEFINE INTERPREATION BASED ON PATHWAYS ####
pathway_interpretation <- gsea %>%
    dplyr::select(design, gsea_tables) %>%
    unnest(everything()) %>%
    mutate(
        interpretation = case_when(
            ID %in% c(
                "GO:0006952",
                "GO:0009605",
                "GO:0009607",
                "GO:0043207",
                "GO:0044419",
                "GO:0051707",
                "GO:0006955",
                "GO:0002376",
                "DOID:2914",
                "DOID:77",
                "DOID:1579",
                "R-HSA-168256"
            ) & design == "Baseline_vs_Week24" ~ "immune response",
            ID %in% c(
                "R-HSA-168256",
                "R-HSA-168249",
                "R-HSA-1643685",
                "R-HSA-5663205"
            ) & design == "Baseline_vs_Week52" ~ "immune-related response to injury",
            ID %in% c("AKI", "R-HSA-109582") ~ "response to injury",
            ID %in% c(
                "R-HSA-9006934",
                "R-HSA-597592",
                "R-HSA-5653656",
                "R-HSA-392499"
            ) ~ "response to injury",
            ID %in% c(
                "R-HSA-2262752",
                "R-HSA-8953897"
            ) ~ "response to injury",
        )
    ) %>%
    dplyr::select(design, ID, Description, interpretation)
pathway_interpretation %>% print(n = "all")


# WRANGLE THE GSEA TABLES ####
gsea_tables <- gsea %>%
    mutate(
        gsea_tables = pmap(
            list(design, gsea_tables),
            function(design, gsea_tables) {
                gsea_tables %>%
                    mutate(
                        design = design,
                        Description = Description %>% as.character()
                    ) %>%
                    left_join(pathway_interpretation, by = c("ID", "Description", "design")) %>%
                    dplyr::select(ID, Description, interpretation, setSize, NES, p.adjust, core_enrichment) %>%
                    distinct(ID, setSize, interpretation, .keep_all = TRUE)
            }
        )
    ) %>%
    unnest(gsea_tables) %>%
    arrange(p.adjust) %>%
    dplyr::filter(p.adjust < 0.001) %>%
    nest(.by = design, gsea_table = c(-design, -genes_gsea)) %>%
    left_join(gsea %>% dplyr::select(design, genes_gsea) %>% distinct(design, .keep_all = TRUE),
        by = "design"
    )
gsea_tables$gsea_table
gsea_tables$genes_gsea


# MAKE FLEXTABLES ####
gsea_flextables <- gsea_tables %>%
    mutate(
        flextable = pmap(
            list(design, gsea_table),
            function(design, gsea_table) {
                cellWidths <- c(2, 3, 7, 4, 1.5, 1.5, 1.5, 13) # for individual tables up or down
                if (design == "Baseline_vs_Week24") {
                    title <- paste("Table Si. Functional enrichment analysis of genes affects by felzartamab treatment between baseline and week24 (by FDR)", sep = "")
                } else if (design == "Week24_vs_Week52") {
                    title <- paste("Table Si. Functional enrichment analysis of genes affects by felzartamab treatment between week24 and week52 (by FDR)", sep = "")
                } else if (design == "Baseline_vs_Week52") {
                    title <- paste("Table Si. Functional enrichment analysis of genes affects by felzartamab treatment between baseline and week52 (by FDR)", sep = "")
                }
                gsea_table %>%
                    mutate(
                        library = db,
                        pathway = Description,
                        "core enrichment genes" = core_enrichment %>% str_replace_all("/", ", "),
                        "n genes" = setSize,
                        NES = NES %>% round(2),
                        FDR = case_when(
                            p.adjust < 0.0001 ~ p.adjust %>% formatC(digits = 0, format = "e"),
                            TRUE ~ p.adjust %>% formatC(digits = 4, format = "f")
                        ),
                    ) %>%
                    dplyr::select(library, ID, pathway, interpretation, "n genes", NES, FDR, "core enrichment genes") %>%
                    flextable::flextable() %>%
                    flextable::add_header_row(top = TRUE, values = rep(title, ncol_keys(.))) %>%
                    flextable::merge_v(part = "header") %>%
                    flextable::merge_h(part = "header") %>%
                    flextable::border_remove() %>%
                    flextable::border(border = officer::fp_border(), part = "all") %>%
                    flextable::padding(padding = 0) %>%
                    flextable::padding(j = "core enrichment genes", padding.left = 5) %>%
                    flextable::align(align = "center", part = "all") %>%
                    flextable::align(j = "core enrichment genes", align = "left", part = "body") %>%
                    flextable::font(fontname = "Arial", part = "all") %>%
                    flextable::fontsize(size = 10, part = "body") %>%
                    flextable::fontsize(size = 12, part = "header") %>%
                    flextable::bg(bg = "white", part = "all") %>%
                    flextable::width(width = cellWidths, unit = "cm")
            }
        )
    )
gsea_flextables$flextable


# SAVE THE COLLATED GSEA RESULTS ####
gsea_summary <- gsea_tables
names(gsea_summary$gsea_table) <- gsea_summary$design
saveDir <- "natmed/results/"
save(gsea_summary, file = paste(saveDir, "gsea_summary.RData", sep = ""))



# PRINT FLETABLES TO POWERPOINT ####
gsea_summary %>%
    dplyr::filter(design == "Baseline_vs_Week24") %>%
    pull(flextable) %>%
    pluck(1) %>%
    print(preview = "pptx")

gsea_summary %>%
    dplyr::filter(design == "Baseline_vs_Week52") %>%
    pull(flextable) %>%
    pluck(1) %>%
    print(preview = "pptx")
