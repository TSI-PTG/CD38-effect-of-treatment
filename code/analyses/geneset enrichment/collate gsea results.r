# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(flextable) # install.packages("flextable") #for simple table outputs
library(officer) # install.packages("officer")
# library(clusterProfiler) # pak::pak("YuLab-SMU/clusterProfiler")
# library(msigdbr) # install.packages("msigdbr")
# library(circlize) # install.packages("circlize") pak::pak("jokergoo/circlize")
# Bioconductor libraries
# library(Biobase) # BiocManager::install("Biobase")
# library(biobroom) # BiocManager::install("biobroom")
# library(org.Hs.eg.db) # BiocManager::install("org.Hs.eg.db")
# library(pRoloc) # BiocManager::install("pRoloc")
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
# load gsea results
loadDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/"
load(paste(loadDir, "felzartamab_gsea_msigdb_baseline_corrected_cortex_corrected_k1208.RData", sep = ""))
load(paste(loadDir, "felzartamab_gsea_wiki_baseline_corrected_cortex_corrected_k1208.RData", sep = ""))
load(paste(loadDir, "felzartamab_gsea_reactome_baseline_corrected_cortex_corrected_k1208.RData", sep = ""))
load(paste(loadDir, "felzartamab_gsea_dose_baseline_corrected_cortex_corrected_k1208.RData", sep = ""))
load(paste(loadDir, "felzartamab_gsea_kegg_baseline_corrected_cortex_corrected_k1208.RData", sep = ""))
load(paste(loadDir, "felzartamab_gsea_go_baseline_corrected_cortex_corrected_k1208.RData", sep = ""))
load(paste(loadDir, "felzartamab_gsea_mesh_baseline_corrected_cortex_corrected_k1208.RData", sep = ""))
load(paste(loadDir, "felzartamab_gsea_aki_baseline_corrected_cortex_corrected_k1208.RData", sep = ""))



# JOIN THE GSEA RESULTS ####
gsea <- purrr::reduce(
    list(
        felzartamab_gsea_go_k1208,
        felzartamab_gsea_msigdb_k1208,
        felzartamab_gsea_wiki_k1208,
        felzartamab_gsea_do_k1208,
        felzartamab_gsea_kegg_k1208,
        felzartamab_gsea_reactome_k1208,
        felzartamab_gsea_aki_k1208
    ), bind_rows
) %>%
    mutate(keep = map_dbl(gsea_tables, nrow)) %>%
    dplyr::filter(keep > 0) %>%
    dplyr::select(db, design, genes_gsea, gsea_tables) %>%
    mutate(gsea_tables = map(gsea_tables, mutate, Description = Description %>% as.character()))


# DEFINE INTERPREATION BASED ON PATHWAYS ####
pathway_interpretation <- gsea %>%
    dplyr::select(gsea_tables) %>%
    unnest(everything()) %>%
    dplyr::select(ID, Description) %>%
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
                "R-HSA-168256",
                "DOID:2914",
                "DOID:77",
                "DOID:1579",
                "R-HSA-168249",
                "R-HSA-168256",
                "R-HSA-1643685",
                "R-HSA-5663205"
            ) ~ "immune response",
            ID == "AKI" ~ "AKI",
            ID %in% c(
                "R-HSA-9006934",
                "R-HSA-597592",
                "R-HSA-5653656",
                "R-HSA-392499"
            ) ~ "cellular metabolism",
            ID %in% c(
                "R-HSA-109582",
                "R-HSA-2262752",
                "R-HSA-8953897"
            ) ~ "cellular stress response",
        )
    )


# WRANGLE THE GSEA TABLES ####
gsea_tables <- gsea %>%
    mutate(
        gsea_tables = map(
            gsea_tables,
            function(gsea_tables) {
                gsea_tables %>%
                    mutate(Description = Description %>% as.character()) %>%
                    left_join(pathway_interpretation) %>%
                    dplyr::select(ID, Description, interpretation, setSize, NES, p.adjust, core_enrichment)
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
felzartamab_gsea_k1208 <- gsea_tables
names(felzartamab_gsea_k1208$gsea_table) <- felzartamab_gsea_k1208$design
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/"
save(felzartamab_gsea_k1208, file = paste(saveDir, "felzartamab_gsea_collated_k1208.RData", sep = ""))






# PRINT FLETABLES TO POWERPOINT ####
gsea_flextables %>%
    dplyr::filter(design == "Baseline_vs_Week24") %>%
    pull(flextable) %>%
    pluck(1) %>%
    print(preview = "pptx")

gsea_flextables %>%
    dplyr::filter(design == "Baseline_vs_Week52") %>%
    pull(flextable) %>%
    pluck(1) %>%
    print(preview = "pptx")
