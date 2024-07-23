# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(flextable) # install.packages("flextable") #for simple table outputs
library(officer) # install.packages("officer")
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
# load GSEA results
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/felzartamab_gsea_baseline_corrected_cortex_corrected_k1208.RData")
# load simplified GO annotations
source("C:/R/CD38-effect-of-treatment/code/data management/functional enrichment annotation/very simplified GO annotation.r")



# FORMAT TABLES FOR GENESET ENRICHMENT ANALYSES (GSEA) ####
gsea_go_tables <- felzartamab_gsea_k1208 %>%
    dplyr::select(design, gsea_go_tables) %>%
    mutate(
        gsea_go_tables = map(
            gsea_go_tables,
            function(gsea_go_tables) {
                gsea_go_tables %>%
                    # dplyr::filter(Description %nin% c("biological_process", "cellular process", "positive regulation of biological process")) %>%
                    mutate(
                        "Simplified description" = case_when(
                            Description %>% str_detect(immune_response) ~ "immune response",
                            Description %>% str_detect(infection_response) ~ "response to infection",
                            Description %>% str_detect(response_to_stimulus) ~ "response to exogenous/endogenous stimulus",
                            Description %>% str_detect(inflammation) ~ "inflammation",
                            Description %>% str_detect(injury) ~ "injury response",
                            Description %>% str_detect(cell_signalling_and_RNA_transcription) ~ "cell signalling and RNA transcription",
                            Description %>% str_detect(cellular_development_and_metabolism) ~ "cell development,\nmobilization,\nand metabolism",
                            TRUE ~ Description
                        ) %>% factor(levels = GO_annotation_levels_truncated),
                        "Core enrichment genes" = core_enrichment %>% str_replace_all("/", ", "),
                        NES = NES %>% round(1),
                        FDR = case_when(
                            p.adjust < 0.0001 ~ p.adjust %>% formatC(digits = 0, format = "e"),
                            TRUE ~ p.adjust %>% formatC(digits = 4, format = "f")
                        ),
                    ) %>%
                    mutate() %>%
                    dplyr::rename(
                        "Go term" = Description,
                        "N genes" = setSize
                    ) %>%
                    dplyr::select(ID, "Go term", "Simplified description", "N genes", NES, FDR, "Core enrichment genes")
            }
        )
    )
gsea_go_tables$gsea_go_tables[[1]]


# MAKE FLEXTABLES SUMMARIZING GENESET ENRICHMENT ANALYSES (GSEA) ####
gsea_go_flextables <- gsea_go_tables %>%
    mutate(
        gsea_go_flextables = pmap(
            list(design, gsea_go_tables),
            function(design, gsea_go_tables) {
                cellWidths <- c(2, 6.1, 2.8, 1.2, 1.2, 1.3, 12)
                title <- paste("Table Si. Geneset enrichment analysis of genes affected by Felzartamab treatment between", design) %>%
                    stringr::str_replace("_vs_", " and ")
                gsea_go_tables %>%
                    dplyr::slice_min(FDR, n = 20) %>%
                    flextable::flextable() %>%
                    flextable::add_header_row(values = rep(title, ncol_keys(.))) %>%
                    flextable::merge_h(part = "header") %>%
                    flextable::border_remove() %>%
                    flextable::border(border = officer::fp_border(), part = "all") %>%
                    flextable::align(align = "center", part = "all") %>%
                    flextable::fontsize(size = 8, part = "all") %>%
                    flextable::fontsize(i = 1, size = 10, part = "header") %>%
                    flextable::bold(part = "header") %>%
                    flextable::bg(bg = "white", part = "all") %>%
                    flextable::width(width = cellWidths, unit = "cm") %>%
                    flextable::padding(padding = 0)
                    # %>%
                    # flextable::width(., width = dim(.)$widths * 26 / (flextable_dim(.)$widths), unit = "cm")
            }
        )
    )
# gsea_go_flextables$gsea_go_flextables[[1]]


# PRINT FLEXTABLES TO POWERPOINT ####
gsea_go_flextables$gsea_go_flextables[[1]] %>% print(preview = "pptx")
gsea_go_flextables$gsea_go_flextables[[2]] %>% print(preview = "pptx")
gsea_go_flextables$gsea_go_flextables[[3]] %>% print(preview = "pptx")
