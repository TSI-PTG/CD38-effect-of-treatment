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



# JOIN THE GSEA RESULTS ####
gsea <- reduce(
    list(
        felzartamab_gsea_go_k1208,
        felzartamab_gsea_msigdb_k1208,
        felzartamab_gsea_wiki_k1208,
        felzartamab_gsea_do_k1208,
        felzartamab_gsea_kegg_k1208,
        felzartamab_gsea_reactome_k1208 
    ), bind_rows
) %>%
    mutate(keep = map_dbl(gsea_tables, nrow)) %>%
    dplyr::filter(keep > 0) %>%
    dplyr::select(db, design, gsea_tables)



# WRANGLE THE GSEA TABLES ####
gsea_tables <- gsea %>%
    mutate(
        gsea_tables = map(
            gsea_tables,
            function(gsea_tables) {
                gsea_tables %>%
                    mutate(Description = Description %>% as.character()) %>%
                    dplyr::select(ID, Description, setSize, NES, p.adjust, core_enrichment)
            }
        )
    ) %>%
    unnest(gsea_tables) %>%
    arrange(p.adjust) %>%
    dplyr::filter(p.adjust < 0.001) %>% 
    nest(.by = design, gsea_table = -design)


gsea_tables$gsea_table
