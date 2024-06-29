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
load(paste(loadDir, "felzartamab_gsea_DO_baseline_corrected_cortex_corrected_k1208.RData", sep = ""))


# WRANGLE THE GSEA RESULTS ####
gsea_reactome <- felzartamab_gsea_reactome_k1208 %>%
    mutate(db = "reactome", .before = 1) %>%
    dplyr::rename(gsea = gsea_reactome, gsea_tables = gsea_reactome_tables) %>%
    dplyr::select(db, design, gsea_tables)

gsea_msigdb <- felzartamab_gsea_msigdb_k1208 %>%
    mutate(db = "msigdb", .before = 1) %>%
    dplyr::rename(gsea = gsea_msigdb, gsea_tables = gsea_msigdb_tables) %>%
    dplyr::select(db, design, gsea_tables)

gsea_wiki <- felzartamab_gsea_wiki_k1208 %>%
    mutate(db = "wiki", .before = 1) %>%
    dplyr::rename(gsea = gsea_wiki, gsea_tables = gsea_wiki_tables) %>%
    dplyr::select(db, design, gsea_tables)

gsea_dose <- felzartamab_gsea_do_k1208 %>%
    mutate(db = "dose", .before = 1) %>%
    dplyr::rename(gsea = gsea_do, gsea_tables = gsea_do_tables) %>%
    dplyr::select(db, design, gsea_tables)




# JOIN THE GSEA RESULTS ####



gsea_reactome$gsea_tables[[1]]

gsea_msigdb$gsea_tables[[1]]
