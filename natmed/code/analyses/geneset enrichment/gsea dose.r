# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(flextable) # install.packages("flextable") 
library(officer) # install.packages("officer")
# Bioconductor libraries
library(Biobase) # BiocManager::install("Biobase")
library(clusterProfiler) # BiocManager::install("clusterProfiler")
# DOSE v 3.30.1 required for reproducibility
install.packages("natmed/data/DOSE_3.30.1.tar.gz", repos = NULL, type = "source")
library(DOSE)
library(org.Hs.eg.db) # BiocManager::install("org.Hs.eg.db")
library(pRoloc) # BiocManager::install("pRoloc")
# Custom operators and functions
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
# load DE data
load("natmed/results/deg.RData")


# WRANGLE GENE DATA FOR GSEA ####
data <- deg %>%
    mutate(
        genes_gsea = map(
            table,
            function(table) {
                table %>%
                    dplyr::filter(p < 0.05) %>%
                    arrange(t %>% dplyr::desc()) %>%
                    dplyr::mutate(
                        ENTREZID = clusterProfiler::bitr(Symb,
                            fromType = "SYMBOL",
                            toType = c("ENTREZID"),
                            OrgDb = org.Hs.eg.db,
                            drop = FALSE
                        ) %>% pull(ENTREZID)
                    ) %>%
                    drop_na(ENTREZID) %>%
                    dplyr::select(t, ENTREZID) %>%
                    pull(t, ENTREZID)
            }
        )
    )


# GENESET ENRICHMENT ANALYSES (GSEA) ####
set.seed(42)
gsea_dose <- data %>%
    mutate(
        gsea = map(
            genes_gsea,
            function(genes_gsea) {
                DOSE::gseDO(
                    gene = genes_gsea,
                    minGSSize = 10, maxGSSize = 200,
                    pvalueCutoff = 0.001, pAdjustMethod = "fdr", seed = TRUE
                ) %>% clusterProfiler::setReadable(OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
            }
        )
    )


# FORMAT TABLES FOR GENESET ENRICHMENT ANALYSES (GSEA) ####
gsea_dose_formatted <- gsea_dose %>%
    mutate(
        gsea_tables = map(
            gsea,
            function(gsea) {
                gsea %>%
                    as_tibble() %>%
                    arrange(pvalue) %>%
                    mutate(sign = case_when(NES < 0 ~ "Down Regulated", NES > 0 ~ "Up Regulated")) 
            }
        )
    )



# FORMAT TABLES FOR GENESET ENRICHMENT ANALYSES (GSEA) ####
gsea_dose_tables <- gsea_dose_formatted %>%
    mutate(
        gsea_flextables = map(
            gsea_tables,
            function(gsea_tables) {
                gsea_tables %>%
                    slice_min(pvalue, n = 20) %>%
                    mutate(
                        NES = NES %>% round(2),
                        pvalue = case_when(
                            pvalue < 0.0001 ~ pvalue %>% formatC(digits = 0, format = "e"),
                            TRUE ~ pvalue %>% formatC(digits = 4, format = "f")
                        ),
                        FDR = case_when(
                            p.adjust < 0.0001 ~ p.adjust %>% formatC(digits = 0, format = "e"),
                            TRUE ~ p.adjust %>% formatC(digits = 4, format = "f")
                        ),
                    ) %>%
                    dplyr::select(Description, setSize, NES, pvalue, FDR,  dplyr::any_of("core_enrichment")) %>%
                    flextable::flextable() %>%
                    flextable::border_remove() %>%
                    flextable::border(border = officer::fp_border(), part = "all") %>%
                    flextable::align(align = "center", part = "all") %>%
                    flextable::align(j = "core_enrichment", align = "left", part = "body") %>%
                    flextable::bg(bg = "white", part = "all") %>%
                    flextable::autofit()
            }
        )
    )
gsea_dose_tables$gsea_tables


# PREPARE THE RESULTS FOR EXPORT ####
gsea_dose <- gsea_dose_tables %>%
    mutate(db = "DOSE", .before = 1)
names(gsea_dose$gsea_flextables) <- gsea_dose$design


# SAVE THE GSEA RESULTS ####
saveDir <- "natmed/results/"
save(gsea_dose, file = paste(saveDir, "gsea_dose.RData", sep = ""))