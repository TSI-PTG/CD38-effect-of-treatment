# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(flextable) # install.packages("flextable") #for simple table outputs
library(officer) # install.packages("officer")
library(clusterProfiler) # pak::pak("YuLab-SMU/clusterProfiler")
library(circlize) # install.packages("circlize") pak::pak("jokergoo/circlize")
# Bioconductor libraries
library(Biobase) # BiocManager::install("Biobase")
library(biobroom) # BiocManager::install("biobroom")
library(org.Hs.eg.db) # BiocManager::install("org.Hs.eg.db")
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
# load affymap
load("Z:/DATA/Datalocks/Other data/affymap219_21Oct2019_1306_JR.RData")
# load DE data
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/IQR_filtered_probes_unique_genes_limma_1208.RData")


# WRANGLE GENE DATA FOR GSEA ####
data <- limma_tables %>%
    mutate(
        genes_gsea = map(
            toptable,
            function(toptable) {
                toptable %>%
                    dplyr::filter(P.Value < 0.05) %>%
                    right_join(affymap219 %>% dplyr::select(AffyID, Symb) %>% tibble(), ., by = "AffyID") %>%
                    arrange(t %>% dplyr::desc()) %>%
                    dplyr::mutate(
                        ENTREZID = bitr(Symb,
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
gsea_go <- data %>%
    mutate(
        gsea_go = map(
            genes_gsea,
            function(genes_gsea) {
                clusterProfiler::gseGO(
                    gene = genes_gsea, ont = "ALL", OrgDb = org.Hs.eg.db,
                    minGSSize = 10, maxGSSize = 100,
                    pvalueCutoff = 0.05, pAdjustMethod = "fdr", seed = TRUE
                ) %>% clusterProfiler::setReadable(OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
            }
        )
    )



# FORMAT TABLES FOR GENESET ENRICHMENT ANALYSES (GSEA) ####
gsea_go_tables <- gsea_go %>%
    mutate(
        gsea_go_tables = map(
            gsea_go,
            function(gsea_go) {
                gsea_go %>%
                    data.frame() %>%
                    tibble() %>%
                    arrange(pvalue) %>%
                    mutate(sign = case_when(NES < 0 ~ "Down Regulated", NES > 0 ~ "Up Regulated")) %>%
                    mutate(Description = factor(Description, levels = Description, ordered = TRUE))
            }
        ),
        gsea_go_flextables = map(
            gsea_go_tables,
            function(gsea_go_tables) {
                gsea_go_tables %>%
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
                    dplyr::select(Description, setSize, NES, pvalue, FDR, core_enrichment) %>%
                    flextable::flextable() %>%
                    flextable::border_remove() %>%
                    flextable::border(border = fp_border(), part = "all") %>%
                    flextable::align(align = "center", part = "all") %>%
                    flextable::align(j = "core_enrichment", align = "left", part = "body") %>%
                    flextable::bg(bg = "white", part = "all") %>%
                    # flextable::bg(i = ~ NES < 0, bg = "grey80", part = "body") %>%
                    flextable::autofit()
            }
        )
    )

gsea_go_tables$gsea_go_flextables[[3]]


# MAKE SOME CLEAN FLEXTABLES ####
GSEA_GO_flextable <- df_GSEA_GO %>%
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
    dplyr::select(Description, setSize, NES, pvalue, FDR, core_enrichment) %>%
    flextable::flextable() %>%
    flextable::border_remove() %>%
    flextable::border(border = fp_border(), part = "all") %>%
    flextable::align(align = "center", part = "all") %>%
    flextable::align(j = "core_enrichment", align = "left", part = "body") %>%
    flextable::bg(bg = "white", part = "all") %>%
    # flextable::bg(i = ~ NES < 0, bg = "grey80", part = "body") %>%
    flextable::autofit()
