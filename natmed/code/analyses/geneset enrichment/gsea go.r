# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(flextable) # install.packages("flextable") 
library(officer) # install.packages("officer")
# Bioconductor libraries
library(Biobase) # BiocManager::install("Biobase")
library(clusterProfiler) # BiocManager::install("clusterProfiler")
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
        gsea = map(
            genes_gsea,
            function(genes_gsea) {
                clusterProfiler::gseGO(
                    gene = genes_gsea, ont = "BP", OrgDb = org.Hs.eg.db,
                    minGSSize = 10, maxGSSize = 200,
                    pvalueCutoff = 0.001, pAdjustMethod = "fdr", seed = TRUE
                ) %>% clusterProfiler::setReadable(OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
            }
        )
    )


# SIMPLIFY GO TERMS ####
set.seed(42)
gsea_go_simplified <- gsea_go %>%
    mutate(gsea_simplified = map(gsea, clusterProfiler::simplify))



# FORMAT TABLES FOR GENESET ENRICHMENT ANALYSES (GSEA) ####
gsea_go_simplified_formatted <- gsea_go_simplified %>%
    mutate(
        gsea_tables = map(
            gsea_simplified,
            function(gsea_simplified) {
                gsea_simplified %>%
                    as_tibble() %>%
                    arrange(pvalue) %>%
                    mutate(
                        sign = case_when(NES < 0 ~ "Down Regulated", NES > 0 ~ "Up Regulated"),
                        Description = factor(Description, levels = Description, ordered = TRUE)
                    ) 
                    # %>%
                    # mutate(
                    #     ID_parent = GO.db::GOBPPARENTS[[ID]][[1]], .by = ID,
                    #     .after = ID,
                    # ) %>%
                    # mutate(
                    #     Description_parent = pRoloc::goIdToTerm(ID_parent),
                    #     .after = Description,
                    # )
            }
        )
    )



# FORMAT TABLES FOR GENESET ENRICHMENT ANALYSES (GSEA) ####
gsea_go_tables <- gsea_go_simplified_formatted %>%
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
                    dplyr::select(Description, setSize, NES, pvalue, FDR, dplyr::any_of("core_enrichment")) %>%
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
gsea_go_tables$gsea_flextables[[1]]


# PREPARE THE RESULTS FOR EXPORT ####
gsea_go <- gsea_go_tables %>%
    mutate(db = "GO", .before = 1)
names(gsea_go$gsea_flextables) <- gsea_go$design


# SAVE THE GSEA RESULTS ####
saveDir <- "natmed/results/"
save(gsea_go, file = paste(saveDir, "gsea_go.RData", sep = ""))

