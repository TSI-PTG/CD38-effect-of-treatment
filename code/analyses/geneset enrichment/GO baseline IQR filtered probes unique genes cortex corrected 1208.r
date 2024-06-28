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
library(pRoloc) # BiocManager::install("pRoloc")
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
# load affymap
load("Z:/DATA/Datalocks/Other data/affymap219_21Oct2019_1306_JR.RData")
# load DE data
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/DEG_at_baseline_limma_1208.RData")


# WRANGLE GENE DATA FOR GSEA ####
data <- table_block_1 %>%
    mutate(design = "baseline") %>%
    nest(.by = design, toptable = everything())  %>% 
    mutate(
        genes_gsea = map(
            toptable,
            function(toptable) {
                toptable %>%
                    dplyr::filter(p < 0.05) %>%
                    dplyr::arrange(t  %>% dplyr::desc())  %>% 
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

data$genes_gsea



# GENESET ENRICHMENT ANALYSES (GSEA) ####
set.seed(42)
gsea_go <- data %>%
    mutate(
        gsea_go = map(
            genes_gsea,
            function(genes_gsea) {
                clusterProfiler::gseGO(
                    gene = genes_gsea, ont = "BP", OrgDb = org.Hs.eg.db,
                    minGSSize = 10, maxGSSize = 200,
                    pvalueCutoff = 0.05, pAdjustMethod = "fdr", seed = TRUE
                ) %>% clusterProfiler::setReadable(OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
            }
        )
    )


# SIMPLIFY GO TERMS ####
set.seed(42)
gsea_go_simplified <- gsea_go %>%
    mutate(gsea_go_simplified = map(gsea_go, clusterProfiler::simplify))



# FORMAT TABLES FOR GENESET ENRICHMENT ANALYSES (GSEA) ####
gsea_go_simplified_formatted <- gsea_go_simplified %>%
    mutate(
        gsea_go_tables = map(
            gsea_go_simplified,
            function(gsea_go_simplified) {
                gsea_go_simplified %>%
                    as_tibble() %>%
                    arrange(pvalue) %>%
                    mutate(
                        sign = case_when(NES < 0 ~ "Down Regulated", NES > 0 ~ "Up Regulated"),
                        Description = factor(Description, levels = Description, ordered = TRUE)
                    ) %>%
                    mutate(
                        ID_parent = GO.db::GOBPPARENTS[[ID]][[1]], .by = ID,
                        .after = ID,
                    ) %>%
                    mutate(
                        Description_parent = pRoloc::goIdToTerm(ID_parent),
                        .after = Description,
                    )
            }
        )
    )



# FORMAT TABLES FOR GENESET ENRICHMENT ANALYSES (GSEA) ####
gsea_go_tables <- gsea_go_simplified_formatted %>%
    mutate(
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
                    flextable::autofit()
            }
        )
    )


# gsea_go_tables$gsea_go_tables[[1]]
gsea_go_tables$gsea_go_flextables[[1]]


# PREPARE THE RESULTS FOR EXPORT ####
felzartamab_gsea_k1208 <- gsea_go_tables
names(felzartamab_gsea_k1208$gsea_go_flextables) <- felzartamab_gsea_k1208$design



# SAVE THE GSEA RESULTS ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/"
save(felzartamab_gsea_k1208, file = paste(saveDir, "baseline_gsea_cortex_corrected_k1208.RData", sep = ""))



# EXPORT THE GSEA RESULTS TO EXCEL FILE ####
saveDir1 <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/output/"
openxlsx::write.xlsx(felzartamab_gsea_k1208$gsea_go_tables,
    asTable = TRUE,
    file = paste(saveDir1, "Baseline_Pathways_IQR_filtered_probes_unique_genes_cortex_corrected_limma_1208_26June24",
        # Sys.Date(),
        # format(Sys.time(), "_%I%M%p"),
        ".xlsx",
        sep = ""
    )
)
