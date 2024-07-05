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
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/IQR_filtered_probes_unique_genes_baseline_corrected_cortex_corrected_limma_1208.RData")
# load gene lists
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/Hinze_injury_markers.RData")



# WRANGLE THE INJURY MARKER DATA ####
# genes_injury_markers %>%
#     unnest(data) %>%
#     nest(.by = c(celltypename, cluster)) %>%
#     print(n = "all")

# genes_injury_markers$data[[1]] %>%
#     print(n=100)

injury_markers <- genes_injury_markers %>%
    unnest(data) %>%
    dplyr::filter(
        cluster %>% str_detect("New"),
        # log2FC > 0.5
    ) %>%
    drop_na(AffyID) %>%
    dplyr::select(celltypename, cluster, AffyID, Symb) %>%
    dplyr::mutate(
        entrez_gene = bitr(Symb,
            fromType = "SYMBOL",
            toType = c("ENTREZID"),
            OrgDb = org.Hs.eg.db,
            drop = FALSE
        ) %>% pull(ENTREZID), .by = "cluster"
    ) %>%
    dplyr::rename(gs_name = cluster) %>%
    drop_na(entrez_gene)


# DEFINE INJURY PATHWAYS ####
map_aki <- injury_markers %>% dplyr::select(gs_name, entrez_gene)
map_aki_joined <- injury_markers %>%
    dplyr::select(gs_name, entrez_gene) %>%
    mutate(gs_name = "AKI")



# WRANGLE GENE DATA FOR GSEA ####
data <- limma_tables %>%
    mutate(
        genes_gsea = map(
            toptable,
            function(toptable) {
                toptable %>%
                    dplyr::select(AffyID, t, `<U+0394><U+0394> p`) %>%
                    dplyr::filter(`<U+0394><U+0394> p` < 0.05) %>%
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
gsea_aki <- data %>%
    mutate(
        gsea = map(
            genes_gsea,
            function(genes_gsea) {
                clusterProfiler::GSEA(
                    gene = genes_gsea, TERM2GENE = map_aki_joined,
                    minGSSize = 5, maxGSSize = Inf,
                    pvalueCutoff = 0.001, pAdjustMethod = "fdr", seed = TRUE
                ) %>% clusterProfiler::setReadable(OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
            }
        )
    )
gsea_aki$gsea[[3]]



# FORMAT TABLES FOR GENESET ENRICHMENT ANALYSES (GSEA) ####
gsea_aki_formatted <- gsea_aki %>%
    mutate(
        gsea_tables = map(
            gsea,
            function(gsea) {
                gsea %>%
                    as_tibble() %>%
                    arrange(pvalue) %>%
                    mutate(
                        sign = case_when(NES < 0 ~ "Down Regulated", NES > 0 ~ "Up Regulated"),
                        Description = "Cellular response to acute kidney injury"
                    )
            }
        )
    )
gsea_aki_formatted$gsea_tables


# FORMAT TABLES FOR GENESET ENRICHMENT ANALYSES (GSEA) ####
gsea_aki_tables <- gsea_aki_formatted %>%
    mutate(
        gsea_flextables = map(
            gsea_tables,
            function(gsea_tables) {
                gsea_tables %>%
                    # slice_min(pvalue, n = 20) %>%
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


# gsea_tables$gsea_tables[[1]]
gsea_aki_tables$gsea_flextables[[3]]


# PREPARE THE RESULTS FOR EXPORT ####
felzartamab_gsea_aki_k1208 <- gsea_aki_tables %>%
    mutate(db = "Hinze", .before = 1)
names(felzartamab_gsea_aki_k1208$gsea_flextables) <- felzartamab_gsea_aki_k1208$design



# SAVE THE GSEA RESULTS ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/"
save(felzartamab_gsea_aki_k1208, file = paste(saveDir, "felzartamab_gsea_aki_baseline_corrected_cortex_corrected_k1208.RData", sep = ""))



# EXPORT THE GSEA RESULTS TO EXCEL FILE ####
saveDir1 <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/output/"
openxlsx::write.xlsx(felzartamab_gsea_aki_k1208$gsea_tables,
    asTable = TRUE,
    file = paste(saveDir1, "aki_pathways_IQR_filtered_probes_unique_genes_baseline_corrected_cortex_corrected_limma_1208_13June24",
        # Sys.Date(),
        # format(Sys.time(), "_%I%M%p"),
        ".xlsx",
        sep = ""
    )
)
