# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(flextable) # install.packages("flextable") #for table outputs
library(officer) # install.packages("officer")
library(openxlsx) # install.packages("openxlsx")
library(readxl) # install.packages("readxl")
# Bioconductor libraries
library(Biobase) # BiocManager::install("Biobase")
library(limma) # BiocManager::install("limma")
library(biobroom) # BiocManager::install("biobroom")
library(genefilter) # BiocManager::install("genefilter")
# Custom operators, functions, and datasets
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
# load reference set
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/data_expressionset_k1208.RData")
# load affymap
load("Z:/DATA/Datalocks/Other data/affymap219_21Oct2019_1306_JR.RData")
# load mean expression by probe in K1208
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/mean_expression_by_probe_1208.RData")
# load in house cell panel
atagc <- read_excel("Z:/MISC/Patrick Gauthier/R/affymap219-CELL-PANEL/backup/UPDATED 2017 ANNOTATIONS - MASTERFILE - U133 HUMAN CELL PANEL - ALL PROBESETS (nonIQR) pfhptg.xlsx")
# load DEG at baseline ####
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/DEG_at_baseline_limma_1208.RData")
# load mean expression by MMDx in K1208
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/mean_expression_K1208_MMDx.RData")


# DEFINE THE SET ####
set00 <- data_expressionset_k1208[, data_expressionset_k1208$Patient %nin% c(15, 18)]



# IQR FILTER THE DATA ####
f1 <- function(x) (IQR(x) > 0.5)
ff <- filterfun(f1)
if (!exists("selected")) {
    selected <- genefilter(set00, ff)
}
set01 <- set00[selected, ]


# DEFINE GENES SIMILAR AT BASELINE ####
genes_baseline <- table_block_1 %>%
    dplyr::filter(p > 0.05) %>%
    pull(AffyID)


# WRANGLE THE MEAN EXPRESSION DATA ####
# means_1208 <- mean_exprs_1208 %>%
#     dplyr::slice_max(mean_expression, by = "Symb")

# mean_exprs_by_probe <- set01 %>%
#     exprs() %>%
#     as_tibble(rownames = "AffyID") %>%
#     right_join(affymap219 %>% dplyr::select(AffyID, Symb) %>% tibble(), ., by = "AffyID") %>%
#     mutate(mean_exprs = set01 %>%
#         exprs() %>% rowMeans(), .after = Symb)

# genes <- mean_exprs_by_probe %>%
#     group_by(Symb) %>%
#     dplyr::slice_max(mean_exprs) %>%
#     dplyr::filter(Symb != "", AffyID %in% genes_baseline) %>%
#     distinct(Symb, .keep_all = TRUE) %>%
#     pull(AffyID)

genes <- mean_exprs_1208 %>%
    dplyr::slice_max(mean_expression, by = "Symb") %>%
    dplyr::filter(Symb != "", AffyID %in% genes_baseline) %>%
    distinct(Symb, .keep_all = TRUE) %>%
    pull(AffyID)

set <- set01[featureNames(set01) %in% genes, ]


# DEFINE SEED ####
seed <- 42


# DEFINE FACTOR FOR CONTRASTS ####
Felzartamab_Followup <- set$Felzartamab_Followup %>% droplevels()
cortex <- set$Cortexprob


# DESIGN ####
design <- model.matrix(~ 0 + Felzartamab_Followup + cortex)


# CONTRAST DESIGN week24 - baseline ####
contrast_block_01 <- makeContrasts(
    "x =  (Felzartamab_FollowupWeek24_Felzartamab-Felzartamab_FollowupBaseline_Felzartamab)/2",
    levels = design
)


# CONTRAST DESIGN week52 - week24 ####
contrast_block_02 <- makeContrasts(
    "x =  (Felzartamab_FollowupWeek52_Felzartamab-Felzartamab_FollowupWeek24_Felzartamab)/2",
    levels = design
)


# CONTRAST DESIGN week52 - baseline####
contrast_block_03 <- makeContrasts(
    "x =  (Felzartamab_FollowupWeek52_Felzartamab-Felzartamab_FollowupBaseline_Felzartamab)/2",
    levels = design
)


# FIT BLOCK week24 - baseline LIMMA MODEL ####
fit_block_1 <- limma::lmFit(set, design)
cfit_block_1 <- limma::contrasts.fit(fit_block_1, contrast_block_01)
ebayes_block_1 <- limma::eBayes(cfit_block_1)
tab_block_1 <- limma::topTable(ebayes_block_1, adjust = "fdr", sort.by = "p", number = "all")
ebayes_block_1 %>% limma::topTable()


# FIT BLOCK week52 - week24 LIMMA MODEL ####
fit_block_2 <- limma::lmFit(set, design)
cfit_block_2 <- limma::contrasts.fit(fit_block_2, contrast_block_02)
ebayes_block_2 <- limma::eBayes(cfit_block_2)
tab_block_2 <- limma::topTable(ebayes_block_2, adjust = "fdr", sort.by = "p", number = "all")
ebayes_block_2 %>% limma::topTable()


# FIT BLOCK week52 - baseline LIMMA MODEL ####
fit_block_3 <- limma::lmFit(set, design)
cfit_block_3 <- limma::contrasts.fit(fit_block_3, contrast_block_03)
ebayes_block_3 <- limma::eBayes(cfit_block_3)
tab_block_3 <- limma::topTable(ebayes_block_3, adjust = "fdr", sort.by = "p", number = "all")
ebayes_block_3 %>% limma::topTable()


# CALCULATE MEAN GENE EXPRESSION FOR EACH PROBE BETWEEN GROUPINGS ####
means_baseline_week24 <- fit_block_1 %>%
    avearrays() %>%
    data.frame() %>%
    rownames_to_column("AffyID") %>%
    tibble() %>%
    mutate_if(is.numeric, ~ 2^. %>% round(0)) %>%
    rename_at(vars(contains("Felz")), ~ str_remove(., "Felzartamab_Followup")) %>%
    dplyr::select(-contains("Week12"), -contains("Placebo"), -any_of(c("cortex")))
    # dplyr::select(-contains("Week52"), -contains("Week12"), -contains("Placebo"), -any_of(c("cortex")))


# means_week24_week52 <- fit_block_2 %>%
#     avearrays() %>%
#     data.frame() %>%
#     rownames_to_column("AffyID") %>%
#     tibble() %>%
#     mutate_if(is.numeric, ~ 2^. %>% round(0)) %>%
#     rename_at(vars(contains("Felz")), ~ str_remove(., "Felzartamab_Followup")) %>%
#     dplyr::select(-contains("Week12"), -contains("Placebo"), -any_of(c("cortex")))
#     # dplyr::select(-contains("Baseline"), -contains("Week12"), -contains("Placebo"), -any_of(c("cortex")))


# means_week52_baseline <- fit_block_3 %>%
#     avearrays() %>%
#     data.frame() %>%
#     rownames_to_column("AffyID") %>%
#     tibble() %>%
#     mutate_if(is.numeric, ~ 2^. %>% round(0)) %>%
#     rename_at(vars(contains("Felz")), ~ str_remove(., "Felzartamab_Followup")) %>%
#     dplyr::select(-contains("Week12"), -contains("Placebo"), -any_of(c("cortex")))
#     # dplyr::select(-contains("Baseline"), -contains("Week12"), -contains("Placebo"), -any_of(c("cortex")))


means_felzartamab <- fit_block_1 %>%
    avearrays() %>%
    data.frame() %>%
    rownames_to_column("AffyID") %>%
    tibble() %>%
    mutate_if(is.numeric, ~ 2^. %>% round(0)) %>%
    rename_at(vars(contains("Felz")), ~ str_remove(., "Felzartamab_Followup")) %>%
    dplyr::select(-contains("Week12"), -contains("Placebo"), -any_of(c("cortex")))




# FORMAT TOPTABLES ####
table_block_1 <- tab_block_1 %>%
    as_tibble(rownames = "AffyID") %>%
    right_join(affymap219 %>% dplyr::select(AffyID, Symb, Gene, PBT), ., by = "AffyID") %>%
    arrange(P.Value) %>%
    mutate_at(c("P.Value", "adj.P.Val"), as.numeric) %>%
    tibble() %>%
    left_join(means_felzartamab, by = "AffyID") %>%
    dplyr::select(
        AffyID, Symb, Gene, PBT,
        all_of(colnames(means_felzartamab)[-1]),
        logFC, P.Value, adj.P.Val,
    ) %>%
    left_join(means_K1208 %>% dplyr::select(-Symb, -Gene, -PBT), by = "AffyID") %>%
    dplyr::rename(
        p = P.Value,
        FDR = adj.P.Val
    )

table_block_2 <- tab_block_2 %>%
    as_tibble(rownames = "AffyID") %>%
    right_join(affymap219 %>% dplyr::select(AffyID, Symb, Gene, PBT), ., by = "AffyID") %>%
    arrange(P.Value) %>%
    mutate_at(c("P.Value", "adj.P.Val"), as.numeric) %>%
    tibble() %>%
    left_join(., means_felzartamab, by = "AffyID") %>%
    dplyr::select(
        AffyID, Symb, Gene, PBT,
        all_of(colnames(means_felzartamab)[-1]),
        logFC, P.Value, adj.P.Val,
    ) %>%
    left_join(means_K1208 %>% dplyr::select(-Symb, -Gene, -PBT), by = "AffyID") %>%
    dplyr::rename(
        p = P.Value,
        FDR = adj.P.Val
    )

table_block_3 <- tab_block_3 %>%
    as_tibble(rownames = "AffyID") %>%
    right_join(affymap219 %>% dplyr::select(AffyID, Symb, Gene, PBT), ., by = "AffyID") %>%
    arrange(P.Value) %>%
    mutate_at(c("P.Value", "adj.P.Val"), as.numeric) %>%
    tibble() %>%
    left_join(., means_felzartamab, by = "AffyID") %>%
    dplyr::select(
        AffyID, Symb, Gene, PBT,
        all_of(colnames(means_felzartamab)[-1]),
        logFC, P.Value, adj.P.Val,
    ) %>%
    left_join(means_K1208 %>% dplyr::select(-Symb, -Gene, -PBT), by = "AffyID") %>%
    dplyr::rename(
        p = P.Value,
        FDR = adj.P.Val
    )

limma_tables <- tibble(
    design = c(
        "Baseline_vs_Week24",
        "Week24_vs_Week52",
        "Baseline_vs_Week52"
    ),
    toptable = list(
        table_block_1 %>%
            dplyr::rename(
                "\u394 felz logFC" = logFC,
                "\u394 felz p" = p,
                "\u394 felz FDR" = FDR
            ),
        table_block_2 %>%
            dplyr::rename(
                "\u394 felz logFC" = logFC,
                "\u394 felz p" = p,
                "\u394 felz FDR" = FDR
            ),
        table_block_3 %>%
            dplyr::rename(
                "\u394 felz logFC" = logFC,
                "\u394 felz p" = p,
                "\u394 felz FDR" = FDR
            )
    ),
    table = list(
        table_block_1,
        table_block_2,
        table_block_3
    )
)
limma_tables$table[[3]]
# tab_block_1 %>% as_tibble(rownames = "AffyID")


# EXPORT THE DATA AS .RData FILE ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/"
names(limma_tables$table) <- limma_tables$design
save(limma_tables, file = paste(saveDir, "ONLY_FELZ_IQR_filtered_probes_unique_genes_baseline_corrected_cortex_corrected_limma_1208.RData", sep = ""))


# EXPORT THE DATA AS AN EXCEL SHEET ####
names(limma_tables$toptable) <- limma_tables$design
saveDir1 <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/output/"
openxlsx::write.xlsx(limma_tables$toptable,
    asTable = TRUE,
    file = paste(saveDir1, "ONLY_FELZ_IQR_filtered_probes_unique_genes_baseline_corrected_cortex_corrected_limma_1208_3July24",
        # Sys.Date(),
        # format(Sys.time(), "_%I%M%p"),
        ".xlsx",
        sep = ""
    )
)


limma_tables %>%
    dplyr::filter(design == "Baseline_vs_Week24") %>%
    pull(table) %>%
    pluck(1) %>%
    dplyr::filter(p < 0.05) %>%
    mutate(direction = ifelse(logFC < 0, "down", "up")) %>%
    nest(.by = direction)

limma_tables %>%
    dplyr::filter(design == "Week24_vs_Week52") %>%
    pull(table) %>%
    pluck(1) %>%
    dplyr::filter(p < 0.05) %>%
    mutate(direction = ifelse(logFC < 0, "down", "up")) %>%
    nest(.by = direction)

limma_tables %>%
    dplyr::filter(design == "Baseline_vs_Week52") %>%
    pull(table) %>%
    pluck(1) %>%
    dplyr::filter(p < 0.05) %>%
    mutate(direction = ifelse(logFC < 0, "down", "up")) %>%
    nest(.by = direction)
