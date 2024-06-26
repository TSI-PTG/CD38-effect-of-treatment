# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(flextable) # install.packages("flextable") #for table outputs
library(officer) # install.packages("officer")
library(openxlsx) # install.packages("openxlsx")
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
# load mean expression in K1208
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


# KEEP UNIQUE GENES (keep probe with highest mean expression) ####
mean_exprs_by_probe <- set01 %>%
    exprs() %>%
    as_tibble(rownames = "AffyID") %>%
    right_join(affymap219 %>% dplyr::select(AffyID, Symb) %>% tibble(), ., by = "AffyID") %>%
    mutate(mean_exprs = set01 %>%
        exprs() %>% rowMeans(), .after = Symb)

genes <- mean_exprs_by_probe %>%
    group_by(Symb) %>%
    dplyr::slice_max(mean_exprs) %>%
    dplyr::filter(Symb != "") %>%
    distinct(Symb, .keep_all = TRUE) %>%
    pull(AffyID)

set <- set01[featureNames(set01) %in% genes, ]


# DEFINE SEED ####
seed <- 42


# DEFINE FACTOR FOR CONTRASTS ####
Felzartamab_Followup <- set$Felzartamab_Followup %>% droplevels()


# BLOCK DESIGN week24 - baseline ####
design_block01 <- model.matrix(~ 0 + Felzartamab_Followup)
contrast_block_01 <- makeContrasts(
    "x =  (Felzartamab_FollowupBaseline_Felzartamab) -(Felzartamab_FollowupBaseline_Placebo)",
    levels = design_block01
)


# FIT BLOCK week24 - baseline LIMMA MODEL ####
fit_block_1 <- limma::lmFit(set, design_block01)
cfit_block_1 <- limma::contrasts.fit(fit_block_1, contrast_block_01)
ebayes_block_1 <- limma::eBayes(cfit_block_1)
tab_block_1 <- limma::topTable(ebayes_block_1, adjust = "fdr", sort.by = "p", number = "all")
ebayes_block_1 %>% limma::topTable()



# CALCULATE MEAN GENE EXPRESSION FOR EACH PROBE BETWEEN GROUPINGS ####
means_baseline <- fit_block_1 %>%
    avearrays() %>%
    data.frame() %>%
    rownames_to_column("AffyID") %>%
    tibble() %>%
    mutate_if(is.numeric, ~ 2^. %>% round(0)) %>%
    rename_at(vars(contains("Felz")), ~ str_remove(., "Felzartamab_Followup")) %>%
    dplyr::select(AffyID, contains("Baseline"))


# FORMAT TOPTABLES ####
table_block_1 <- tab_block_1 %>%
    as_tibble(rownames = "AffyID") %>%
    right_join(affymap219 %>% dplyr::select(AffyID, Symb, Gene, PBT), ., by = "AffyID") %>%
    arrange(P.Value) %>%
    mutate_at(c("P.Value", "adj.P.Val"), as.numeric) %>%
    tibble() %>%
    left_join(means_baseline, by = "AffyID") %>%
    dplyr::select(
        AffyID, Symb, Gene, PBT,
        Baseline_Placebo, Baseline_Felzartamab,
        t, logFC, P.Value, adj.P.Val,
    ) %>%
    mutate(
        FC = 2^logFC,
        .after = logFC
    ) %>%
    dplyr::rename(p = P.Value, FDR = adj.P.Val) 


table_block_1 %>%
    dplyr::filter(abs(logFC) > 0.5)


# table_block_1 %>%
#     dplyr::filter(Symb == "HSP90B1")


# EXPORT THE DATA AS .RData FILE ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/"
save(table_block_1, file = paste(saveDir, "DEG_at_baseline_limma_1208.RData", sep = ""))


# EXPORT THE DATA AS AN EXCEL SHEET ####
saveDir1 <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/output/"
openxlsx::write.xlsx(table_block_1,
    asTable = TRUE,
    file = paste(saveDir1, "DEG_at_baseline_limma_1208_13Jun24",
        # Sys.Date(),
        # format(Sys.time(), "_%I%M%p"),
        ".xlsx",
        sep = ""
    )
)
