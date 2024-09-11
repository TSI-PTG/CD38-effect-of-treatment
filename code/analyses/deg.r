# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(flextable) # install.packages("flextable")
library(officer) # install.packages("officer")
# Bioconductor libraries
library(Biobase) # BiocManager::install("Biobase")
library(limma) # BiocManager::install("limma")
library(genefilter) # BiocManager::install("genefilter")
# Custom operators and functions
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
# load reference set
load("data/cd38.RData")
# load affymap
load("data/affymap219.RData")
# load DEG at baseline
load("results/baseline_deg.RData")


# DEFINE SEED ####
set.seed(42)


# DEFINE THE SET ####
set00 <- set[, set$Patient %nin% c(15, 18)]


# IQR FILTER THE DATA ####
f1 <- function(x) (IQR(x) > 0.5)
ff <- filterfun(f1)
if (!exists("selected")) {
    selected <- genefilter(set00, ff)
}
set01 <- set00[selected, ]


# DEFINE GENES SIMILAR AT BASELINE ####
genes_baseline <- baseline_deg %>%
    dplyr::filter(p > 0.05) %>%
    pull(AffyID)


# WRANGLE THE MEAN EXPRESSION DATA ####
mean_exprs_by_probe <- set01 %>%
    exprs() %>%
    as_tibble(rownames = "AffyID") %>%
    right_join(affymap219 %>% dplyr::select(AffyID, Symb) %>% tibble(), ., by = "AffyID") %>%
    mutate(
        mean_exprs = set01 %>%
            exprs() %>%
            rowMeans(),
        .after = Symb
    )

genes <- mean_exprs_by_probe %>%
    group_by(Symb) %>%
    dplyr::slice_max(mean_exprs) %>%
    dplyr::filter(Symb != "", AffyID %in% genes_baseline) %>%
    distinct(Symb, .keep_all = TRUE) %>%
    pull(AffyID)

set02 <- set01[featureNames(set01) %in% genes, ]


# DEFINE FACTOR FOR CONTRASTS ####
Felzartamab_Followup <- set02$Felzartamab_Followup %>% droplevels()
cortex <- set02$Cortexprob


# DESIGN ####
design <- model.matrix(~ 0 + Felzartamab_Followup + cortex)


# CONTRAST DESIGN week24 - baseline ####
contrast_interaction_01 <- makeContrasts(
    "x =  (Felzartamab_FollowupWeek24_Felzartamab-Felzartamab_FollowupBaseline_Felzartamab)/2 -(Felzartamab_FollowupWeek24_Placebo-Felzartamab_FollowupBaseline_Placebo)/2",
    levels = design
)
contrast_felzartamab_01 <- makeContrasts(
    "x =  (Felzartamab_FollowupWeek24_Felzartamab-Felzartamab_FollowupBaseline_Felzartamab)/2",
    levels = design
)
contrast_placebo_01 <- makeContrasts(
    "x =  (Felzartamab_FollowupWeek24_Placebo-Felzartamab_FollowupBaseline_Placebo)/2",
    levels = design
)


# CONTRAST DESIGN week52 - week24 ####
contrast_interaction_02 <- makeContrasts(
    "x =  (Felzartamab_FollowupWeek52_Felzartamab-Felzartamab_FollowupWeek24_Felzartamab)/2 - (Felzartamab_FollowupWeek52_Placebo-Felzartamab_FollowupWeek24_Placebo)/2",
    levels = design
)
contrast_felzartamab_02 <- makeContrasts(
    "x =  (Felzartamab_FollowupWeek52_Felzartamab-Felzartamab_FollowupWeek24_Felzartamab)/2",
    levels = design
)
contrast_placebo_02 <- makeContrasts(
    "x =  (Felzartamab_FollowupWeek52_Placebo-Felzartamab_FollowupWeek24_Placebo)/2",
    levels = design
)


# CONTRAST DESIGN week52 - baseline####
contrast_interaction_03 <- makeContrasts(
    "x =  (Felzartamab_FollowupWeek52_Felzartamab-Felzartamab_FollowupBaseline_Felzartamab)/2 - (Felzartamab_FollowupWeek52_Placebo-Felzartamab_FollowupBaseline_Placebo)/2",
    levels = design
)
contrast_felzartamab_03 <- makeContrasts(
    "x =  (Felzartamab_FollowupWeek52_Felzartamab-Felzartamab_FollowupBaseline_Felzartamab)/2",
    levels = design
)
contrast_placebo_03 <- makeContrasts(
    "x =  (Felzartamab_FollowupWeek52_Placebo-Felzartamab_FollowupBaseline_Placebo)/2",
    levels = design
)


# FIT BLOCK week24 - baseline LIMMA MODEL ####
fit_interaction_1 <- limma::lmFit(set02, design)
cfit_interaction_1 <- limma::contrasts.fit(fit_interaction_1, contrast_interaction_01)
ebayes_interaction_1 <- limma::eBayes(cfit_interaction_1)
tab_interaction_1 <- limma::topTable(ebayes_interaction_1, adjust = "fdr", sort.by = "p", number = "all") %>%
    as_tibble(rownames = "AffyID") %>%
    mutate(se = (ebayes_interaction_1$stdev.unscaled * sqrt(ebayes_interaction_1$s2.post)) %>% as.vector(), .after = logFC) %>%
    right_join(affymap219 %>% dplyr::select(AffyID, Symb, Gene, PBT), ., by = "AffyID")

fit_felzartamab_1 <- limma::lmFit(set02, design)
cfit_felzartamab_1 <- limma::contrasts.fit(fit_felzartamab_1, contrast_felzartamab_01)
ebayes_felzartamab_1 <- limma::eBayes(cfit_felzartamab_1)
tab_felzartamab_1 <- limma::topTable(ebayes_felzartamab_1, adjust = "fdr", sort.by = "p", number = "all") %>%
    as_tibble(rownames = "AffyID") %>%
    mutate(se = (ebayes_felzartamab_1$stdev.unscaled * sqrt(ebayes_felzartamab_1$s2.post)) %>% as.vector(), .after = logFC) %>%
    right_join(affymap219 %>% dplyr::select(AffyID, Symb, Gene, PBT), ., by = "AffyID")

fit_placebo_1 <- limma::lmFit(set02, design)
cfit_placebo_1 <- limma::contrasts.fit(fit_placebo_1, contrast_placebo_01)
ebayes_placebo_1 <- limma::eBayes(cfit_placebo_1)
tab_placebo_1 <- limma::topTable(ebayes_placebo_1, adjust = "fdr", sort.by = "p", number = "all") %>%
    as_tibble(rownames = "AffyID") %>%
    mutate(se = (ebayes_placebo_1$stdev.unscaled * sqrt(ebayes_placebo_1$s2.post)) %>% as.vector(), .after = logFC) %>%
    right_join(affymap219 %>% dplyr::select(AffyID, Symb, Gene, PBT), ., by = "AffyID")


# FIT BLOCK week52 - week24 LIMMA MODEL ####
fit_interaction_2 <- limma::lmFit(set02, design)
cfit_interaction_2 <- limma::contrasts.fit(fit_interaction_2, contrast_interaction_02)
ebayes_interaction_2 <- limma::eBayes(cfit_interaction_2)
tab_interaction_2 <- limma::topTable(ebayes_interaction_2, adjust = "fdr", sort.by = "p", number = "all") %>%
    as_tibble(rownames = "AffyID") %>%
    mutate(se = (ebayes_interaction_2$stdev.unscaled * sqrt(ebayes_interaction_2$s2.post)) %>% as.vector(), .after = logFC) %>%
    right_join(affymap219 %>% dplyr::select(AffyID, Symb, Gene, PBT), ., by = "AffyID")

fit_felzartamab_2 <- limma::lmFit(set02, design)
cfit_felzartamab_2 <- limma::contrasts.fit(fit_felzartamab_2, contrast_felzartamab_02)
ebayes_felzartamab_2 <- limma::eBayes(cfit_felzartamab_2)
tab_felzartamab_2 <- limma::topTable(ebayes_felzartamab_2, adjust = "fdr", sort.by = "p", number = "all") %>%
    as_tibble(rownames = "AffyID") %>%
    mutate(se = (ebayes_felzartamab_2$stdev.unscaled * sqrt(ebayes_felzartamab_2$s2.post)) %>% as.vector(), .after = logFC) %>%
    right_join(affymap219 %>% dplyr::select(AffyID, Symb, Gene, PBT), ., by = "AffyID")

fit_placebo_2 <- limma::lmFit(set02, design)
cfit_placebo_2 <- limma::contrasts.fit(fit_placebo_2, contrast_placebo_02)
ebayes_placebo_2 <- limma::eBayes(cfit_placebo_2)
tab_placebo_2 <- limma::topTable(ebayes_placebo_2, adjust = "fdr", sort.by = "p", number = "all") %>%
    as_tibble(rownames = "AffyID") %>%
    mutate(se = (ebayes_placebo_2$stdev.unscaled * sqrt(ebayes_placebo_2$s2.post)) %>% as.vector(), .after = logFC) %>%
    right_join(affymap219 %>% dplyr::select(AffyID, Symb, Gene, PBT), ., by = "AffyID")


# FIT BLOCK week52 - baseline LIMMA MODEL ####
fit_interaction_3 <- limma::lmFit(set02, design)
cfit_interaction_3 <- limma::contrasts.fit(fit_interaction_3, contrast_interaction_03)
ebayes_interaction_3 <- limma::eBayes(cfit_interaction_3)
tab_interaction_3 <- limma::topTable(ebayes_interaction_3, adjust = "fdr", sort.by = "p", number = "all") %>%
    as_tibble(rownames = "AffyID") %>%
    mutate(se = (ebayes_interaction_3$stdev.unscaled * sqrt(ebayes_interaction_3$s2.post)) %>% as.vector(), .after = logFC) %>%
    right_join(affymap219 %>% dplyr::select(AffyID, Symb, Gene, PBT), ., by = "AffyID")

fit_felzartamab_3 <- limma::lmFit(set02, design)
cfit_felzartamab_3 <- limma::contrasts.fit(fit_felzartamab_3, contrast_felzartamab_03)
ebayes_felzartamab_3 <- limma::eBayes(cfit_felzartamab_3)
tab_felzartamab_3 <- limma::topTable(ebayes_felzartamab_3, adjust = "fdr", sort.by = "p", number = "all") %>%
    as_tibble(rownames = "AffyID") %>%
    mutate(se = (ebayes_felzartamab_3$stdev.unscaled * sqrt(ebayes_felzartamab_3$s2.post)) %>% as.vector(), .after = logFC) %>%
    right_join(affymap219 %>% dplyr::select(AffyID, Symb, Gene, PBT), ., by = "AffyID")

fit_placebo_3 <- limma::lmFit(set02, design)
cfit_placebo_3 <- limma::contrasts.fit(fit_placebo_3, contrast_placebo_03)
ebayes_placebo_3 <- limma::eBayes(cfit_placebo_3)
tab_placebo_3 <- limma::topTable(ebayes_placebo_3, adjust = "fdr", sort.by = "p", number = "all") %>%
    as_tibble(rownames = "AffyID") %>%
    mutate(se = (ebayes_placebo_3$stdev.unscaled * sqrt(ebayes_placebo_3$s2.post)) %>% as.vector(), .after = logFC) %>%
    right_join(affymap219 %>% dplyr::select(AffyID, Symb, Gene, PBT), ., by = "AffyID")


# MERGE BLOCKS ####
tab_block_1 <- tab_interaction_1 %>%
    left_join(
        tab_felzartamab_1 %>%
            dplyr::select(AffyID, logFC) %>%
            rename(flogFC = logFC),
        by = "AffyID"
    ) %>%
    left_join(
        tab_placebo_1 %>%
            dplyr::select(AffyID, logFC) %>%
            rename(plogFC = logFC),
        by = "AffyID"
    ) %>%
    relocate(plogFC, flogFC, .before = logFC) %>%
    arrange(P.Value)

tab_block_2 <- tab_interaction_2 %>%
    left_join(
        tab_felzartamab_2 %>%
            dplyr::select(AffyID, logFC) %>%
            rename(flogFC = logFC),
        by = "AffyID"
    ) %>%
    left_join(
        tab_placebo_2 %>%
            dplyr::select(AffyID, logFC) %>%
            rename(plogFC = logFC),
        by = "AffyID"
    ) %>%
    relocate(plogFC, flogFC, .before = logFC) %>%
    arrange(P.Value)

tab_block_3 <- tab_interaction_3 %>%
    left_join(
        tab_felzartamab_3 %>%
            dplyr::select(AffyID, logFC) %>%
            rename(flogFC = logFC),
        by = "AffyID"
    ) %>%
    left_join(
        tab_placebo_3 %>%
            dplyr::select(AffyID, logFC) %>%
            rename(plogFC = logFC),
        by = "AffyID"
    ) %>%
    relocate(plogFC, flogFC, .before = logFC) %>%
    arrange(P.Value)


# CALCULATE MEAN GENE EXPRESSION FOR EACH PROBE BETWEEN GROUPINGS ####
means <- fit_interaction_1 %>%
    avearrays() %>%
    as_tibble(rownames = "AffyID") %>%
    mutate_if(is.numeric, ~ 2^. %>% round(0)) %>%
    rename_at(vars(contains("Felz")), ~ str_remove(., "Felzartamab_Followup")) %>%
    dplyr::select(-contains("Week12"), -any_of(c("cortex")))


# FORMAT TOPTABLES ####
table_interaction_1 <- tab_block_1 %>%
    arrange(P.Value) %>%
    mutate_at(c("P.Value", "adj.P.Val"), as.numeric) %>%
    tibble() %>%
    left_join(means, by = "AffyID") %>%
    dplyr::select(
        AffyID, Symb, Gene, PBT,
        t, plogFC, flogFC, logFC, se, P.Value, adj.P.Val,
        all_of(colnames(means))
    ) %>%
    dplyr::rename(
        p = P.Value,
        FDR = adj.P.Val
    )

table_interaction_2 <- tab_block_2 %>%
    arrange(P.Value) %>%
    mutate_at(c("P.Value", "adj.P.Val"), as.numeric) %>%
    tibble() %>%
    left_join(means, by = "AffyID") %>%
    dplyr::select(
        AffyID, Symb, Gene, PBT,
        t, plogFC, flogFC, logFC, se, P.Value, adj.P.Val,
        all_of(colnames(means))
    ) %>%
    dplyr::rename(
        p = P.Value,
        FDR = adj.P.Val
    )

table_interaction_3 <- tab_block_3 %>%
    arrange(P.Value) %>%
    mutate_at(c("P.Value", "adj.P.Val"), as.numeric) %>%
    tibble() %>%
    left_join(means, by = "AffyID") %>%
    dplyr::select(
        AffyID, Symb, Gene, PBT,
        t, plogFC, flogFC, logFC, se, P.Value, adj.P.Val,
        all_of(colnames(means))
    ) %>%
    dplyr::rename(
        p = P.Value,
        FDR = adj.P.Val
    )

deg <- tibble(
    design = c(
        "Baseline_vs_Week24",
        "Week24_vs_Week52",
        "Baseline_vs_Week52"
    ),
    table = list(
        table_interaction_1,
        table_interaction_2,
        table_interaction_3
    ),
    table_renamed = list(
        table_interaction_1 %>%
            dplyr::rename(
                "\u394 placebo logFC" = plogFC,
                "\u394 felz logFC" = flogFC,
                "\u394\u394 logFC" = logFC,
                "\u394\u394 p" = p,
                "\u394\u394 FDR" = FDR
            ),
        table_interaction_2 %>%
            dplyr::rename(
                "\u394 placebo logFC" = plogFC,
                "\u394 felz logFC" = flogFC,
                "\u394\u394 logFC" = logFC,
                "\u394\u394 p" = p,
                "\u394\u394 FDR" = FDR
            ),
        table_interaction_3 %>%
            dplyr::rename(
                "\u394 placebo logFC" = plogFC,
                "\u394 felz logFC" = flogFC,
                "\u394\u394 logFC" = logFC,
                "\u394\u394 p" = p,
                "\u394\u394 FDR" = FDR
            )
    )
)


# EXPORT THE DATA AS .RData FILE ####
saveDir <- "results/"
names(deg$table) <- deg$design
save(deg, file = paste(saveDir, "deg.RData", sep = ""))
