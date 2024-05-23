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
# Custom operators, functions, and datasets
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
# load reference set
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/data_expressionset_k1208.RData")
# load affymap
load("Z:/DATA/Datalocks/Other data/affymap219_21Oct2019_1306_JR.RData")



# DEFINE SEED ####
seed <- 42


# DEFINE THE SET ####
set <- data_expressionset_k1208
set %>%
    pData() %>%
    colnames()


# DEFINE FACTOR FOR CONTRASTS ####
Felzartamab_Followup <- set$Felzartamab_Followup %>% droplevels()
design <- model.matrix(~ 0 + Felzartamab_Followup)



# NOISE DESIGN1 - Baseline Placebo vs Baseline Felzartamab - week24 - baseline ####
contrast_noise_01 <- makeContrasts(
    "x =  (Felzartamab_FollowupBaseline_Felzartamab) - (Felzartamab_FollowupBaseline_Placebo)",
    levels = design
)


# NOISE DESIGN2 - (Baseline Placebo + Week24 Placebo) vs Week24 Felzartamab - week24 - baseline ####
contrast_noise_02 <- makeContrasts(
    "x =  (Felzartamab_FollowupWeek24_Felzartamab) - (Felzartamab_FollowupWeek24_Placebo + Felzartamab_FollowupBaseline_Placebo)/2",
    levels = design
)


# FIT NOISE1 week24 - baseline LIMMA MODEL ####
fit_noise_1 <- limma::lmFit(set, design)
cfit_noise_1 <- limma::contrasts.fit(fit_noise_1, contrast_noise_01)
ebayes_noise_1 <- limma::eBayes(cfit_noise_1)
tab_noise_1 <- limma::topTable(ebayes_noise_1, adjust = "fdr", sort.by = "p", number = "all")
ebayes_noise_1 %>% limma::topTable()


# FIT NOISE2 week24 - baseline LIMMA MODEL ####
fit_noise_2 <- limma::lmFit(set, design)
cfit_noise_2 <- limma::contrasts.fit(fit_noise_2, contrast_noise_02)
ebayes_noise_2 <- limma::eBayes(cfit_noise_2)
tab_noise_2 <- limma::topTable(ebayes_noise_2, adjust = "fdr", sort.by = "p", number = "all")
ebayes_noise_2 %>% limma::topTable()



# CALCULATE MEAN GENE EXPRESSION FOR EACH PROBE BETWEEN GROUPINGS ####
means_baseline_week24 <- fit_noise_1 %>%
    avearrays() %>%
    data.frame() %>%
    rownames_to_column("AffyID") %>%
    tibble() %>%
    mutate_if(is.numeric, ~ 2^. %>% round(0)) %>%
    rename_at(vars(contains("Felz")), ~ str_remove(., "Felzartamab_Followup")) %>%
    dplyr::select(-contains("Week52"), -contains("Week12"))

# means_week24_week52 <- fit_block_2 %>%
#     avearrays() %>%
#     data.frame() %>%
#     rownames_to_column("AffyID") %>%
#     tibble() %>%
#     mutate_if(is.numeric, ~ 2^. %>% round(0)) %>%
#     rename_at(vars(contains("Felz")), ~ str_remove(., "Felzartamab_Followup")) %>%
#     dplyr::select(-contains("Baseline"), -contains("Week12"))

# means_week52_baseline <- fit_block_3 %>%
#     avearrays() %>%
#     data.frame() %>%
#     rownames_to_column("AffyID") %>%
#     tibble() %>%
#     mutate_if(is.numeric, ~ 2^. %>% round(0)) %>%
#     rename_at(vars(contains("Felz")), ~ str_remove(., "Felzartamab_Followup")) %>%
#     dplyr::select(-contains("Week24"), -contains("Week12"))


# FORMAT TOPTABLES ####
table_noise_1 <- tab_noise_1 %>%
    as_tibble(rownames = "AffyID") %>%
    right_join(affymap219 %>% dplyr::select(AffyID, Symb, Gene, PBT), ., by = "AffyID") %>%
    arrange(logFC) %>%
    mutate_at(c("P.Value", "adj.P.Val"), as.numeric) %>%
    tibble() %>%
    left_join(., means_baseline_week24, by = "AffyID") %>%
    dplyr::select(
        AffyID, Symb, PBT,
        all_of(colnames(means_baseline_week24)[-1]),
        logFC, P.Value
    ) %>%
    dplyr::rename(
        logFC_noise_1 = logFC,
        p_noise_1 = P.Value
    )

table_noise_2 <- tab_noise_2 %>%
    as_tibble(rownames = "AffyID") %>%
    right_join(affymap219 %>% dplyr::select(AffyID, Symb, Gene, PBT), ., by = "AffyID") %>%
    arrange(P.Value) %>%
    mutate_at(c("P.Value", "adj.P.Val"), as.numeric) %>%
    tibble() %>%
    left_join(., means_baseline_week24, by = "AffyID") %>%
    dplyr::select(
        AffyID, Symb, PBT,
        all_of(colnames(means_baseline_week24)[-1]),
        logFC, P.Value
    ) %>%
    dplyr::rename(
        logFC_noise_2 = logFC,
        p_noise_2 = P.Value
    )


table_noise <- left_join(table_noise_1, table_noise_2)


# ISOLATE NOISE PROBES ####
# these gense are significantly different between Placebo and Felzartamab Baseline
probes_noise1 <- table_noise %>%
    dplyr::filter(
        logFC_noise_2 %>% abs() < 0.75, logFC_noise_1 %>% abs() > 1
    ) %>%
    pull(AffyID)

probes_noise <- table_noise %>%
    dplyr::filter(
        p_noise_1 < 0.05, p_noise_2 > 0.05
    ) %>%
    pull(AffyID)

table_noise %>%
    dplyr::filter(AffyID %in% probes_noise1)
table_noise %>%
    dplyr::filter(AffyID %in% probes_noise)


# EXPORT THE DATA AS .RData FILE ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/"
save(probes_noise, file = paste(saveDir, "probes_noise_1208.RData", sep = ""))

