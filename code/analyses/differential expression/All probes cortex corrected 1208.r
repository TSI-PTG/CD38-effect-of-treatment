# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(flextable) # install.packages("flextable")  
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
affymap219 <- affymap219 %>% tibble()
# load mean expression by probe in K1208
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/mean_expression_by_probe_1208.RData")
# load in house cell panel
atagc <- read_excel("Z:/MISC/Patrick Gauthier/R/affymap219-CELL-PANEL/backup/UPDATED 2017 ANNOTATIONS - MASTERFILE - U133 HUMAN CELL PANEL - ALL PROBESETS (nonIQR) pfhptg.xlsx")
# load mean expression by MMDx in K1208
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/mean_expression_K1208_MMDx.RData")


# DEFINE SEED ####
seed <- 42
set.seed(seed)


# DEFINE THE SET ####
set01 <- data_expressionset_k1208[, data_expressionset_k1208$Patient %nin% c(15, 18)]


# WRANGLE THE CELL PANEL DATA ####
cell_panel <- atagc %>%
    dplyr::select(`Affy Probeset ID`, `Index`, `Unstim HUVEC`, `HUVEC + IFNg`, `Unstim RPTEC`, `RPTEC + IFNg`) %>%
    dplyr::rename(
        AffyID_U133 = `Affy Probeset ID`,
        Symb = Index,
        `HUVEC (unstimulated)` = `Unstim HUVEC`,
        `HUVEC (IFNg stimulated)` = `HUVEC + IFNg`,
        `RPTEC (unstimulated)` = `Unstim RPTEC`,
        `RPTEC (IFNg stimulated)` = `RPTEC + IFNg`
    ) %>%
    dplyr::slice_max(`HUVEC (unstimulated)`, by = "Symb", with_ties = FALSE)

means_K1208_cell_panel <- means_K1208 %>%
    left_join(cell_panel %>% dplyr::select(-AffyID_U133))


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
    # distinct(Symb, .keep_all = TRUE) %>%
    pull(AffyID)

set <- set01[featureNames(set01) %in% genes, ]
set <- set01


# DEFINE FACTOR FOR CONTRASTS ####
Felzartamab_Followup <- set$Felzartamab_Followup %>% droplevels()
cortex <- set$Cortexprob



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
fit_interaction_1 <- limma::lmFit(set, design)
cfit_interaction_1 <- limma::contrasts.fit(fit_interaction_1, contrast_interaction_01)
ebayes_interaction_1 <- limma::eBayes(cfit_interaction_1)
tab_interaction_1 <- limma::topTable(ebayes_interaction_1, adjust = "fdr", sort.by = "p", number = "all") %>%
    as_tibble(rownames = "AffyID") %>%
    right_join(affymap219 %>% dplyr::select(AffyID, Symb, Gene, PBT), ., by = "AffyID")

fit_felzartamab_1 <- limma::lmFit(set, design)
cfit_felzartamab_1 <- limma::contrasts.fit(fit_felzartamab_1, contrast_felzartamab_01)
ebayes_felzartamab_1 <- limma::eBayes(cfit_felzartamab_1)
tab_felzartamab_1 <- limma::topTable(ebayes_felzartamab_1, adjust = "fdr", sort.by = "p", number = "all") %>%
    as_tibble(rownames = "AffyID") %>%
    right_join(affymap219 %>% dplyr::select(AffyID, Symb, Gene, PBT), ., by = "AffyID")

fit_placebo_1 <- limma::lmFit(set, design)
cfit_placebo_1 <- limma::contrasts.fit(fit_placebo_1, contrast_placebo_01)
ebayes_placebo_1 <- limma::eBayes(cfit_placebo_1)
tab_placebo_1 <- limma::topTable(ebayes_placebo_1, adjust = "fdr", sort.by = "p", number = "all") %>%
    as_tibble(rownames = "AffyID") %>%
    right_join(affymap219 %>% dplyr::select(AffyID, Symb, Gene, PBT), ., by = "AffyID")


# FIT BLOCK week52 - week24 LIMMA MODEL ####
fit_interaction_2 <- limma::lmFit(set, design)
cfit_interaction_2 <- limma::contrasts.fit(fit_interaction_2, contrast_interaction_02)
ebayes_interaction_2 <- limma::eBayes(cfit_interaction_2)
tab_interaction_2 <- limma::topTable(ebayes_interaction_2, adjust = "fdr", sort.by = "p", number = "all") %>%
    as_tibble(rownames = "AffyID") %>%
    right_join(affymap219 %>% dplyr::select(AffyID, Symb, Gene, PBT), ., by = "AffyID")

fit_felzartamab_2 <- limma::lmFit(set, design)
cfit_felzartamab_2 <- limma::contrasts.fit(fit_felzartamab_2, contrast_felzartamab_02)
ebayes_felzartamab_2 <- limma::eBayes(cfit_felzartamab_2)
tab_felzartamab_2 <- limma::topTable(ebayes_felzartamab_2, adjust = "fdr", sort.by = "p", number = "all") %>%
    as_tibble(rownames = "AffyID") %>%
    right_join(affymap219 %>% dplyr::select(AffyID, Symb, Gene, PBT), ., by = "AffyID")

fit_placebo_2 <- limma::lmFit(set, design)
cfit_placebo_2 <- limma::contrasts.fit(fit_placebo_2, contrast_placebo_02)
ebayes_placebo_2 <- limma::eBayes(cfit_placebo_2)
tab_placebo_2 <- limma::topTable(ebayes_placebo_2, adjust = "fdr", sort.by = "p", number = "all") %>%
    as_tibble(rownames = "AffyID") %>%
    right_join(affymap219 %>% dplyr::select(AffyID, Symb, Gene, PBT), ., by = "AffyID")


# FIT BLOCK week52 - baseline LIMMA MODEL ####
fit_interaction_3 <- limma::lmFit(set, design)
cfit_interaction_3 <- limma::contrasts.fit(fit_interaction_3, contrast_interaction_03)
ebayes_interaction_3 <- limma::eBayes(cfit_interaction_3)
tab_interaction_3 <- limma::topTable(ebayes_interaction_3, adjust = "fdr", sort.by = "p", number = "all") %>%
    as_tibble(rownames = "AffyID") %>%
    right_join(affymap219 %>% dplyr::select(AffyID, Symb, Gene, PBT), ., by = "AffyID")

fit_felzartamab_3 <- limma::lmFit(set, design)
cfit_felzartamab_3 <- limma::contrasts.fit(fit_felzartamab_3, contrast_felzartamab_03)
ebayes_felzartamab_3 <- limma::eBayes(cfit_felzartamab_3)
tab_felzartamab_3 <- limma::topTable(ebayes_felzartamab_3, adjust = "fdr", sort.by = "p", number = "all") %>%
    as_tibble(rownames = "AffyID") %>%
    right_join(affymap219 %>% dplyr::select(AffyID, Symb, Gene, PBT), ., by = "AffyID")

fit_placebo_3 <- limma::lmFit(set, design)
cfit_placebo_3 <- limma::contrasts.fit(fit_placebo_3, contrast_placebo_03)
ebayes_placebo_3 <- limma::eBayes(cfit_placebo_3)
tab_placebo_3 <- limma::topTable(ebayes_placebo_3, adjust = "fdr", sort.by = "p", number = "all") %>%
    as_tibble(rownames = "AffyID") %>%
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
    left_join(means, by = "AffyID") %>%
    dplyr::select(
        AffyID, Symb, Gene, PBT,
        t, plogFC, flogFC, logFC, P.Value, adj.P.Val,
        all_of(colnames(means))
    ) %>%
    left_join(means_K1208_cell_panel %>% dplyr::select(-Symb, -Gene, -PBT), by = "AffyID") %>%
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
        t, plogFC, flogFC, logFC, P.Value, adj.P.Val,
        all_of(colnames(means))
    ) %>%
    left_join(means_K1208_cell_panel %>% dplyr::select(-Symb, -Gene, -PBT), by = "AffyID") %>%
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
        t, plogFC, flogFC, logFC, P.Value, adj.P.Val,
        all_of(colnames(means))
    ) %>%
    left_join(means_K1208_cell_panel %>% dplyr::select(-Symb, -Gene, -PBT), by = "AffyID") %>%
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
    ),
    table = list(
        table_interaction_1,
        table_interaction_2,
        table_interaction_3
    )
)




# EXPORT THE DATA AS .RData FILE ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/"
names(limma_tables$table) <- limma_tables$design
# save(limma_tables, file = paste(saveDir, "All_probes_cortex_corrected_limma_1208.RData", sep = ""))


# EXPORT THE DATA AS AN EXCEL SHEET ####
saveDir1 <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/output/"
openxlsx::write.xlsx(limma_tables$table,
    asTable = TRUE,
    file = paste(saveDir1, "All_probes_cortex_corrected_limma_1208_15Aug24",
        # Sys.Date(),
        # format(Sys.time(), "_%I%M%p"),
        ".xlsx",
        sep = ""
    )
)
