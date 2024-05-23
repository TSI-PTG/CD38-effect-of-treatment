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



# BLOCK DESIGN week24 - baseline ####
design_block01 <- model.matrix(~ 0 + Felzartamab_Followup)
contrast_block_01 <- makeContrasts(
    "x =  (Felzartamab_FollowupWeek24_Felzartamab-Felzartamab_FollowupBaseline_Felzartamab)/2 -(Felzartamab_FollowupWeek24_Placebo-Felzartamab_FollowupBaseline_Placebo)/2",
    levels = design_block01
)


# BLOCK DESIGN week52 - week24 ####
design_block02 <- model.matrix(~ 0 + Felzartamab_Followup)
contrast_block_02 <- makeContrasts(
    "x =  (Felzartamab_FollowupWeek52_Felzartamab-Felzartamab_FollowupWeek24_Felzartamab)/2 - (Felzartamab_FollowupWeek52_Placebo-Felzartamab_FollowupWeek24_Placebo)/2",
    levels = design_block02
)


# BLOCK DESIGN week52 - baseline####
design_block03 <- model.matrix(~ 0 + Felzartamab_Followup)
contrast_block_03 <- makeContrasts(
    "x =  (Felzartamab_FollowupWeek52_Felzartamab-Felzartamab_FollowupBaseline_Felzartamab)/2 - (Felzartamab_FollowupWeek52_Placebo-Felzartamab_FollowupBaseline_Placebo)/2",
    levels = design_block03
)




# FIT BLOCK week24 - baseline LIMMA MODEL ####
fit_block_1 <- limma::lmFit(set, design_block01)
cfit_block_1 <- limma::contrasts.fit(fit_block_1, contrast_block_01)
ebayes_block_1 <- limma::eBayes(cfit_block_1)
tab_block_1 <- limma::topTable(ebayes_block_1, adjust = "fdr", sort.by = "p", number = "all")
ebayes_block_1 %>% limma::topTable()


# FIT BLOCK week52 - week24 LIMMA MODEL ####
fit_block_2 <- limma::lmFit(set, design_block02)
cfit_block_2 <- limma::contrasts.fit(fit_block_2, contrast_block_02)
ebayes_block_2 <- limma::eBayes(cfit_block_2)
tab_block_2 <- limma::topTable(ebayes_block_2, adjust = "fdr", sort.by = "p", number = "all")
ebayes_block_2 %>% limma::topTable()


# FIT BLOCK week52 - baseline LIMMA MODEL ####
fit_block_3 <- limma::lmFit(set, design_block03)
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
    dplyr::select(-contains("Week52"), -contains("Week12"))

means_week24_week52 <- fit_block_2 %>%
    avearrays() %>%
    data.frame() %>%
    rownames_to_column("AffyID") %>%
    tibble() %>%
    mutate_if(is.numeric, ~ 2^. %>% round(0)) %>%
    rename_at(vars(contains("Felz")), ~ str_remove(., "Felzartamab_Followup")) %>%
    dplyr::select(-contains("Baseline"), -contains("Week12"))

means_week52_baseline <- fit_block_3 %>%
    avearrays() %>%
    data.frame() %>%
    rownames_to_column("AffyID") %>%
    tibble() %>%
    mutate_if(is.numeric, ~ 2^. %>% round(0)) %>%
    rename_at(vars(contains("Felz")), ~ str_remove(., "Felzartamab_Followup")) %>%
    dplyr::select(-contains("Week24"), -contains("Week12"))


# FORMAT TOPTABLES ####
table_block_1 <- tab_block_1 %>%
    as_tibble(rownames = "AffyID") %>%
    right_join(affymap219 %>% dplyr::select(AffyID, Symb, Gene, PBT), ., by = "AffyID") %>%
    arrange(P.Value) %>%
    mutate_at(c("P.Value", "adj.P.Val"), as.numeric) %>%
    tibble() %>%
    left_join(., means_baseline_week24, by = "AffyID") %>%
    dplyr::select(
        AffyID, Symb, Gene, PBT,
        all_of(colnames(means_baseline_week24)[-1]),
        logFC, P.Value, adj.P.Val,
    ) %>%
    mutate(
        pFC = 2^(log2(Week24_Placebo) - log2(Baseline_Placebo)) %>% round(2),
        fFC = 2^(log2(Week24_Felzartamab) - log2(Baseline_Felzartamab)) %>% round(2),
        FC = 2^logFC,
        .after = logFC
    )

table_block_2 <- tab_block_2 %>%
    as_tibble(rownames = "AffyID") %>%
    right_join(affymap219 %>% dplyr::select(AffyID, Symb, Gene, PBT), ., by = "AffyID") %>%
    arrange(P.Value) %>%
    mutate_at(c("P.Value", "adj.P.Val"), as.numeric) %>%
    tibble() %>%
    left_join(., means_week24_week52, by = "AffyID") %>%
    dplyr::select(
        AffyID, Symb, Gene, PBT,
        all_of(colnames(means_week24_week52)[-1]),
        logFC, P.Value, adj.P.Val,
    ) %>%
    mutate(
        pFC = 2^(log2(Week52_Placebo) - log2(Week24_Placebo)) %>% round(2),
        fFC = 2^(log2(Week52_Felzartamab) - log2(Week24_Felzartamab)) %>% round(2),
        FC = 2^logFC,
        .after = logFC
    )

table_block_3 <- tab_block_3 %>%
    as_tibble(rownames = "AffyID") %>%
    right_join(affymap219 %>% dplyr::select(AffyID, Symb, Gene, PBT), ., by = "AffyID") %>%
    arrange(P.Value) %>%
    mutate_at(c("P.Value", "adj.P.Val"), as.numeric) %>%
    tibble() %>%
    left_join(., means_week52_baseline, by = "AffyID") %>%
    dplyr::select(
        AffyID, Symb, Gene, PBT,
        all_of(colnames(means_week52_baseline)[-1]),
        logFC, P.Value, adj.P.Val,
    ) %>%
    mutate(
        pFC = 2^(log2(Week52_Placebo) - log2(Baseline_Placebo)) %>% round(2),
        fFC = 2^(log2(Week52_Felzartamab) - log2(Baseline_Felzartamab)) %>% round(2),
        FC = 2^logFC,
        .after = logFC
    )


limma_tables <- tibble(
    design = c(
        "Baseline_vs_Week24",
        "Week24_vs_Week52",
        "Baseline_vs_Week52"
    ),
    toptable = list(
        tab_block_1 %>%
            as_tibble(rownames = "AffyID"),
        tab_block_2 %>%
            as_tibble(rownames = "AffyID"),
        tab_block_3 %>%
            as_tibble(rownames = "AffyID")
    ),
    table = list(
        table_block_1 %>%
            relocate(logFC, .before = "FC") %>%
            dplyr::rename(
                "\u394 placebo FC" = pFC,
                "\u394 felz FC" = fFC,
                "\u394\u394 logFC" = logFC,
                "\u394\u394 FC" = FC,
                "\u394\u394 p" = P.Value,
                "\u394\u394 FDR" = adj.P.Val
            ),
        table_block_2 %>%
            relocate(logFC, .before = "FC") %>%
            dplyr::rename(
                "\u394 placebo FC" = pFC,
                "\u394 felz FC" = fFC,
                "\u394\u394 logFC" = logFC,
                "\u394\u394 FC" = FC,
                "\u394\u394 p" = P.Value,
                "\u394\u394 FDR" = adj.P.Val
            ),
        table_block_3 %>%
            relocate(logFC, .before = "FC") %>%
            dplyr::rename(
                "\u394 placebo FC" = pFC,
                "\u394 felz FC" = fFC,
                "\u394\u394 logFC" = logFC,
                "\u394\u394 FC" = FC,
                "\u394\u394 p" = P.Value,
                "\u394\u394 FDR" = adj.P.Val
            )
    )
)


limma_tables$table[[1]]
tab_block_1 %>% as_tibble(rownames = "AffyID")



# EXPORT THE DATA AS .RData FILE ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/"
names(limma_tables$table) <- limma_tables$design
save(limma_tables, file = paste(saveDir, "all probes limma 1208.RData", sep = ""))




# EXPORT THE DATA AS AN EXCEL SHEET ####
saveDir1 <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/output/"
openxlsx::write.xlsx(limma_tables$table,
    asTable = TRUE,
    file = paste(saveDir1, "all_probes_limma_1208_23May24",
        # Sys.Date(),
        # format(Sys.time(), "_%I%M%p"),
        ".xlsx",
        sep = ""
    )
)
