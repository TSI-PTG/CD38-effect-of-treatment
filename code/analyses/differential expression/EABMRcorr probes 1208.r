# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(flextable) # install.packages("flextable") #for table outputs
library(officer) # install.packages("officer")
library(openxlsx) # install.packages("openxlsx")
library(readxl) # install.packages("readxl")
# Bioconductor libraries
library(Biobase) # BiocManager::install("Biobase")
library(genefilter) # BiocManager::install("genefilter")
library(limma) # BiocManager::install("limma")
# Custom operators, functions, and datasets
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
# load reference set
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/data_expressionset_k1208.RData")
# load affymap
load("Z:/DATA/Datalocks/Other data/affymap219_21Oct2019_1306_JR.RData")
# load SCC data
simplefile <- read_excel("Z:/MISC/Phil/AA All papers in progress/A GC papers/0000 simple XL files/Kidney 5086/MASTER COPY K5086 SimpleCorrAAInjRej 5AAInjNR 7AARej.xlsx")


# DEFINE SEED ####
seed <- 42


# IQR FILTER THE DATA ####
selected <- simplefile %>%
    dplyr::select(Affy, SYMB, "pvalRej7AA4-EABMR") %>%
    arrange(`pvalRej7AA4-EABMR`) %>%
    distinct(SYMB, .keep_all = TRUE) %>%
    slice(1:20) %>%
    pull(Affy)


# DEFINE THE SET ####
set <- data_expressionset_k1208[
    selected,
    data_expressionset_k1208$Patient %nin% c(15, 18)
]
# set <- data_expressionset_k1208[,data_expressionset_k1208$Patient %nin% c(15, 18)]
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
    as.data.frame() %>%
    rownames_to_column(var = "AffyID") %>%
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
    as.data.frame() %>%
    rownames_to_column(var = "AffyID") %>%
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
# log2(148) - log2(166)
# log2(132) - log2(67)
# 2^(log2(148) - log2(166))
# 2^(log2(132) - log2(67))


table_block_3 <- tab_block_3 %>%
    as.data.frame() %>%
    rownames_to_column(var = "AffyID") %>%
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


# GLOBAL PARAMETERS FOR FLEXTABLES ####
header2 <- c(
    # "AffyID",
    "Symb", "Gene", "PBT",
    rep("Mean expression", 4),
    rep("Differential expression", 5)
)
header3 <- c(
    # "AffyID",
    "Symb", "Gene", "PBT",
    rep("Placebo", 2), rep("Felzartamab", 2),
    "\u394\nFC\n(Placebo)", "\u394\nFC\n(Felzartamab)", rep("Effect of treatment", 3)
)
header4_Baseline_week24 <- c(
    # "AffyID",
    "Symb", "Gene", "PBT",
    "Baseline\n(N=10)", "Week24\n(N=10)", "Baseline\n(N=10)", "Week24\n(N=10)",
    "\u394\nFC\n(Placebo)", "\u394\nFC\n(Felzartamab)", "\u394\u394\nFC", "\u394\u394\nP", "\u394\u394\nFDR"
)
header4_week24_week52 <- c(
    # "AffyID",
    "Symb", "Gene", "PBT",
    "Week24\n(N=10)", "Week52\n(N=10)", "Week24\n(N=10)", "Week52\n(N=10)",
    "\u394\nFC\n(Placebo)", "\u394\nFC\n(Felzartamab)", "\u394\u394\nFC", "\u394\u394\nP", "\u394\u394\nFDR"
)
header4_week52_Baseline <- c(
    # "AffyID",
    "Symb", "Gene", "PBT",
    "Week24\n(N=10)", "Week52\n(N=10)", "Week24\n(N=10)", "Week52\n(N=10)",
    "\u394\nFC\n(Placebo)", "\u394\nFC\n(Felzartamab)", "\u394\u394\nFC", "\u394\u394\nP", "\u394\u394\nFDR"
)

title1 <- paste("Table i. Effect of Felzartamab treatment on expression of top 20 genes correlated with EABMR in Baseline vs Week24 (by P-value)", sep = "")
title2 <- paste("Table i. Effect of Felzartamab treatment on expression of top 20 genes correlated with EABMR in Week24 vs Week52 (by P-value)", sep = "")
title3 <- paste("Table i. Effect of Felzartamab treatment on expression of top 20 genes correlated with EABMR in Baseline vs Week52 (by P-value)", sep = "")

cellWidths <- c(1, 4, 2, 1, 1, 1, 1, 1, rep(1.1, 4)) # for individual tables up or down
cellWidths %>% length()


# FORMAT FLEXTABLES ####
flextable_block_1 <- table_block_1 %>%
    dplyr::rename(FDR = adj.P.Val, Symbol = Symb) %>%
    dplyr::select(-AffyID, -logFC) %>%
    mutate(
        pFC = pFC %>% round(2),
        fFC =fFC %>% round(2),
        FC = FC %>% round(2)
    ) %>%
    mutate_at(
        vars(contains("p."), FDR),
        ~ ifelse(
            . < 0.01,
            formatC(., digits = 0, format = "e"),
            formatC(., digits = 3, format = "f")
        )
    ) %>%
    mutate(
        PBT = PBT %>% str_remove_all(",RAT|,Rej-RAT|,GRIT1|,GRIT2|,cIRIT"),
        Gene = Gene %>% str_remove("///.*")
    ) %>%
    flextable::flextable() %>%
    flextable::delete_part("header") %>%
    flextable::add_header_row(top = TRUE, values = header4_Baseline_week24) %>%
    flextable::add_header_row(top = TRUE, values = header3) %>%
    flextable::add_header_row(top = TRUE, values = header2) %>%
    # flextable::add_header_row(top = TRUE, values = header1) %>%
    flextable::add_header_row(values = rep(title1, ncol_keys(.))) %>%
    flextable::merge_v(part = "header") %>%
    flextable::merge_h(part = "header") %>%
    flextable::bg(bg = "white", part = "all") %>%
    flextable::border_remove() %>%
    flextable::border(border = fp_border(), part = "all") %>%
    flextable::align(part = "header", align = "center") %>%
    flextable::align(part = "body", align = "center") %>%
    flextable::valign(i = 2:4, valign = "bottom", part = "header") %>%
    flextable::padding(padding.left = 3, padding.bottom = 0, padding.top = 0) %>%
    flextable::font(fontname = "Arial", part = "all") %>%
    flextable::fontsize(size = 8, part = "body") %>%
    flextable::fontsize(size = 8, part = "header") %>%
    flextable::fontsize(i = 1, size = 12, part = "header") %>%
    flextable::bold(part = "header") %>%
    flextable::width(width = cellWidths, unit = "cm") %>%
    flextable::width(., width = dim(.)$widths * 33 / (flextable_dim(.)$widths), unit = "cm")

flextable_block_2 <- table_block_2 %>%
    dplyr::rename(FDR = adj.P.Val, Symbol = Symb) %>%
    dplyr::select(-AffyID, -logFC) %>%
    mutate(
        pFC = pFC %>% round(2),
        fFC =fFC %>% round(2),
        FC = FC %>% round(2)
    ) %>%
    mutate_at(
        vars(contains("p."), FDR),
        ~ ifelse(
            . < 0.01,
            formatC(., digits = 0, format = "e"),
            formatC(., digits = 3, format = "f")
        )
    ) %>%
    mutate(
        PBT = PBT %>% str_remove_all(",RAT|,Rej-RAT|,GRIT1|,GRIT2|,cIRIT"),
        Gene = Gene %>% str_remove("///.*")
    ) %>%
    flextable::flextable() %>%
    flextable::delete_part("header") %>%
    flextable::add_header_row(top = TRUE, values = header4_week24_week52) %>%
    flextable::add_header_row(top = TRUE, values = header3) %>%
    flextable::add_header_row(top = TRUE, values = header2) %>%
    # flextable::add_header_row(top = TRUE, values = header1) %>%
    flextable::add_header_row(values = rep(title2, ncol_keys(.))) %>%
    flextable::merge_v(part = "header") %>%
    flextable::merge_h(part = "header") %>%
    flextable::bg(bg = "white", part = "all") %>%
    flextable::border_remove() %>%
    flextable::border(border = fp_border(), part = "all") %>%
    flextable::align(part = "header", align = "center") %>%
    flextable::align(part = "body", align = "center") %>%
    flextable::valign(i = 2:4, valign = "bottom", part = "header") %>%
    flextable::padding(padding.left = 3, padding.bottom = 0, padding.top = 0) %>%
    flextable::font(fontname = "Arial", part = "all") %>%
    flextable::fontsize(size = 8, part = "body") %>%
    flextable::fontsize(size = 8, part = "header") %>%
    flextable::fontsize(i = 1, size = 12, part = "header") %>%
    flextable::bold(part = "header") %>%
    flextable::width(width = cellWidths, unit = "cm") %>%
    flextable::width(., width = dim(.)$widths * 33 / (flextable_dim(.)$widths), unit = "cm")

flextable_block_3 <- table_block_3 %>%
    dplyr::rename(FDR = adj.P.Val, Symbol = Symb) %>%
    dplyr::select(-AffyID, -logFC) %>%
    mutate(
        pFC = pFC %>% round(2),
        fFC =fFC %>% round(2),
        FC = FC %>% round(2)
    ) %>%
    mutate_at(
        vars(contains("p."), FDR),
        ~ ifelse(
            . < 0.01,
            formatC(., digits = 0, format = "e"),
            formatC(., digits = 3, format = "f")
        )
    ) %>%
    mutate(
        PBT = PBT %>% str_remove_all(",RAT|,Rej-RAT|,GRIT1|,GRIT2|,cIRIT"),
        Gene = Gene %>% str_remove("///.*")
    ) %>%
    flextable::flextable() %>%
    flextable::delete_part("header") %>%
    flextable::add_header_row(top = TRUE, values = header4_week52_Baseline) %>%
    flextable::add_header_row(top = TRUE, values = header3) %>%
    flextable::add_header_row(top = TRUE, values = header2) %>%
    # flextable::add_header_row(top = TRUE, values = header1) %>%
    flextable::add_header_row(values = rep(title3, ncol_keys(.))) %>%
    flextable::merge_v(part = "header") %>%
    flextable::merge_h(part = "header") %>%
    flextable::bg(bg = "white", part = "all") %>%
    flextable::border_remove() %>%
    flextable::border(border = fp_border(), part = "all") %>%
    flextable::align(part = "header", align = "center") %>%
    flextable::align(part = "body", align = "center") %>%
    flextable::valign(i = 2:4, valign = "bottom", part = "header") %>%
    flextable::padding(padding.left = 3, padding.bottom = 0, padding.top = 0) %>%
    flextable::font(fontname = "Arial", part = "all") %>%
    flextable::fontsize(size = 8, part = "body") %>%
    flextable::fontsize(size = 8, part = "header") %>%
    flextable::fontsize(i = 1, size = 12, part = "header") %>%
    flextable::bold(part = "header") %>%
    flextable::width(width = cellWidths, unit = "cm") %>%
    flextable::width(., width = dim(.)$widths * 33 / (flextable_dim(.)$widths), unit = "cm")



# PRINT THE FLEXTABLES ####
# flextable_block_1 %>% print(preview = "pptx")
# flextable_block_2 %>% print(preview = "pptx")
# flextable_block_3 %>% print(preview = "pptx")



# MERGE TABLES FOR SIMPLE EXPORT ####
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
        table_block_1,
        table_block_2,
        table_block_3
    )
)


# EXPORT THE DATA AS AN EXCEL SHEET ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/"
names(limma_tables$table) <- limma_tables$design
save(limma_tables, file = paste(saveDir, "EABMRcorr probes limma 1208.RData", sep = ""))

