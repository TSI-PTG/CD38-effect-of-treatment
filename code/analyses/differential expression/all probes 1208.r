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
set  %>% pData  %>% colnames


# WRANGLE THE GROUPING NAMES ####

set$Felzartamab_Followup  %>% factor



# BLOCK DESIGN #1 ####
Felzartamab_Followup <- set$Felzartamab_Followup
design_block <- model.matrix(~ 0 + Felzartamab_Followup)
contrast_block_01 <- makeContrasts(
    "x =  (Felzartamab_FollowupBaseline_Placebo-Felzartamab_FollowupWeek24_Placebo)/2 - (Felzartamab_FollowupBaseline_Felzartamab-Felzartamab_FollowupWeek24_Felzartamab)/2",
    levels = design_block
)



# BLOCK DESIGN #2 ####
Felzartamab_Followup <- set$Felzartamab_Followup
design_block <- model.matrix(~ 0 + Felzartamab_Followup)
contrast_block_2 <- makeContrasts(
    "x =  Felzartamab_FollowupFU1_Felz - (Felzartamab_FollowupIndex_Felz + Felzartamab_FollowupIndex_NoFelz + Felzartamab_FollowupFU1_NoFelz)/3",
    levels = design_block
)







# FACTORIAL DESIGN ####
Group <- set$Group
Felz <- set$Felz
Patient <- set$Patient
contrasts(Group) <- contr.sum(2)
contrasts(Felz) <- contr.sum(2)
design_factorial <- model.matrix(~ Group * Felz)


# FIT FACTORIAL LIMMA MODEL ####
fit_factorial <- limma::lmFit(set, design_factorial)
ebayes_factorial <- limma::eBayes(fit_factorial)
tab_factorial <- limma::topTable(ebayes_factorial, coef = 3, adjust = "fdr", sort.by = "p", number = "all")
ebayes_factorial %>% limma::topTable()


# FIT BLOCK #1 LIMMA MODEL ####
fit_block_1 <- limma::lmFit(set, design_block)
cfit_block_1 <- limma::contrasts.fit(fit_block_1, contrast_block_01)
ebayes_block_1 <- limma::eBayes(cfit_block_1)
tab_block_1 <- limma::topTable(ebayes_block_1, adjust = "fdr", sort.by = "p", number = "all")
ebayes_block_1 %>% limma::topTable()


# FIT BLOCK #1 LIMMA MODEL ####
fit_block_2 <- limma::lmFit(set, design_block)
cfit_block_2 <- limma::contrasts.fit(fit_block_2, contrast_block_2)
ebayes_block_2 <- limma::eBayes(cfit_block_2)
tab_block_2 <- limma::topTable(ebayes_block_2, adjust = "fdr", sort.by = "p", number = "all")


# CALCULATE MEAN GENE EXPRESSION FOR EACH PROBE BETWEEN GROUPINGS ####
means <- fit_block_1 %>%
    avearrays() %>%
    data.frame() %>%
    rownames_to_column("AffyID") %>%
    tibble() %>%
    mutate_if(is.numeric, ~ 2^. %>% round(0)) %>%
    rename_at(vars(contains("Felz")), ~str_remove(., "Felzartamab_Followup"))


# FORMAT TOPTABLES ####
table_block_1 <- tab_block_1 %>%
    as.data.frame() %>%
    rownames_to_column(var = "AffyID") %>%
    right_join(affymap219 %>% dplyr::select(AffyID, Symb, Gene, PBT), ., by = "AffyID") %>%
    arrange(P.Value) %>%
    # mutate_at(c("logFC", "t"), ~ round(., 1)) %>%
    # mutate_at(c("P.Value", "adj.P.Val"), ~ format(., digits = 1)) %>%
    mutate_at(c("P.Value", "adj.P.Val"), as.numeric) %>%
    tibble() %>%
    left_join(., means, by = "AffyID") %>%
    # dplyr::relocate(c("BLADpos Mean", "BLADneg Mean"), .before = t) %>%
    dplyr::select(-AveExpr, -B) %>%
    mutate(FC = 2^logFC, .after = logFC)

table_block_2 <- tab_block_2 %>%
    as.data.frame() %>%
    rownames_to_column(var = "AffyID") %>%
    right_join(affymap219 %>% dplyr::select(AffyID, Symb, Gene, PBT), ., by = "AffyID") %>%
    arrange(P.Value) %>%
    # mutate_at(c("logFC", "t"), ~ round(., 1)) %>%
    # mutate_at(c("P.Value", "adj.P.Val"), ~ format(., digits = 1)) %>%
    mutate_at(c("P.Value", "adj.P.Val"), as.numeric) %>%
    tibble() %>%
    left_join(., means, by = "AffyID") %>%
    # dplyr::relocate(c("BLADpos Mean", "BLADneg Mean"), .before = t) %>%
    dplyr::select(-AveExpr, -B) %>%
    mutate(FC = 2^logFC, .after = logFC)

table_factorial <- tab_factorial %>%
    as.data.frame() %>%
    rownames_to_column(var = "AffyID") %>%
    right_join(affymap219 %>% dplyr::select(AffyID, Symb, Gene, PBT), ., by = "AffyID") %>%
    arrange(P.Value) %>%
    # mutate_at(c("logFC", "t"), ~ round(., 1)) %>%
    # mutate_at(c("P.Value", "adj.P.Val"), ~ format(., digits = 1)) %>%
    mutate_at(c("P.Value", "adj.P.Val"), as.numeric) %>%
    tibble() %>%
    left_join(., means, by = "AffyID") %>%
    # dplyr::relocate(c("BLADpos Mean", "BLADneg Mean"), .before = t) %>%
    dplyr::select(-AveExpr, -B) %>%
    mutate(FC = 2^logFC, .after = logFC)


# MERGE TABLES FOR SIMPLE EXPORT ####
limma_tables <- tibble(
    design = c("interaction", "absolute difference"),
    table = list(table_block_1, table_block_2) 
)



# EXPORT THE DATA AS AN EXCEL SHEET ####
savedir1 <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0Georg Felz CD38 Vienna/G_Rstuff/output/"

names(limma_tables$table) <- limma_tables$design
openxlsx::write.xlsx(limma_tables$table,
    asTable = TRUE,
    file = paste(savedir1, "all_probes_limma_Vienna44_16Nov23 ",
        # Sys.Date(),
        # format(Sys.time(), "_%I%M%p"),
        ".xlsx",
        sep = ""
    )
)
