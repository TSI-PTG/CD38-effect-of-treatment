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
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0Georg Felz CD38 Vienna/G_Rstuff/data/Vienna44_18Oct23.RData")
# load affymap
load("Z:/DATA/Datalocks/Other data/affymap219_21Oct2019_1306_JR.RData")


# DEFINE SEED ####
seed <- 42


# DEFINE THE SET ####
set <- Vienna44


# WRANGLE THE PHENOTYPE DATA ####
pData(set) <- set %>%
    pData() %>%
    dplyr::rename(
        Patient = STUDY_EVALUATION_ID,
        Felz = Felzartamab_presumed
    ) %>%
    mutate(
        Patient = Patient %>% factor(),
        Group = Group %>% factor(levels = c("Index", "FU1")),
        Felz = Felz %>% factor(labels = c("NoFelz", "Felz")),
        TxBx = TxBx %>% as.numeric(),
        Group_Felz = paste(Group, Felz, sep = ":") %>%
            factor(
                levels = c("Index:NoFelz", "FU1:NoFelz", "Index:Felz", "FU1:Felz"),
                labels = c("Index_NoFelz", "FU1_NoFelz", "Index_Felz", "FU1_Felz")
            )
    ) %>%
    dplyr::select(CEL, Patient, Center, Group, Felz, Group_Felz, TxBx)



# BLOCK DESIGN #1 ####
Group_Felz <- set$Group_Felz
design_block <- model.matrix(~ 0 + Group_Felz)
contrast_block_01 <- makeContrasts(
    "x =  (Group_FelzIndex_NoFelz-Group_FelzFU1_NoFelz)/2 - (Group_FelzIndex_Felz-Group_FelzFU1_Felz)/2",
    levels = design_block
)


# BLOCK DESIGN #2 ####
Group_Felz <- set$Group_Felz
design_block <- model.matrix(~ 0 + Group_Felz)
contrast_block_2 <- makeContrasts(
    "x =  Group_FelzFU1_Felz - (Group_FelzIndex_Felz + Group_FelzIndex_NoFelz + Group_FelzFU1_NoFelz)/3",
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
    rename_at(vars(contains("Felz")), ~str_remove(., "Group_Felz"))


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
