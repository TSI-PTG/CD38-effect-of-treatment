# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(ggbeeswarm) # install.packages("ggbeeswarm")
library(ggpubr) # install.packages("ggpubr")
library(ggh4x) # install.packages("ggh4x")
library(broom) # install.packages("broom") #for tabular model object transformations
library(flextable) # install.packages("flextable") #for table outputs
library(officer) # install.packages("officer")
# library(emmeans) # install.packages("emmeans") #for post-hoc testing and CLD
# library(multcomp) # install.packages("multcomp") #for for CLD
library(furrr) # install.packages("furrr")
library(progressr) # install.packages("progressr")
library(parallel) # install.packages("parallel")
library(openxlsx) # install.packages("openxlsx")
# Bioconductor libraries
library(Biobase) # BiocManager::install("Biobase")
library(genefilter) # BiocManager::install("genefilter")
library(limma) # BiocManager::install("limma")
library(biobroom) # BiocManager::install("biobroom")
# Custom operators, functions, and datasets
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
se <- function(x) sd(x) / sqrt(length((x)))
log10zero <- scales::trans_new(
    name = "log10zero",
    transform = function(x) log10(x + 0.001),
    inverse = function(x) 10^x - 0.001
)
# SET UP PROGRESS BAR HANDLER
plan(sequential)
handlers(global = TRUE)
handlers(list(handler_progress(
    format = ":message [:current/:total] [:bar] [:percent in :elapsed] [:eta remaining]",
    clear = FALSE,
    intrusiveness = 1,
    complete = "+"
)))
# Suppress pesky dplyr summarise info
options(dplyr.summarise.inform = FALSE)
# laod reference set
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0Georg CD38 Vienna/G_Rstuff/data/Vienna44_18Oct23.RData")
# load affymap
load("Z:/DATA/Datalocks/Other data/affymap219_21Oct2019_1306_JR.RData")


# IQR FILTER THE DATA ####
f1 <- function(x) (IQR(x) > 0.5)
ff <- filterfun(f1)
if (!exists("selected")) {
    selected <- genefilter(Vienna44, ff)
}
set00 <- Vienna44[selected, ]


# DEFINE SEED ####
seed <- 42


# DEFINE THE SET ####
set <- set00


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
    "x =  (Group_FelzIndex_Felz + Group_FelzFU1_Felz)/2 - (Group_FelzIndex_NoFelz + Group_FelzFU1_NoFelz)/2",
    levels = design_block
)


# BLOCK DESIGN #2 ####
Group_Felz <- set$Group_Felz
design_block <- model.matrix(~ 0 + Group_Felz)
contrast_block_2 <- makeContrasts(
    "x =  (Group_FelzIndex_Felz + Group_FelzIndex_NoFelz + Group_FelzFU1_NoFelz)/3 - Group_FelzFU1_Felz",
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
tab_factorial <- limma::topTable(ebayes_factorial, coef = 3, adjust = "fdr", sort.by = "p", resort.by = "p", number = "all")
ebayes_factorial %>% limma::topTable()



# FIT BLOCK #1 LIMMA MODEL ####
fit_block_1 <- limma::lmFit(set, design_block)
cfit_block_1 <- limma::contrasts.fit(fit_block_1, contrast_block_01)
ebayes_block_1 <- limma::eBayes(cfit_block_1)
tab_block_1 <- limma::topTable(ebayes_block_1, adjust = "fdr", sort.by = "p", number = "all")


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


# # EXPORT THE DATA AS AN EXCEL SHEET ####
# savedir1 <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0Georg CD38 Vienna/G_Rstuff/output/"

# names(res_art_table$data) <- res_art_table$direction
# openxlsx::write.xlsx(res_art_table$data,
#     asTable = TRUE,
#     file = paste(savedir1, "all_probes_ANOVAs_Vienna44_18Oct23 ",
#         # Sys.Date(),
#         # format(Sys.time(), "_%I%M%p"),
#         ".xlsx",
#         sep = ""
#     )
# )
