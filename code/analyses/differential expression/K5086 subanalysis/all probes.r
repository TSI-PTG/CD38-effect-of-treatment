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
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/K5086Jul12_2024_PTG.RData")
# load affymap
load("Z:/DATA/Datalocks/Other data/affymap219_21Oct2019_1306_JR.RData")
affymap219 <- affymap219 %>% tibble()
# load in house cell panel
atagc <- read_excel("Z:/MISC/Patrick Gauthier/R/affymap219-CELL-PANEL/backup/UPDATED 2017 ANNOTATIONS - MASTERFILE - U133 HUMAN CELL PANEL - ALL PROBESETS (nonIQR) pfhptg.xlsx")


# WRANGLE THE PHENOTYPE DATA ####
pData(K5086) <- K5086 %>%
    pData() %>%
    tibble() %>%
    mutate(
        RejAA7 = case_when(
            RejAA7Clust == 1 ~ "NR",
            RejAA7Clust == 2 ~ "TCMR1",
            RejAA7Clust == 3 ~ "TCMR2",
            RejAA7Clust == 4 ~ "EABMR",
            RejAA7Clust == 5 ~ "FABMR",
            RejAA7Clust == 6 ~ "LABMR",
            RejAA7Clust == 7 ~ "Minor"
        ) %>%
            factor(levels = c("NR", "EABMR", "FABMR", "LABMR", "TCMR1", "TCMR2", "Minor"))
    )


# DEFINE SEED ####
seed <- 42
set.seed(seed)


# DEFINE THE SET ####
set00 <- K5086


# IQR FILTER THE DATA ####
# f1 <- function(x) (IQR(x) > 0.5)
# ff <- filterfun(f1)
# if (!exists("selected")) {
#     selected <- genefilter(set00, ff)
# }
# set01 <- set00[selected, ]
set01 <- set00


# WRANGLE THE CELL PANEL DATA ####
cell_panel <- atagc %>%
    dplyr::select(
        `Affy Probeset ID`, `Index`,
        `Mcrphg unstim`, `Mcrphg + IFNg`,
        `Unstim HUVEC`, `HUVEC + IFNg`,
        `Unstim RPTEC`, `RPTEC + IFNg`
    ) %>%
    dplyr::rename(
        AffyID_U133 = `Affy Probeset ID`,
        Symb = Index,
        `Macrophage (unstimulated)` = `Mcrphg unstim`,
        `Macrophage (IFNg stimulated)` = `Mcrphg + IFNg`,
        `HUVEC (unstimulated)` = `Unstim HUVEC`,
        `HUVEC (IFNg stimulated)` = `HUVEC + IFNg`,
        `RPTEC (unstimulated)` = `Unstim RPTEC`,
        `RPTEC (IFNg stimulated)` = `RPTEC + IFNg`
    ) %>%
    dplyr::slice_max(`HUVEC (unstimulated)`, by = "Symb", with_ties = FALSE)

# means_K1208_cell_panel <- means_K1208 %>%
#     left_join(
#         cell_panel %>%
#             dplyr::select(-AffyID_U133),
#         by = "Symb"
#     )


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
    dplyr::filter(Symb != "") %>%
    distinct(Symb, .keep_all = TRUE) %>%
    pull(AffyID)

# genes <- mean_exprs_1208 %>%
#     dplyr::slice_max(mean_expression, by = "Symb") %>%
#     dplyr::filter(Symb != "", AffyID %in% genes_baseline) %>%
#     distinct(Symb, .keep_all = TRUE) %>%
#     pull(AffyID)

# set <- set01[featureNames(set01) %in% genes, ]
set <- set01


# DEFINE FACTOR FOR CONTRASTS ####
RejAA7 <- K5086$RejAA7 %>% droplevels()


# DESIGN ####
design <- model.matrix(~ 0 + RejAA7)


# CONTRAST DESIGN FABMR ####
contrast_FABMRvEE <- makeContrasts(
    "x =  RejAA7FABMR - (RejAA7NR+RejAA7TCMR1+RejAA7TCMR2+RejAA7EABMR+RejAA7LABMR+RejAA7Minor)/6",
    levels = design
)
contrast_FABMRvNONABMR <- makeContrasts(
    "x =  RejAA7FABMR - (RejAA7NR+RejAA7TCMR1+RejAA7TCMR2+RejAA7Minor)/4",
    levels = design
)
contrast_FABMRvNR <- makeContrasts(
    "x =  RejAA7FABMR - RejAA7NR",
    levels = design
)


# CONTRAST DESIGN EABMR ####
contrast_EABMRvEE <- makeContrasts(
    "x =  RejAA7EABMR - (RejAA7NR+RejAA7TCMR1+RejAA7TCMR2+RejAA7FABMR+RejAA7LABMR+RejAA7Minor)/6",
    levels = design
)
contrast_EABMRvNONABMR <- makeContrasts(
    "x =  RejAA7EABMR - (RejAA7NR+RejAA7TCMR1+RejAA7TCMR2+RejAA7Minor)/4",
    levels = design
)
contrast_EABMRvNR <- makeContrasts(
    "x =  RejAA7EABMR - RejAA7NR",
    levels = design
)


# CONTRAST DESIGN ALL ABMR ####
contrast_ABMRvNONABMR <- makeContrasts(
    "x =  (RejAA7EABMR+RejAA7FABMR+RejAA7LABMR)/3 - (RejAA7NR+RejAA7TCMR1+RejAA7TCMR2+RejAA7Minor)/4",
    levels = design
)
contrast_ABMRvNR <- makeContrasts(
    "x =   (RejAA7EABMR+RejAA7FABMR+RejAA7LABMR)/3 - RejAA7NR",
    levels = design
)


# FIT LIMMA MODEL ####
fit <- limma::lmFit(set, design)


# MAKE FABMR CONTRASTS ####
cfit_FABMRvNR <- limma::contrasts.fit(fit, contrast_FABMRvNR)
ebayes_FABMRvNR <- limma::eBayes(cfit_FABMRvNR)
tab_FABMRvNR <- limma::topTable(ebayes_FABMRvNR, adjust = "fdr", sort.by = "p", number = "all") %>%
    as_tibble(rownames = "AffyID") %>%
    right_join(affymap219 %>% dplyr::select(AffyID, Symb, Gene, PBT), ., by = "AffyID") %>%
        arrange(P.Value, logFC %>% desc())

cfit_FABMRvNONABMR <- limma::contrasts.fit(fit, contrast_FABMRvNONABMR)
ebayes_FABMRvNONABMR <- limma::eBayes(cfit_FABMRvNONABMR)
tab_FABMRvNONABMR <- limma::topTable(ebayes_FABMRvNONABMR, adjust = "fdr", sort.by = "p", number = "all") %>%
    as_tibble(rownames = "AffyID") %>%
    right_join(affymap219 %>% dplyr::select(AffyID, Symb, Gene, PBT), ., by = "AffyID") %>%
        arrange(P.Value, logFC %>% desc())


# MAKE EABMR CONTRASTS ####
cfit_EABMRvNR <- limma::contrasts.fit(fit, contrast_EABMRvNR)
ebayes_EABMRvNR <- limma::eBayes(cfit_EABMRvNR)
tab_EABMRvNR <- limma::topTable(ebayes_EABMRvNR, adjust = "fdr", sort.by = "p", number = "all") %>%
    as_tibble(rownames = "AffyID") %>%
    right_join(affymap219 %>% dplyr::select(AffyID, Symb, Gene, PBT), ., by = "AffyID") %>%
        arrange(P.Value, logFC %>% desc())

cfit_EABMRvNONABMR <- limma::contrasts.fit(fit, contrast_EABMRvNONABMR)
ebayes_EABMRvNONABMR <- limma::eBayes(cfit_EABMRvNONABMR)
tab_EABMRvNONABMR <- limma::topTable(ebayes_EABMRvNONABMR, adjust = "fdr", sort.by = "p", number = "all") %>%
    as_tibble(rownames = "AffyID") %>%
    right_join(affymap219 %>% dplyr::select(AffyID, Symb, Gene, PBT), ., by = "AffyID") %>%
        arrange(P.Value, logFC %>% desc())


# MAKE ALL ABMR CONTRASTS ####
cfit_ABMRvNR <- limma::contrasts.fit(fit, contrast_ABMRvNR)
ebayes_ABMRvNR <- limma::eBayes(cfit_ABMRvNR)
tab_ABMRvNR <- limma::topTable(ebayes_ABMRvNR, adjust = "fdr", sort.by = "p", number = "all") %>%
    as_tibble(rownames = "AffyID") %>%
    right_join(affymap219 %>% dplyr::select(AffyID, Symb, Gene, PBT), ., by = "AffyID") %>%
        arrange(P.Value, logFC %>% desc())

cfit_ABMRvNONABMR <- limma::contrasts.fit(fit, contrast_ABMRvNONABMR)
ebayes_ABMRvNONABMR <- limma::eBayes(cfit_ABMRvNR)
tab_ABMRvNONABMR <- limma::topTable(ebayes_ABMRvNONABMR, adjust = "fdr", sort.by = "p", number = "all") %>%
    as_tibble(rownames = "AffyID") %>%
    right_join(affymap219 %>% dplyr::select(AffyID, Symb, Gene, PBT), ., by = "AffyID") %>%
        arrange(P.Value, logFC %>% desc())


# CALCULATE MEAN GENE EXPRESSION FOR EACH PROBE BETWEEN GROUPINGS ####
means <- fit %>%
    avearrays() %>%
    as_tibble(rownames = "AffyID") %>%
    mutate_if(is.numeric, ~ 2^. %>% round(0)) %>%
    rename_at(vars(contains("RejAA7")), ~ str_remove(., "RejAA7")) %>%
    dplyr::select(-contains("Week12"), -any_of(c("cortex")))


# FORMAT FABMR TOPTABLES ####
table_FABMRvNR <- tab_FABMRvNR %>%
    mutate_at(c("P.Value", "adj.P.Val"), as.numeric) %>%
    left_join(means, by = "AffyID") %>%
    dplyr::select(
        AffyID, Symb, Gene, PBT,
        t, logFC, P.Value, adj.P.Val,
        all_of(colnames(means))
    ) %>%
    left_join(cell_panel %>% dplyr::select(-AffyID_U133), by = "Symb") %>%
    dplyr::rename(
        p = P.Value,
        FDR = adj.P.Val
    )
table_FABMRvNONABMR <- tab_FABMRvNONABMR %>%
    mutate_at(c("P.Value", "adj.P.Val"), as.numeric) %>%
    left_join(means, by = "AffyID") %>%
    dplyr::select(
        AffyID, Symb, Gene, PBT,
        t, logFC, P.Value, adj.P.Val,
        all_of(colnames(means))
    ) %>%
    left_join(cell_panel %>% dplyr::select(-AffyID_U133), by = "Symb") %>%
    dplyr::rename(
        p = P.Value,
        FDR = adj.P.Val
    )


# FORMAT EABMR TOPTABLES ####
table_EABMRvNR <- tab_EABMRvNR %>%
    mutate_at(c("P.Value", "adj.P.Val"), as.numeric) %>%
    left_join(means, by = "AffyID") %>%
    dplyr::select(
        AffyID, Symb, Gene, PBT,
        t, logFC, P.Value, adj.P.Val,
        all_of(colnames(means))
    ) %>%
    left_join(cell_panel %>% dplyr::select(-AffyID_U133), by = "Symb") %>%
    dplyr::rename(
        p = P.Value,
        FDR = adj.P.Val
    )
table_EABMRvNONABMR <- tab_EABMRvNONABMR %>%
    mutate_at(c("P.Value", "adj.P.Val"), as.numeric) %>%
    left_join(means, by = "AffyID") %>%
    dplyr::select(
        AffyID, Symb, Gene, PBT,
        t, logFC, P.Value, adj.P.Val,
        all_of(colnames(means))
    ) %>%
    left_join(cell_panel %>% dplyr::select(-AffyID_U133), by = "Symb") %>%
    dplyr::rename(
        p = P.Value,
        FDR = adj.P.Val
    )


# FORMAT ALL ABMR TOPTABLES ####
table_ABMRvNR <- tab_ABMRvNR %>%
    mutate_at(c("P.Value", "adj.P.Val"), as.numeric) %>%
    left_join(means, by = "AffyID") %>%
    dplyr::select(
        AffyID, Symb, Gene, PBT,
        t, logFC, P.Value, adj.P.Val,
        all_of(colnames(means))
    ) %>%
    left_join(cell_panel %>% dplyr::select(-AffyID_U133), by = "Symb") %>%
    dplyr::rename(
        p = P.Value,
        FDR = adj.P.Val
    )
table_ABMRvNONABMR <- tab_ABMRvNONABMR %>%
    mutate_at(c("P.Value", "adj.P.Val"), as.numeric) %>%
    left_join(means, by = "AffyID") %>%
    dplyr::select(
        AffyID, Symb, Gene, PBT,
        t, logFC, P.Value, adj.P.Val,
        all_of(colnames(means))
    ) %>%
    left_join(cell_panel %>% dplyr::select(-AffyID_U133), by = "Symb") %>%
    dplyr::rename(
        p = P.Value,
        FDR = adj.P.Val
    )


limma_tables <- tibble(
    design = c(
        "FABMRvNR",
        "FABMRvNONABMR",
        "EABMRvNR",
        "EABMRvNONABMR",
        "ABMRvNR",
        "ABMRvNONABMR"
    ),
    table = list(
        table_FABMRvNR,
        table_FABMRvNONABMR,
        table_EABMRvNR,
        table_EABMRvNONABMR,
        table_ABMRvNR,
        table_ABMRvNONABMR
    )
)
limma_tables$table[[3]]



# EXPORT THE DATA AS .RData FILE ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/"
names(limma_tables$table) <- limma_tables$design
save(limma_tables, file = paste(saveDir, "limma_K5086.RData", sep = ""))



# EXPORT THE DATA AS AN EXCEL SHEET ####
names(limma_tables$table) <- limma_tables$design
saveDir1 <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/output/"
openxlsx::write.xlsx(limma_tables$table,
    asTable = TRUE,
    file = paste(saveDir1, "limma_K5086_19July24",
        # Sys.Date(),
        # format(Sys.time(), "_%I%M%p"),
        ".xlsx",
        sep = ""
    )
)




limma_tables %>%
    dplyr::filter(design == "ABMRvNONABMR") %>%
    pull(table) %>%
    pluck(1) %>%
    dplyr::filter(p < 0.05) %>%
    mutate(direction = ifelse(logFC < 0, "down", "up")) %>%
    nest(.by = direction)

