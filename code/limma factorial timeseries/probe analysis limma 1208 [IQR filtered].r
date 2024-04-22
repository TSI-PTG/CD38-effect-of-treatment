# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(flextable) # install.packages("flextable") #for table outputs
library(officer) # install.packages("officer")
library(openxlsx) # install.packages("openxlsx")
# Bioconductor libraries
library(Biobase) # BiocManager::install("Biobase")
library(genefilter) # BiocManager::install("genefilter")
library(limma) # BiocManager::install("limma")
# Custom operators, functions, and datasets
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
# load reference set
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0Georg Felz CD38 Vienna/G_Rstuff/data/vienna_1208_6MAr24.RData")
# load affymap
load("Z:/DATA/Datalocks/Other data/affymap219_21Oct2019_1306_JR.RData")
# load reference expression data
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0Georg Felz CD38 Vienna/G_Rstuff/data/mean_expression_K5086_MMDx.RData")


# IQR FILTER THE DATA ####
f1 <- function(x) (IQR(x) > 0.5)
ff <- filterfun(f1)
if (!exists("selected")) {
    selected <- genefilter(vienna_1208, ff)
}
set00 <- vienna_1208[selected, vienna_1208$STUDY_EVALUATION_ID %nin% c(15, 18)]


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
        Group = Group %>% factor(levels = c("Index", "FU1", "FU2")),
        Felz = Felz %>% factor(labels = c("NoFelz", "Felz")),
        TxBx = TxBx %>% as.numeric(),
        Group_Felz = paste(Group, Felz, sep = "_")
         %>%
            factor(
                levels = c("Index_NoFelz", "FU1_NoFelz","FU2_NoFelz", "Index_Felz", "FU1_Felz", "FU2_Felz"),
                labels = c("Index_NoFelz", "FU1_NoFelz", "FU2_NoFelz", "Index_Felz", "FU1_Felz", "FU2_Felz")
            )
    ) %>%
    dplyr::select(CEL, Patient, Center, Group, Felz, Group_Felz, TxBx)



# DEFINE SUBSETS FOR PAIRWISE CONTRASTS ####
set_Index_FU1 <- set[, set$Group != "FU2"]
set_FU1_FU2 <- set[, set$Group != "Index"]
set_Index_FU2 <- set[, set$Group != "FU1"]


# BLOCK DESIGN FU1 - Index ####
Group_Felz <- set_Index_FU1$Group_Felz  %>% droplevels
design_block01 <- model.matrix(~ 0 + Group_Felz)
contrast_block_01 <- makeContrasts(
    "x =  (Group_FelzFU1_Felz-Group_FelzIndex_Felz)/2 - (Group_FelzFU1_NoFelz-Group_FelzIndex_NoFelz)/2",
    levels = design_block01
)


# BLOCK DESIGN FU2 - FU1 ####
Group_Felz <- set_FU1_FU2$Group_Felz %>% droplevels()
design_block02 <- model.matrix(~ 0 + Group_Felz)
contrast_block_02 <- makeContrasts(
    "x =  (Group_FelzFU2_Felz-Group_FelzFU1_Felz)/2 - (Group_FelzFU2_NoFelz -Group_FelzFU1_NoFelz)/2",
    levels = design_block02
)


# BLOCK DESIGN FU2 - Index ####
Group_Felz <- set_Index_FU2$Group_Felz %>% droplevels()
design_block03 <- model.matrix(~ 0 + Group_Felz)
contrast_block_03 <- makeContrasts(
    "x =  (Group_FelzFU2_Felz-Group_FelzIndex_Felz)/2 - (Group_FelzFU2_NoFelz -Group_FelzIndex_NoFelz)/2",
    levels = design_block03
)


# FIT BLOCK FU1 - Index LIMMA MODEL ####
fit_block_1 <- limma::lmFit(set_Index_FU1, design_block01)
cfit_block_1 <- limma::contrasts.fit(fit_block_1, contrast_block_01)
ebayes_block_1 <- limma::eBayes(cfit_block_1)
tab_block_1 <- limma::topTable(ebayes_block_1, adjust = "fdr", sort.by = "p", number = "all")
ebayes_block_1 %>% limma::topTable()


# FIT BLOCK FU2 - FU1 LIMMA MODEL ####
fit_block_2 <- limma::lmFit(set_FU1_FU2, design_block02)
cfit_block_2 <- limma::contrasts.fit(fit_block_2, contrast_block_02)
ebayes_block_2 <- limma::eBayes(cfit_block_2)
tab_block_2 <- limma::topTable(ebayes_block_2, adjust = "fdr", sort.by = "p", number = "all")
ebayes_block_2 %>% limma::topTable()


# FIT BLOCK FU2 - Index LIMMA MODEL ####
fit_block_3 <- limma::lmFit(set_Index_FU2, design_block03)
cfit_block_3 <- limma::contrasts.fit(fit_block_3, contrast_block_03)
ebayes_block_3 <- limma::eBayes(cfit_block_3)
tab_block_3 <- limma::topTable(ebayes_block_3, adjust = "fdr", sort.by = "p", number = "all")
ebayes_block_3 %>% limma::topTable()


# CALCULATE MEAN GENE EXPRESSION FOR EACH PROBE BETWEEN GROUPINGS ####
means_Index_FU1 <- fit_block_1 %>%
    avearrays() %>%
    data.frame() %>%
    rownames_to_column("AffyID") %>%
    tibble() %>%
    mutate_if(is.numeric, ~ 2^. %>% round(0)) %>%
    rename_at(vars(contains("Felz")), ~ str_remove(., "Group_Felz"))

means_FU2_FU1 <- fit_block_2 %>%
    avearrays() %>%
    data.frame() %>%
    rownames_to_column("AffyID") %>%
    tibble() %>%
    mutate_if(is.numeric, ~ 2^. %>% round(0)) %>%
    rename_at(vars(contains("Felz")), ~ str_remove(., "Group_Felz"))

means_FU2_Index <- fit_block_3 %>%
    avearrays() %>%
    data.frame() %>%
    rownames_to_column("AffyID") %>%
    tibble() %>%
    mutate_if(is.numeric, ~ 2^. %>% round(0)) %>%
    rename_at(vars(contains("Felz")), ~ str_remove(., "Group_Felz"))



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
    left_join(., means_Index_FU1, by = "AffyID") %>%
    # dplyr::relocate(c("BLADpos Mean", "BLADneg Mean"), .before = t) %>%
    dplyr::select(-AveExpr, -B) %>%
    mutate(
        FC = 2^logFC, 
        FCfelz = FU1_Felz/Index_Felz,
        .after = logFC) %>%
    left_join(
        .,
        means_K5086 %>% dplyr::select(-Symb, -Gene, -PBT),
        by = "AffyID"
    )

table_block_2 <- tab_block_2 %>%
    as.data.frame() %>%
    rownames_to_column(var = "AffyID") %>%
    right_join(affymap219 %>% dplyr::select(AffyID, Symb, Gene, PBT), ., by = "AffyID") %>%
    arrange(P.Value) %>%
    # mutate_at(c("logFC", "t"), ~ round(., 1)) %>%
    # mutate_at(c("P.Value", "adj.P.Val"), ~ format(., digits = 1)) %>%
    mutate_at(c("P.Value", "adj.P.Val"), as.numeric) %>%
    tibble() %>%
    left_join(., means_FU2_FU1, by = "AffyID") %>%
    # dplyr::relocate(c("BLADpos Mean", "BLADneg Mean"), .before = t) %>%
    dplyr::select(-AveExpr, -B) %>%
    mutate(FC = 2^logFC, .after = logFC) %>%
    left_join(
        .,
        means_K5086 %>% dplyr::select(-Symb, -Gene, -PBT),
        by = "AffyID"
    )

table_block_3 <- tab_block_3 %>%
    as.data.frame() %>%
    rownames_to_column(var = "AffyID") %>%
    right_join(affymap219 %>% dplyr::select(AffyID, Symb, Gene, PBT), ., by = "AffyID") %>%
    arrange(P.Value) %>%
    # mutate_at(c("logFC", "t"), ~ round(., 1)) %>%
    # mutate_at(c("P.Value", "adj.P.Val"), ~ format(., digits = 1)) %>%
    mutate_at(c("P.Value", "adj.P.Val"), as.numeric) %>%
    tibble() %>%
    left_join(., means_FU2_Index, by = "AffyID") %>%
    # dplyr::relocate(c("BLADpos Mean", "BLADneg Mean"), .before = t) %>%
    dplyr::select(-AveExpr, -B) %>%
    mutate(FC = 2^logFC, .after = logFC) %>%
    left_join(
        .,
        means_K5086 %>% dplyr::select(-Symb, -Gene, -PBT),
        by = "AffyID"
    )

# GLOBAL PARAMETERS FOR FLEXTABLES ####
header1 <- c(
    # "AffyID",
    "Symb", "Gene", "PBT", "FC", "P", "FDR",
    rep("Felzartamab study", 4),
    rep("K5086 reference set", 6)
)
header2 <- c(
    # "AffyID",
    "Symb", "Gene", "PBT", "FC", "P", "FDR",
    rep("Mean expression by group", 4),
    rep("Mean expression by MMDx", 6)
)
header3 <- c(
    # "AffyID",
    "Symb", "Gene", "PBT", "FC", "P", "FDR",
    rep("NoFelz", 2), rep("Felz", 2),
    rep("Mean expression by MMDx", 6)
)
header4_Index_FU1 <- c(
    # "AffyID",
    "Symb", "Gene", "PBT", "FC", "P", "FDR",
    "Index\n(N=10)", "FU1\n(N=10)", "Index\n(N=10)", "FU1\n(N=10)",
    "NR\n(N=2476)", "Mixed\n(N=327)", "TCMR\n(N=326)", "pTCMR\n(N=74)", "ABMR\n(N=1242)", "pABMR\n(N=209)"
)
header4_FU2_FU1 <- c(
    # "AffyID",
    "Symb", "Gene", "PBT", "FC", "P", "FDR",
    "FU1\n(N=10)", "FU2\n(N=10)", "FU1\n(N=10)", "FU2\n(N=10)",
    "NR\n(N=2476)", "Mixed\n(N=327)", "TCMR\n(N=326)", "pTCMR\n(N=74)", "ABMR\n(N=1242)", "pABMR\n(N=209)"
)
header4_FU2_Index <- c(
    # "AffyID",
    "Symb", "Gene", "PBT", "FC", "P", "FDR",
    "FU1\n(N=10)", "FU2\n(N=10)", "FU1\n(N=10)", "FU2\n(N=10)",
    "NR\n(N=2476)", "Mixed\n(N=327)", "TCMR\n(N=326)", "pTCMR\n(N=74)", "ABMR\n(N=1242)", "pABMR\n(N=209)"
)

title <- paste("Table i. Top 20 differentially expressed transcripts in biopsies from treated vs untreated patients (by P-value)", sep = "")

cellWidths <- c(1.5, 5, 3, 1, 1, 1, rep(1.1, 10)) # for individual tables up or down
cellWidths %>% length()


table_block_1 %>%
    arrange(FCfelz)



# FORMAT FLEXTABLES ####
flextable_block_1 <- table_block_1 %>%
    arrange(FCfelz) %>%
    dplyr::rename(FDR = adj.P.Val, Symbol = Symb) %>%
    dplyr::select(-AffyID, -logFC, -t, -FCfelz) %>%
    # distinct(Symbol, .keep_all = T) %>%
    dplyr::slice(1:20) %>%
    mutate(FC = FC %>% round(2)) %>%
    mutate_at(
        vars(contains("p."), FDR),
        ~ ifelse(
            . < 0.01,
            formatC(., digits = 0, format = "e"),
            formatC(., digits = 3, format = "f")
        )
    ) %>%
    flextable::flextable() %>%
    flextable::delete_part("header") %>%
    flextable::add_header_row(top = TRUE, values = header4_Index_FU1) %>%
    flextable::add_header_row(top = TRUE, values = header3) %>%
    flextable::add_header_row(top = TRUE, values = header2) %>%
    flextable::add_header_row(top = TRUE, values = header1) %>%
    flextable::add_header_row(values = rep(title, ncol_keys(.))) %>%
    flextable::merge_v(part = "header") %>%
    flextable::merge_h(part = "header") %>%
    # flextable::bg(bg = "grey90", part = "header") %>%
    flextable::bg(bg = "white", part = "all") %>%
    flextable::border_remove() %>%
    flextable::border(border = fp_border(), part = "all") %>%
    flextable::align(part = "header", align = "center") %>%
    flextable::align(part = "body", align = "center") %>%
    flextable::valign(i = 2:4, valign = "bottom", part = "header") %>%
    flextable::padding(padding.left = 3, padding.bottom = 0, padding.top = 0) %>%
    flextable::font(fontname = "Arial", part = "all") %>%
    flextable::fontsize(size = 7, part = "body") %>%
    flextable::fontsize(size = 8, part = "header") %>%
    flextable::fontsize(i = 1, size = 12, part = "header") %>%
    flextable::bold(part = "header") %>%
    flextable::width(width = cellWidths, unit = "cm") %>%
    flextable::width(., width = dim(.)$widths * 33 / (flextable_dim(.)$widths), unit = "cm")

flextable_block_2 <- table_block_2 %>%
    dplyr::rename(FDR = adj.P.Val, Symbol = Symb) %>%
    dplyr::select(-AffyID, -logFC, -t) %>%
    # distinct(Symbol, .keep_all = T) %>%
    dplyr::slice(1:20) %>%
    mutate(FC = FC %>% round(2)) %>%
    mutate_at(
        vars(contains("p."), FDR),
        ~ ifelse(
            . < 0.01,
            formatC(., digits = 0, format = "e"),
            formatC(., digits = 3, format = "f")
        )
    ) %>%
    flextable::flextable() %>%
    flextable::delete_part("header") %>%
    flextable::add_header_row(top = TRUE, values = header4_FU2_FU1) %>%
    flextable::add_header_row(top = TRUE, values = header3) %>%
    flextable::add_header_row(top = TRUE, values = header2) %>%
    flextable::add_header_row(top = TRUE, values = header1) %>%
    flextable::add_header_row(values = rep(title, ncol_keys(.))) %>%
    flextable::merge_v(part = "header") %>%
    flextable::merge_h(part = "header") %>%
    # flextable::bg(bg = "grey90", part = "header") %>%
    flextable::bg(bg = "white", part = "all") %>%
    flextable::border_remove() %>%
    flextable::border(border = fp_border(), part = "all") %>%
    flextable::align(part = "header", align = "center") %>%
    flextable::align(part = "body", align = "center") %>%
    flextable::valign(i = 2:4, valign = "bottom", part = "header") %>%
    flextable::padding(padding.left = 3, padding.bottom = 0, padding.top = 0) %>%
    flextable::font(fontname = "Arial", part = "all") %>%
    flextable::fontsize(size = 7, part = "body") %>%
    flextable::fontsize(size = 8, part = "header") %>%
    flextable::fontsize(i = 1, size = 12, part = "header") %>%
    flextable::bold(part = "header") %>%
    flextable::width(width = cellWidths, unit = "cm") %>%
    flextable::width(., width = dim(.)$widths * 33 / (flextable_dim(.)$widths), unit = "cm")

flextable_block_3 <- table_block_3 %>%
    dplyr::rename(FDR = adj.P.Val, Symbol = Symb) %>%
    dplyr::select(-AffyID, -logFC, -t) %>%
    # distinct(Symbol, .keep_all = T) %>%
    dplyr::slice(1:20) %>%
    mutate(FC = FC %>% round(2)) %>%
    mutate_at(
        vars(contains("p."), FDR),
        ~ ifelse(
            . < 0.01,
            formatC(., digits = 0, format = "e"),
            formatC(., digits = 3, format = "f")
        )
    ) %>%
    flextable::flextable() %>%
    flextable::delete_part("header") %>%
    flextable::add_header_row(top = TRUE, values = header4_FU2_Index) %>%
    flextable::add_header_row(top = TRUE, values = header3) %>%
    flextable::add_header_row(top = TRUE, values = header2) %>%
    flextable::add_header_row(top = TRUE, values = header1) %>%
    flextable::add_header_row(values = rep(title, ncol_keys(.))) %>%
    flextable::merge_v(part = "header") %>%
    flextable::merge_h(part = "header") %>%
    # flextable::bg(bg = "grey90", part = "header") %>%
    flextable::bg(bg = "white", part = "all") %>%
    flextable::border_remove() %>%
    flextable::border(border = fp_border(), part = "all") %>%
    flextable::align(part = "header", align = "center") %>%
    flextable::align(part = "body", align = "center") %>%
    flextable::valign(i = 2:4, valign = "bottom", part = "header") %>%
    flextable::padding(padding.left = 3, padding.bottom = 0, padding.top = 0) %>%
    flextable::font(fontname = "Arial", part = "all") %>%
    flextable::fontsize(size = 7, part = "body") %>%
    flextable::fontsize(size = 8, part = "header") %>%
    flextable::fontsize(i = 1, size = 12, part = "header") %>%
    flextable::bold(part = "header") %>%
    flextable::width(width = cellWidths, unit = "cm") %>%
    flextable::width(., width = dim(.)$widths * 33 / (flextable_dim(.)$widths), unit = "cm")

# PRINT THE FLEXTABLES ####
flextable_block_1 %>% print(preview = "pptx")
flextable_block_2 %>% print(preview = "pptx")
flextable_block_3 %>% print(preview = "pptx")



# MERGE TABLES FOR SIMPLE EXPORT ####
# limma_tables <- tibble(
#     design = c("interaction", "groupwise"),
#     table = list(table_block_1, table_block_2)
# )
limma_tables <- tibble(
    design = c(
        "Index_vs_FU1",
        "FU1_vs_FU2",
        "Index_vs_FU2"
        ),
    table = list(
        table_block_1,
        table_block_2,
        table_block_3
    )
)


# EXPORT THE DATA AS AN EXCEL SHEET ####
savedir1 <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0Georg Felz CD38 Vienna/G_Rstuff/output/"

names(limma_tables$table) <- limma_tables$design
# openxlsx::write.xlsx(limma_tables$table,
#     asTable = TRUE,
#     file = paste(savedir1, "IQR_filtered_probes_limma_vienna_1208_7Mar24",
#         # Sys.Date(),
#         # format(Sys.time(), "_%I%M%p"),
#         ".xlsx",
#         sep = ""
#     )
# )
