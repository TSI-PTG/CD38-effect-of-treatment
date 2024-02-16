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
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0Georg Felz CD38 Vienna/G_Rstuff/data/vienna_1208_16Feb24.RData")
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
set00 <- vienna_1208[selected, ]


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
    "x =  (Group_FelzFU1_Felz-Group_FelzIndex_Felz)/2 -(Group_FelzFU1_NoFelz -Group_FelzIndex_NoFelz)/2",
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
    left_join(., means, by = "AffyID") %>%
    # dplyr::relocate(c("BLADpos Mean", "BLADneg Mean"), .before = t) %>%
    dplyr::select(-AveExpr, -B) %>%
    mutate(FC = 2^logFC, .after = logFC) %>%
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
    left_join(., means, by = "AffyID") %>%
    # dplyr::relocate(c("BLADpos Mean", "BLADneg Mean"), .before = t) %>%
    dplyr::select(-AveExpr, -B) %>%
    mutate(FC = 2^logFC, .after = logFC) %>%
    left_join(
        .,
        means_K5086 %>% dplyr::select(-Symb, -Gene, -PBT),
        by = "AffyID"
    )

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
header4 <- c(
    # "AffyID",
    "Symb", "Gene", "PBT", "FC", "P", "FDR",
    "Index\n(N=11)", "FU1\n(N=11)", "Index\n(N=11)", "FU1\n(N=11)",
    "NR\n(N=2476)", "Mixed\n(N=327)", "TCMR\n(N=326)", "pTCMR\n(N=74)", "ABMR\n(N=1242)", "pABMR\n(N=209)"
)

title <- paste("Table i. Top 20 differentially expressed transcripts in biopsies from treated vs untreated patients (by P-value)", sep = "")

cellWidths <- c(1.5, 5, 3, 1, 1, 1, rep(1.1, 10)) # for individual tables up or down
cellWidths %>% length()


# FORMAT FLEXTABLES ####
flextable_block_1 <- table_block_1 %>%
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
    flextable::add_header_row(top = TRUE, values = header4) %>%
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
    flextable::add_header_row(top = TRUE, values = header4) %>%
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



# MERGE TABLES FOR SIMPLE EXPORT ####
# limma_tables <- tibble(
#     design = c("interaction", "groupwise"),
#     table = list(table_block_1, table_block_2)
# )
limma_tables <- tibble(
    design = c("interaction"),
    table = list(table_block_1)
)


# EXPORT THE DATA AS AN EXCEL SHEET ####
savedir1 <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0Georg Felz CD38 Vienna/G_Rstuff/output/"

names(limma_tables$table) <- limma_tables$design
openxlsx::write.xlsx(limma_tables$table,
    asTable = TRUE,
    file = paste(savedir1, "IQR_filtered_probes_limma_vienna_1208_17Nov23",
        # Sys.Date(),
        # format(Sys.time(), "_%I%M%p"),
        ".xlsx",
        sep = ""
    )
)
