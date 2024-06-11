# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(readxl) # install.packages("readxl")
library(circlize) # install.packages("circlize")
# load AFFYMAP
load("Z:/DATA/Datalocks/Other data/affymap219_21Oct2019_1306_JR.RData") # for labeling genes
# load KTB18 SIGNATURE MATRIX
KTB18 <- read.delim("Z:/MISC/Phil/AA All papers in progress/A GC papers/A a kidney cibersort/data/CIBERSORTx Signature Matrix/KTB18.txt")
# load KTB18 SIGNATURE MATRIX
LM22 <- read.delim("Z:/MISC/Phil/AA All papers in progress/A GC papers/A a kidney cibersort/data/CIBERSORTx Signature Matrix/LM22.txt")
# load in house cell panel
atagc <- read_excel("Z:/MISC/Patrick Gauthier/R/affymap219-CELL-PANEL/backup/UPDATED 2017 ANNOTATIONS - MASTERFILE - U133 HUMAN CELL PANEL - ALL PROBESETS (nonIQR) pfhptg.xlsx") # TBB with TBBParaRef
# load k1208 data
# if (!exists("K1208")) {
#     load("Z:/DATA/Datalocks/Kidney/K1208_Feb28_2022_JR.RData")
# }
# load mean expression in K1208
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/mean_expression_K1208_MMDx.RData")



# MEAN EXPRESSION BY PROBE ####
mean_exprs_by_probe <- means_K1208 %>%
    mutate(
        mean_exprs = means_K1208 %>%
            dplyr::select(-AffyID:-PBT) %>%
            rowMeans()
    ) %>%
    dplyr::slice_max(mean_exprs, by = Symb)

affy_genes <- mean_exprs_by_probe %>%
    dplyr::filter(Symb != "") %>%
    distinct(Symb, .keep_all = TRUE) %>%
    dplyr::select(AffyID, Symb, Gene, PBT)





# WRANGLE GENE EXPRESSION PROFILES ####
gep_KTB18 <- affy_genes %>%
    right_join(
        KTB18 %>%
            as_tibble() %>%
            dplyr::rename(Symb = NAME),
        by = "Symb"
    )

gep_LM22 <- affy_genes %>%
    right_join(
        LM22 %>%
            as_tibble() %>%
            dplyr::rename(Symb = Gene.symbol),
        by = "Symb"
    )

gep_atagc <- affy_genes %>%
    right_join(
        atagc %>%
            mutate(
                mean_exprs = atagc %>%
                    dplyr::select(-`Affy Probeset ID`:-Index) %>%
                    rowMeans()
            ) %>%
            dplyr::slice_max(mean_exprs, by = Index) %>%
            dplyr::select(-`Affy Probeset ID`, -`Gene Symbol`) %>%
            dplyr::rename(
                Symb = `Index`,
                NKcell = `NK cell`
            ),
        by = "Symb"
    )




# DEFINE QUANTILES FOR NK CELL AND TCELL SELECTIVITY ####
nk_quantile_KTB18 <- 0.95
tcell_quantile_KTB18 <- 0.33

nk_quantile_LM22 <- 0.90
tcell_quantile_LM22 <- 0.5

nk_quantile_atagc <- 0.95
tcell_quantile_atagc <- 0.05


# ISOLATE NK SELECTIVE GENES ####
nk_genes_atagc <- gep_atagc %>%
    dplyr::filter(
        NKcell > quantile(NKcell, nk_quantile_atagc, na.rm = TRUE)
    ) %>%
    dplyr::filter(
        CD4 < quantile(CD4, tcell_quantile_atagc, na.rm = TRUE),
        CD8 < quantile(CD8, tcell_quantile_atagc, na.rm = TRUE),
    ) %>%
    dplyr::select(
        AffyID, Symb,
        NKcell, CD4, CD8
    ) %>%
    left_join(
        mean_exprs_by_probe %>% dplyr::select(-AffyID, -Gene, -PBT, -mean_exprs),
        by = "Symb"
    ) %>%
    dplyr::filter(AffyID != "") %>%
    distinct(Symb, .keep_all = TRUE) %>%
    mutate(
        NKcell_inclusion = paste("NKcell > ", nk_quantile_atagc * 100, "th percentile", sep = ""),
        Tcell_exclusion = paste("CD4 and CD8 < ", tcell_quantile_atagc * 100, "th percentile", sep = ""),
    )

nk_genes_KTB18 <- gep_KTB18 %>%
    dplyr::filter(
        NKFCGR3Alow > quantile(NKFCGR3Alow, nk_quantile_KTB18) |
            NKFCGR3Ahigh > quantile(NKFCGR3Ahigh, nk_quantile_KTB18)
    ) %>%
    dplyr::filter(
        CD4naiveTcells < quantile(CD4naiveTcells, tcell_quantile_KTB18),
        CD4memoryTcells < quantile(CD4memoryTcells, tcell_quantile_KTB18),
        CD8effectorTcells < quantile(CD8effectorTcells, tcell_quantile_KTB18),
        CD8effmemTcells < quantile(CD8effmemTcells, tcell_quantile_KTB18),
        CD8temraTcells < quantile(CD8temraTcells, tcell_quantile_KTB18)
    ) %>%
    dplyr::select(
        AffyID, Symb,
        NKFCGR3Alow, NKFCGR3Ahigh,
        CD4naiveTcells, CD4memoryTcells, CD8effectorTcells, CD8effmemTcells, CD8temraTcells
    ) %>%
    left_join(
        mean_exprs_by_probe %>% dplyr::select(-AffyID, -Gene, -PBT, -mean_exprs),
        by = "Symb"
    ) %>%
    dplyr::filter(AffyID != "") %>%
    distinct(Symb, .keep_all = TRUE) %>%
    mutate(
        NKcell_inclusion = paste("NKFCGR3Alow or NKFCGR3Ahigh > ", nk_quantile_KTB18 * 100, "th percentile", sep = ""),
        Tcell_exclusion = paste("all CD4 and CD8 < ", tcell_quantile_KTB18 * 100, "th percentile", sep = ""),
    )


nk_genes_LM22 <- gep_LM22 %>%
    dplyr::filter(
        NK.cells.activated > quantile(NK.cells.activated, nk_quantile_LM22) |
            NK.cells.resting > quantile(NK.cells.resting, nk_quantile_LM22)
    ) %>%
    dplyr::filter(
        T.cells.CD8 < quantile(T.cells.CD8, tcell_quantile_LM22),
        T.cells.CD4.naive < quantile(T.cells.CD4.naive, tcell_quantile_LM22),
        T.cells.CD4.memory.resting < quantile(T.cells.CD4.memory.resting, tcell_quantile_LM22),
        T.cells.CD4.memory.activated < quantile(T.cells.CD4.memory.activated, tcell_quantile_LM22),
        # T.cells.follicular.helper < quantile(T.cells.follicular.helper, tcell_quantile_LM22),
        # T.cells.regulatory..Tregs. < quantile(T.cells.regulatory..Tregs., tcell_quantile_LM22),
        # T.cells.gamma.delta < quantile(T.cells.gamma.delta, tcell_quantile_LM22)
    ) %>%
    dplyr::select(
        AffyID, Symb,
        NK.cells.activated, NK.cells.resting,
        T.cells.CD4.naive, T.cells.CD4.memory.resting, T.cells.CD4.memory.activated,
        T.cells.CD8, T.cells.follicular.helper, T.cells.regulatory..Tregs., T.cells.gamma.delta
    ) %>%
    left_join(
        mean_exprs_by_probe %>% dplyr::select(-AffyID, -Gene, -PBT, -mean_exprs),
        by = "Symb"
    ) %>%
    dplyr::filter(AffyID != "") %>%
    distinct(Symb, .keep_all = TRUE) %>%
    mutate(
        NKcell_inclusion = paste("NK.cells.activated or NK.cells.resting > ", nk_quantile_LM22 * 100, "th percentile", sep = ""),
        Tcell_exclusion = paste("all CD4 and CD8 < ", tcell_quantile_LM22 * 100, "th percentile", sep = ""),
    )



# JOIN THE NK GENES FROM THE THREE PANELS ####
nk_genes <- tibble(
    panel = c("ATAGC_U133", "KTB18_RNAseq", "LM22_U133"),
    data = list(nk_genes_atagc, nk_genes_KTB18, nk_genes_LM22),
    thresholds = list(
        tibble(nk_quantile_atagc, tcell_quantile_atagc),
        tibble(nk_quantile_KTB18, tcell_quantile_KTB18),
        tibble(nk_quantile_LM22, tcell_quantile_LM22)
    )
)
names(nk_genes$data) <- nk_genes$panel


nk_genes %>% unnest(everything())



# DEFINE COLOR GRADIENT FOR HEATMAP ####
# max <- gep %>%
#     dplyr::filter(Symb %in% nk_genes) %>%
#     select_if(is.numeric) %>%
#     max()

# col_fun <- colorRamp2(c(0, max), c("#c9e3d8", "red"))


# CREATE HEATMAP ####
# gep %>%
#     dplyr::filter(Symb %in% nk_genes) %>%
#     data.frame() %>%
#     column_to_rownames("Symb") %>%
#     dplyr::select(-AffyID, -Gene) %>%
#     ComplexHeatmap::Heatmap(
#         col = col_fun,
#         cluster_columns = FALSE,
#         cluster_rows = FALSE,
#         heatmap_legend_param = list(
#             legend_height = unit(6, "cm"),
#             title = "expression"
#         )
#     )


# EXPORT THE GENE SET ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/"
save(nk_genes, file = paste(saveDir, "NK cell selective genes.RData", sep = ""))


# EXPORT THE DATA AS AN EXCEL SHEET ####
saveDir1 <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/output/"
openxlsx::write.xlsx(nk_genes$data,
    asTable = TRUE,
    file = paste(saveDir1, "NK_genes_11Jun24",
        # Sys.Date(),
        # format(Sys.time(), "_%I%M%p"),
        ".xlsx",
        sep = ""
    )
)
