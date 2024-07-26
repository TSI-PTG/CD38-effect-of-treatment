# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(flextable) # install.packages("flextable")
library(officer) # install.packages("officer")
library(openxlsx) # install.packages("openxlsx")
library(readxl) # install.packages("readxl")
# Custom operators, functions, and datasets
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
# load affymap
load("Z:/DATA/Datalocks/Other data/affymap219_21Oct2019_1306_JR.RData")
# load gene lists
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/Hinze_injury_markers.RData")
# load K4502 injury simple file
simplefile_path <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/0000 simple XL files/Kidney 4502/MASTER COPY K4502 Injset SimpleCorrAAInjRej 5AAInj 7AARej.xlsx"
simplefile <- read_excel(path = simplefile_path, sheet = "simpleCorrAAInjRejInjset")


# DEFINE GENES FROM FUNCTIONAL ENRICHMENT ####
genes_injury <- Hmisc::.q(
    SPTBN1, NRP1, KLF6, 'HLA-DRB1', DOCK11, DUSP1, CHSY1, LDHA,
    CDC42SE2, TPM4, PTBP3, ARHGAP29, SYNE2, SERP1, CORO1C, LRRFIP1, PPFIBP1, CANX, JUN,
    SCOC, RNF213, SPP1, PIK3AP1, VIM, SVIL, TPM3, RAN, ANXA5, ADAM10, HSPA1A, CEBPB,
    SPARC, ACTB, HNRNPK, SLFN5, RHOA
)





genes_injury_markers %>%
    # mutate(data = map(data, filter, abs(log2FC) > 1)) %>%
    unnest(data) %>%
    drop_na(AffyID) %>%
    filter(Symb == "SLFN5") %>%
    dplyr::select(AffyID, Symb, cluster, log2FC)


# COLLATE ALL INJURY MARKER DATA ####
injury_markers_all <- genes_injury_markers %>%
    # mutate(data = map(data, filter, abs(log2FC) > 1)) %>%
    unnest(data) %>%
    drop_na(AffyID) %>%
    # dplyr::filter(cluster %>% str_detect("New")) %>%
    dplyr::select(cluster, AffyID, Symb, Gene, PBT) %>%
    nest(.by = AffyID) %>%
    mutate(
        "cellular expression in AKI" = map_chr(
            data,
            function(data) {
                data %>%
                    pull(cluster) %>%
                    paste(collapse = ",")
            }
        )
    ) %>%
    unnest(data) %>%
    distinct(AffyID, .keep_all = TRUE) %>%
    dplyr::select(-cluster)



# WRANGLE THE INJURY MARKER DATA ####
injury_markers <- genes_injury_markers %>%
    # mutate(data = map(data, filter, abs(log2FC) > 1)) %>%
    unnest(data) %>%
    drop_na(AffyID) %>%
    dplyr::filter(cluster %>% str_detect("New")) %>%
    dplyr::select(cluster, AffyID, Symb, Gene, PBT) %>%
    nest(.by = AffyID) %>%
    mutate(
        "cellular expression" = map_chr(
            data,
            function(data) {
                data %>%
                    pull(cluster) %>%
                    paste(collapse = ",")
            }
        )
    ) %>%
    unnest(data) %>%
    # distinct(Symb, .keep_all = TRUE) %>%
    dplyr::select(-cluster)


# WRANGLE THE SIMPLE FILE DATA ####
simplefile %>%
    dplyr::filter(SYMB == "VIM") %>%
    dplyr::select(Affy, SYMB)


K4502 <- simplefile %>%
    dplyr::rename(
        AffyID = Affy,
        Symb = SYMB
    ) %>%
    dplyr::select(
        AffyID,
        `corrInjAA3-AKI1`, `pvalInjAA3-AKI1`, `corrInjAA4-AKI2`, `pvalInjAA4-AKI2`,
        corrInjPCA1, pvalInjPCA1, corrInjPCA2, pvalInjPCA2, corrInjPCA3, pvalInjPCA3,
    )






# JOIN THE INJURY MARKER DATA WITH THE SIMPLEFILE DATA ####
join <- injury_markers %>%
    left_join(K4502, by = "AffyID") %>%
    dplyr::select(-AffyID) %>%
    dplyr::slice_max(`corrInjAA4-AKI2`, by = "Symb") %>%
    distinct(Symb, .keep_all = TRUE)





# DEFINE THE INJURY GENES BY AKI1 AND AKI2 ASSOCIATION ####
genes_injury <- join %>%
    dplyr::filter(
        # `corrInjAA3-AKI1` > 0,
        # `pvalInjAA3-AKI1` < 0.05,
        # `corrInjAA4-AKI2` > 0,
        # `pvalInjAA4-AKI2` < 0.05
    ) %>%
    mutate(
        AKI1_rank = rank(`corrInjAA3-AKI1` %>% desc(), ties.method = "average"),
        AKI2_rank = rank(`corrInjAA4-AKI2` %>% desc(), ties.method = "average"),
        PC1_rank = rank(corrInjPCA1 %>% desc(), ties.method = "average"),
        PC2_rank = rank(corrInjPCA2, ties.method = "average"),
        rank = rank(AKI2_rank + PC2_rank, ties.method = "average"),
        .before = 1
    ) %>%
    arrange(rank)



genes_injury %>%
    print(n = 100)



genes_injury %>%
    # dplyr::select(AKI2_rank, Symb:`pvalInjAA4-AKI2`) %>%
    dplyr::filter(Symb %in% Hmisc::.q(
        SLFN5, CYR61, SPARC, CEBPB, ADAM10,
        VIM, LRRFIP1, TPM4, DUSP1, DOCK11, "HLA-DRB1",
        KLF6, NRP1, RARRES3, ARL6IP5, FGFR1, MMP7, TSC22D1, WWC1, "HLA-E"
    ))




genes_injury <- join %>%
    # slice_min(`corrInjAA3-AKI1`, n = 100) %>%
    slice_min(`pvalInjAA3-AKI1`, n = 20)





# SAVE THE DATA ####
# saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/"
# save(genes_injury, file = paste(saveDir, "injury_genes.RData", sep = ""))


# MAKE A FLEXTABLE OF THE DATA ####
flextable <- join %>%
    flextable::flextable() %>%
    flextable::border(part = "header", border = officer::fp_border()) %>%
    flextable::border(part = "body", border = officer::fp_border()) %>%
    flextable::align(align = "center") %>%
    flextable::align(align = "center", part = "header") %>%
    flextable::font(fontname = "Arial", part = "all") %>%
    flextable::fontsize(size = 8, part = "all") %>%
    flextable::fontsize(size = 8, part = "footer") %>%
    flextable::fontsize(i = 1, size = 12, part = "header") %>%
    flextable::bold(part = "header") %>%
    flextable::bg(bg = "white", part = "all") %>%
    flextable::padding(padding = 0, part = "all")

flextable %>% print(preview = "pptx")



# EXPORT DATA FOR ALL MARKER GENES ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/"
openxlsx::write.xlsx(injury_markers_all,
    asTable = TRUE,
    file = paste(saveDir, "Injury_marker_genes_26Jul24",
        # Sys.Date(),
        # format(Sys.time(), "_%I%M%p"),
        ".xlsx",
        sep = ""
    )
)
