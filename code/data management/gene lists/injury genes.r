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



genes_injury_markers %>%
    mutate(data = map(data, filter, abs(log2FC) > 1)) %>%
    unnest(data) %>%
    drop_na(AffyID) %>%
    filter(Symb == "VIM") %>%
    dplyr::select(AffyID, Symb, cluster)




# WRANGLE THE INJURY MARKER DATA ####
injury_markers <- genes_injury_markers %>%
    mutate(data = map(data, filter, abs(log2FC) > 1)) %>%
    unnest(data) %>%
    drop_na(AffyID) %>%
    dplyr::filter(cluster %>% str_detect("New")) %>%
    dplyr::select(cluster, AffyID, Symb, Gene, PBT) %>%
    nest(.by = AffyID) %>%
    # mutate(
    #     "cellular expression" = map_chr(
    #         data,
    #         function(data) {
    #             data %>%
    #                 pull(cluster) %>%
    #                 paste(collapse = ",")
    #         }
    #     )
    # ) %>%
    unnest(data) %>%
    distinct(Symb, .keep_all = TRUE) %>%
    dplyr::select(-cluster)


injury_markers %>%
    left_join(K4502, by = "AffyID") %>%
    dplyr::filter(Symb == "VIM") %>%
    dplyr::select(AffyID, Symb, cluster)



# WRANGLE THE SIMPLE FILE DATA ####
simplefile %>%
    dplyr::filter(SYMB == "VIM") %>%
    dplyr::select(Affy, SYMB)


K4502 <- simplefile %>%
    dplyr::rename(
        AffyID = Affy,
        Symb = SYMB
    ) %>%
    # dplyr::slice_min(`pvalInjAA4-AKI2`, by = "SYMB") %>%
    dplyr::select(
        AffyID,
        `corrInjAA3-AKI1`, `pvalInjAA3-AKI1`, `corrInjAA4-AKI2`, `pvalInjAA4-AKI2`,
        corrInjPCA1, pvalInjPCA1, corrInjPCA2, pvalInjPCA2, corrInjPCA3, pvalInjPCA3,
    )






# JOIN THE INJURY MARKER DATA WITH THE SIMPLEFILE DATA ####
injury_markers %>%
    left_join(K4502, by = "AffyID") %>%
        dplyr::filter(Symb == "VIM") %>%
        dplyr::select(AffyID, Symb)




join <- injury_markers %>%
    left_join(K4502, by = "AffyID") %>%
    distinct(Symb, .keep_all = TRUE) %>%
    dplyr::select(-AffyID)


injury_markers %>%
    dplyr::filter(Symb == "VIM") %>%
    dplyr::select(AffyID, Symb)



join %>%
    dplyr::filter(Symb == "VIM") %>%
    dplyr::select(AffyID, Symb)

# DEFINE THE INJURY GENES BY AKI1 AND AKI2 ASSOCIATION ####
genes_injury <- join %>%
    dplyr::filter(`corrInjAA3-AKI1` < 0 & `corrInjAA4-AKI2` < 0)



genes_injury %>% print(n = "all")


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
