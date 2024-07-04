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


# DEFINE CELL STATES OF INTEREST ####
cell_states <- c(
    "PT-New1",
    "PT-New2",
    "PT-New3",
    "PT-New4",
    "TAL-New1",
    "TAL-New2",
    "TAL-New3",
    "TAL-New4",
    "DCT-New1",
    "DCT-New2",
    "DCT-New3",
    "DCT-New4",
    "CNT-New1",
    "CNT-New2",
    "CNT-New3",
    "CD-PC-New1",
    "CD-PC-New2",
    "CD-IC-New1",
    "CD-IC-New2"
)


# DEFINE INJURY SELECTIVE GENES ####
genes_injury <- Hmisc::.q(
    IFITM3,
    LRP2,
    MET,
    MYO5B,
    NQO1,
    SLC2A1,
    SLC12A1,
    VCAM1,
    IGFBP7,
    ALDOB,
    SERPINA1,
    SPP1,
    LCN2,
    HAVCR1,
    NRF2,
    HIF1A,
    MYC,
    JUN
)


# WRANGLE THE INJURY MARKER DATA ####
injury_markers <- genes_injury_markers %>%
    unnest(data) %>%
    dplyr::filter(Symb %in% genes_injury) %>%
    dplyr::select(cluster:PBT) %>%
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
    dplyr::select(-cluster)



# WRANGLE THE SIMPLE FILE DATA ####
K4502 <- simplefile %>%
    dplyr::select(Affy, corrInjPCA1, corrInjPCA2, corrInjPCA3) %>%
    dplyr::rename(AffyID = Affy)

join <- injury_markers %>%
    left_join(K4502, by = "AffyID") %>%
    distinct(Symb, .keep_all = TRUE) %>%
    dplyr::select(-AffyID)



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

flextable  %>% print(preview = "pptx")
