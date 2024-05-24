# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(ggrepel) # install.packages("ggrepel")
library(ggpubr) # install.packages("ggpubr")
library(patchwork) # install.packages("patchwork")
library(readxl) # install.packages("readxl")
# Custom operators, functions, and datasets
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
# load reference data
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/IQR_filtered_probes_limma_1208.RData")
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/felzartamab_ARTanova.RData")
# load SCC data
simplefile <- read_excel("Z:/MISC/Phil/AA All papers in progress/A GC papers/0000 simple XL files/Kidney 5086/MASTER COPY K5086 SimpleCorrAAInjRej 5AAInjNR 7AARej.xlsx")


# CHECK THE CATEGORIES ####
felzartamab_ARTanova$category


# DEFINE THE ABMR SCORES ####
scores_ABMR <- felzartamab_ARTanova %>%
    dplyr::filter(category %in% c("cfDNA", "ABMR", "injury", "archetypes", "rejectionPC", "injuryPC")) %>%
    dplyr::select(score, variable, medians_delta) %>%
    unnest(everything()) %>%
    dplyr::filter(Followup_pairwise == "Baseline - Week52") %>%
    distinct(variable, .keep_all = TRUE) %>%
    dplyr::select(score, variable, Followup_pairwise, median_delta_delta)

scores_ABMR$median_delta_delta[[1]]
