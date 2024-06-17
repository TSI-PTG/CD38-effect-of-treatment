# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(gghalves) # install.packages("gghalves")
library(ggrepel)
# Bioconductor libraries
library(Biobase) # BiocManager::install("Biobase")
library(limma) # BiocManager::install("limma")
# Custom operators, functions, and datasets
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
# load reference set
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/data_expressionset_k1208.RData")
# load affymap
load("Z:/DATA/Datalocks/Other data/affymap219_21Oct2019_1306_JR.RData")
# load gene lists
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/genes_NK_GEP.RData")
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/NK_genes_L765.RData")
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/ABMR_endothelial_genes.RData")
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/ABMR_activity_genes.RData")
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/IFNG_genes.RData")


# DEFINE THE PROBES ####
probes <- genes_ABMR_endothelial %>%
    slice_max(means, by = "Symb") %>%
    pull(AffyID)



# DEFINE THE SET ####
set <- data_expressionset_k1208[, data_expressionset_k1208$Patient %nin% c(15, 18)]


# WRANGLE THE PHENOTYPE DATA ####
data_phenotype <- set %>%
    pData() %>%
    dplyr::select(CEL, Patient, Felzartamab, Followup, Felzartamab_Followup)


# WRANGLE THE EXPRESSION DATA ####
data_exprs <- set %>%
    exprs() %>%
    as_tibble(rownames = "AffyID") %>%
    dplyr::filter(AffyID %in% probes) %>%
    pivot_longer(-AffyID, names_to = "CEL") %>%
    pivot_wider(names_from = AffyID, values_from = value)



# JOIN THE EXPRESSION DATA WITH PHENOTYPE DATA ####
data <- data_phenotype %>%
    left_join(data_exprs, by = "CEL")



# FILTER DATA BY TREATMENT ARM ####
data_felz <- data %>%
    dplyr::filter(Felzartamab == "Felzartamab") %>%
    arrange(Patient) %>%
    dplyr::select(-CEL) %>%
    nest(.by = c("Patient", "Felzartamab", "Followup", "Felzartamab_Followup"))

data_felz$data[[1]]  %>% flatten_dbl()


# CALCULATE MEAN EXPRESSION ACROSS ALL SELECTED PROBES ####
data_felz_means <- data_felz %>%
    mutate(
        mean = map_dbl(
            data,
            function(data) {
                data %>%
                    flatten_dbl() %>%
                    mean()
            }
        ),
        geomean = 2^mean
    )


# EXPLORE THE RESULTS ####
data_felz_means %>%
    dplyr::filter(Followup == "Baseline") %>%
    dplyr::select(Patient, geomean)


data_felz_means %>%
    dplyr::filter(Followup == "Week24") %>%
    dplyr::select(Patient, geomean)


data_felz_means %>%
    dplyr::filter(Followup == "Week52") %>%
    dplyr::select(Patient, geomean)
