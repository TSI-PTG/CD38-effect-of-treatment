# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
# Bioconductor libraries
library(Biobase) # BiocManager::install("Biobase")
# Custom operators, functions, and datasets
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
# load affymap
load("Z:/DATA/Datalocks/Other data/affymap219_21Oct2019_1306_JR.RData")
# load reference set
load("Z:/DATA/Datalocks/Kidney/K5086Mar13_2024_JR.RData")


# WRANGLE THE EXPRESSION SET DATA ####
data_exprs <- K5086 %>%
    exprs() %>%
    as_tibble(rownames = "AffyID") %>%
    left_join(affymap219 %>% dplyr::select(AffyID, Symb), by = "AffyID") %>%
    relocate(Symb, .after = AffyID)


# CALCULATE MEAN EXPRESSION ####
mean_exprs_5086 <- data_exprs %>%
    mutate(
        mean_expression = data_exprs %>% dplyr::select(-AffyID, -Symb) %>% rowMeans(),
        .after = Symb
    ) %>%
    dplyr::select(AffyID:mean_expression)


# EXPORT THE DATA AS .RData FILE ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/"
save(mean_exprs_5086, file = paste(saveDir, "mean_expression_by_probe_5086.RData", sep = ""))
