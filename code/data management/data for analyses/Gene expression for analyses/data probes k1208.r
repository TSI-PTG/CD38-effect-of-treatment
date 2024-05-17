# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
# Bioconductor libraries
library(Biobase) # BiocManager::install("Biobase")
# Custom operators, functions, and datasets
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
source("C:/R/CD38-effect-of-treatment/code/functions/complex_pivot.R")
# Load new reference set
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0Georg Felz CD38 Vienna/G_Rstuff/data/vienna_1208_6Mar24.RData")
# load reference data
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/data_scores_k1208.RData")


# DEFINE CLEAN EXPRESSION SET ####
set_exprs <- vienna_1208
pData(set_exprs) <- set_exprs %>%
    pData() %>%
    tibble %>% 
    dplyr::select(CEL)


# DEFINE CELS IN PHENOTYPE DATA ####
cels <- data_scores_k1208 %>% pull(CEL)


# MAKE SURE NO MISSING CELS ####
set_exprs$CEL[set_exprs$CEL %nin% cels]


# TRIM SET TO MATCH PHENOTYPE DATA BY CEL ####
set_felzartamab <- set_exprs[, which(set_exprs$CEL %in% c(cels))]
set_felzartamab %>% pData()


# ALIGN PHENOTYPE DATA TO EXPRESSION SET ####
set_pdata <- set_felzartamab %>%
    pData() %>%
    left_join(data_scores_k1208, by = "CEL")


# ALIGN PHENOTYPE DATA TO EXPRESSION SET ####
pData(set_felzartamab) <- set_pdata


# SAVE THE DATA ####
data_expressionset_k1208 <- set_felzartamab
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/"
save(data_expressionset_k1208, file = paste(saveDir, "data_expressionset_k1208.RData", sep = ""))
