# HOUSEKEEPING ####
# CRAN packages
library(tidyverse) # install.packages("tidyverse")
library(quantreg) # install.packages("quantreg")
library(broom) # install.packages("broom")
# Custom operators, functions, and datasets
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
# Custom operators, functions, and datasets
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
se <- function(x) sd(x) / sqrt(length((x)))
corPvalueStudent <- function(cor, nSamples) {
    T <- sqrt(nSamples - 2) * cor / sqrt(1 - cor^2)
    2 * pt(abs(T), nSamples - 2, lower.tail = FALSE)
}
# load reference data
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/data_expressionset_k1208.RData")


# DEFINE THE SET ####
set <- data_expressionset_k1208[, data_expressionset_k1208$Patient %nin% c(15, 18)]


# DEFINE PBT SCORES ####
data_scores <- set %>%
    pData() %>%
    dplyr::select(CEL, ABMRact, ABMRifng, ABMRnk, ABMRendo) %>%
    drop_na


# DEFINE GENE EXPRESSION DATA ####
data_probes <- set %>%
    exprs() %>%
    t() %>%
    as_tibble(rownames = "CEL")










# END ####
