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
# load affymap
load("Z:/DATA/Datalocks/Other data/affymap219_21Oct2019_1306_JR.RData")
# load reference data
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/K5086Jul12_2024_PTG.RData")


# DEFINE THE SET ####
set <- K5086


# DEFINE PBT SCORES ####
data_scores <- set %>%
    pData() %>%
    dplyr::select(ABMRact, ABMRifng, ABMRnk, ABMRendo) %>%
    drop_na() %>%
    tibble

# DEFINE GENE EXPRESSION DATA ####
data_probes <- set %>%
    exprs() %>%
    t() %>%
    as_tibble()



# CALCULATE SPEARMAN CORRELATIONS ####
.Machine$double.eps <- 1e-300

SCC_probes <- cor(data_probes, data_scores, use = "p", method = "spearman") %>%
    dplyr::as_tibble(rownames = "AffyID") %>%
    dplyr::rename_at(dplyr::vars(data_scores %>% colnames()), ~ paste(., "SCC", sep = "_"))



# WRANGLE THE SPEARMAN CORRELATION RESULTS ####
data_scc <- SCC_probes %>%
    right_join(
        affymap219 %>%
            dplyr::select(AffyID, Symb, Gene, PBT) %>%
            tibble(),
        .,
        by = "AffyID"
    ) %>%
    pivot_longer(contains("_SCC"), names_to = "score", values_to = "SCC") %>%
    nest(.by = score) %>%
    mutate(
        data = map(
            data,
            function(data) {
                data %>%
                    mutate(p = corPvalueStudent(SCC, 60), .after = SCC) %>%
                    arrange(p, SCC)
            }
        )
    )
data_scc$data[[1]]


# EXPORT THE SCC TO AN .RData FILE ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/"
save(data_scc, file = paste(saveDir, "ABMRactivity_geneset_probe_scc_K5086.RData", sep = ""))



# END ####
