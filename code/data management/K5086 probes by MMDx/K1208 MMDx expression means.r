# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
# Bioconductor libraries
library(Biobase) # BiocManager::install("Biobase")
library(limma) # BiocManager::install("limma")
# Custom operators and functions
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
# load data
if (!exists("K1208")) {
    load("Z:/DATA/Datalocks/Kidney/K1208_Feb28_2022_JR.RData")
}
# load affymap
load("Z:/DATA/Datalocks/Other data/affymap219_21Oct2019_1306_JR.RData")


# DEFINE SET ####
set00 <- K1208


# WRANGLE THE PHENOTYPE DATA ####
pData(set00) <- set00 %>%
    pData() %>%
    mutate(
        MMDx = PFHMMDxmod %>%
            trimws() %>%
            str_replace("\\s*\\(.*?\\)\\s*", "") %>%
            str_replace("\\bABMR\\b", "FABMR") %>% 
            factor(
                levels = c(
                    "NR",
                    "Mixed",
                    "TCMR",
                    "EABMR",
                    "FABMR",
                    "LABMR"
                )
            )
    )
set00 %>%
    pData() %>%
    dplyr::select(PFHMMDxmod, MMDx)



# TRIM SET TO EXCLUDE MISSING MMDx ####
set <- set00[, set00$MMDx %>% complete.cases()]


# DEFINE DESIGN FOR MEAN EXPRESSION CALCULATIONS ####
MMDx <- set$MMDx
design <- model.matrix(~ 0 + MMDx)


# FIT LIMMMA ###
fit_MMDx <- limma::lmFit(set, design)


# GET MEAN EXPRESSION LIMMMA ###
means <- fit_MMDx %>%
    avearrays() %>%
    data.frame() %>%
    rownames_to_column("AffyID") %>%
    tibble() %>%
    mutate_if(is.numeric, ~ 2^. %>% round(0))


# ADD ANNOTATION TO MEAN DATA ####
means_K1208 <- left_join(
    affymap219 %>% dplyr::select(AffyID, Symb, Gene, PBT),
    means,
    by = "AffyID"
) %>% tibble()


# SAVE THE MEAN EXPRESSION DATA ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/"
save(means_K1208, file = paste(saveDir, "mean_expression_K1208_MMDx.RData", sep = ""))
