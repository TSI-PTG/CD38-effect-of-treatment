# HOUSEKEEPING ####
library(tidyverse) # install.packages("tidyverse")
library(readr) # install.packages("readr")
library(haven) # install.packages("haven")
# Custom operators, functions, and datasets
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
# load data
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0Georg Felz CD38 Vienna/G_Rstuff/data/data_K1208.RData")
# load reference set
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0Georg Felz CD38 Vienna/G_Rstuff/data/vienna_1208_6Mar24.RData")



# WRANGLE THE REFERENCE SET DATA ####
data_reference <- vienna_1208 %>%
    pData() %>%
    tibble() %>%
    dplyr::select(
        Center, STUDY_EVALUATION_ID, Felzartamab_presumed, CEL, Group,
        ABMRpm, ggt0, ptcgt0, NKB, DSAST,
        TCMRt, tgt1, igt1, TCB, QCAT
    ) %>%
    dplyr::rename(Trial_Center = Center, Felzartamab = Felzartamab_presumed)



# JOIN THE SPSS REFERENCE SET DATA ####
data_joined <- data_K1208 %>%
    left_join(data_reference,
        by = c("CEL", "Trial_Center", "STUDY_EVALUATION_ID", "Felzartamab", "Group"),
        suffix = c("_spss", "_reference"),
    ) %>%
    dplyr::select(contains(
        Hmisc::.q(
            ABMRpm, ggt0, ptcgt0, NKB, DSAST,
            TCMRt, tgt1, igt1, TCB, QCAT
        )
    ))



# data_joined %>%
#     print(n="all")
