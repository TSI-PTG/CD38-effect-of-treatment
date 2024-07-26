# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
# Custom operators, functions, and datasets
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
# load gene lists
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/Hinze_injury_markers.RData")
# load pbt lists
load("Z:/DATA/Datalocks/Other data/PBTlist219_14Dec2020_JR.RData")


# ISOLATE 'NEW' INJURY MARKER GENES BY CELL STATE ####
injury_markers_new <- genes_injury_markers %>%
    unnest(data) %>%
    drop_na(AffyID) %>%
    dplyr::filter(cluster %>% str_detect("New")) %>%
    dplyr::select(cluster, AffyID) %>%
    mutate(cluster = cluster %>% str_replace_all("-", "_")) %>%
    nest(.by = cluster) %>%
    mutate(data = map(data, pull, AffyID))


# TRANSFORM ISOLATED INJURY MARKERS TO A PBTLIST ####
PBTlist219_AKIinduced <- injury_markers_new$data %>% as.list()
names(PBTlist219_AKIinduced) <- injury_markers_new$cluster


# SAVE THE PBT LIST ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/"
save(PBTlist219_AKIinduced, file = paste(saveDir, "PBTlist219_AKIinduced_26Jul2024_PTG.RData", sep = ""))

