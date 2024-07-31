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
    mutate(
        data = map(data, pull, AffyID),
        n = map_dbl(data, length)
    ) %>%
    dplyr::filter(n > 1)


# ISOLATE 'NEW' INJURY MARKER GENES BY CELL STATE ####
injury_markers_new_pct40 <- genes_injury_markers %>%
    unnest(data) %>%
    drop_na(AffyID) %>%
    mutate(diff.pct = pct.1 - pct.2) %>% 
    dplyr::filter(
        cluster %>% str_detect("New"),
        diff.pct > 0.4
        ) %>%
    dplyr::select(cluster, AffyID) %>%
    mutate(cluster = cluster %>% str_replace_all("-", "_")) %>%
    nest(.by = cluster) %>%
    mutate(
        data = map(data, pull, AffyID),
        n = map_dbl(data, length)
    ) %>%
    dplyr::filter(n > 1)


# ISOLATE 'NEW' INJURY MARKER GENES BY CELL STATE ####
injury_markers_new_pct30 <- genes_injury_markers %>%
    unnest(data) %>%
    drop_na(AffyID) %>%
    mutate(diff.pct = pct.1 - pct.2) %>% 
    dplyr::filter(
        cluster %>% str_detect("New"),
        diff.pct > 0.3
        ) %>%
    dplyr::select(cluster, AffyID) %>%
    mutate(cluster = cluster %>% str_replace_all("-", "_")) %>%
    nest(.by = cluster) %>%
    mutate(
        data = map(data, pull, AffyID),
        n = map_dbl(data, length)
    ) %>%
    dplyr::filter(n > 1)


# ISOLATE 'NEW' INJURY MARKER GENES BY CELL STATE ####
injury_markers_new_pct20 <- genes_injury_markers %>%
    unnest(data) %>%
    drop_na(AffyID) %>%
    mutate(diff.pct = pct.1 - pct.2) %>%
    dplyr::filter(
        cluster %>% str_detect("New"),
        diff.pct > 0.2
    ) %>%
    dplyr::select(cluster, AffyID) %>%
    mutate(cluster = cluster %>% str_replace_all("-", "_")) %>%
    nest(.by = cluster) %>%
    mutate(
        data = map(data, pull, AffyID),
        n = map_dbl(data, length)
    ) %>%
    dplyr::filter(n > 1)


# TRANSFORM ISOLATED INJURY MARKERS TO A PBTLIST ####
PBTlist219_AKIinduced <- injury_markers_new$data %>% as.list()
names(PBTlist219_AKIinduced) <- injury_markers_new$cluster

PBTlist219_AKIinduced_diffpct40 <- injury_markers_new_pct40$data %>% as.list()
names(PBTlist219_AKIinduced_diffpct40) <- injury_markers_new_pct40$cluster

PBTlist219_AKIinduced_diffpct30 <- injury_markers_new_pct30$data %>% as.list()
names(PBTlist219_AKIinduced_diffpct30) <- injury_markers_new_pct30$cluster

PBTlist219_AKIinduced_diffpct20 <- injury_markers_new_pct20$data %>% as.list()
names(PBTlist219_AKIinduced_diffpct20) <- injury_markers_new_pct20$cluster


# SAVE THE PBT LIST ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/"
# save(PBTlist219_AKIinduced, file = paste(saveDir, "PBTlist219_AKIinduced_26Jul2024_PTG.RData", sep = ""))
save(PBTlist219_AKIinduced_diffpct40, file = paste(saveDir, "PBTlist219_AKIinduced_diffpct40_31Jul2024_PTG.RData", sep = ""))
save(PBTlist219_AKIinduced_diffpct30, file = paste(saveDir, "PBTlist219_AKIinduced_diffpct30_31Jul2024_PTG.RData", sep = ""))
save(PBTlist219_AKIinduced_diffpct20, file = paste(saveDir, "PBTlist219_AKIinduced_diffpct20_31Jul2024_PTG.RData", sep = ""))

