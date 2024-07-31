# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
# Custom operators, functions, and datasets
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
pbt.generate <- function(data_CEL, reference_set, pbtlist) {
    for (i in 1:length(names(pbtlist))) {
        Biobase::pData(data_CEL)[names(pbtlist)[i]] <- NA
        conavg <- apply(Biobase::exprs(reference_set)[pbtlist[[i]], ], 1, mean)
        smallset <- data_CEL[pbtlist[[i]], ]
        dog <- sweep(Biobase::exprs(smallset), 1, conavg)
        dog <- apply(dog, 2, mean)
        sel <- which(names(Biobase::pData(data_CEL)) == names(pbtlist)[i])
        Biobase::pData(data_CEL)[, sel] <- dog
    }
    return(data_CEL)
}
# load data
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/K5086_29jul2024.RData")
# load nephrectomies
load("Z:/Genome-Archive/RefData/KidneyReports2021REDCap/data/Conset.RData") # 4 Controls used for PBT calculations
# load new PBT list
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/PBTlist219_AKIinduced_26Jul2024_PTG.RData")
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/PBTlist219_AKIinduced_diffpct40_31Jul2024_PTG.RData")
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/PBTlist219_AKIinduced_diffpct30_31Jul2024_PTG.RData")
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/PBTlist219_AKIinduced_diffpct20_31Jul2024_PTG.RData")


# WRANGLE REFERENCE SET ####
K4502 <- K5086[, which(K5086$Cort > 0.1)]
pData(K4502) <- K4502 %>%
    pData() %>%
    tibble() %>%
    mutate(
        injAA = case_when(
            InjAA5Clust == 1 ~ "MildCKD",
            InjAA5Clust == 2 ~ "CKDAKI",
            InjAA5Clust == 3 ~ "AKI1",
            InjAA5Clust == 4 ~ "AKI2",
            InjAA5Clust == 5 ~ "Normal"
        )
    ) %>%
    dplyr::select(injAA)


# GENERATE PBT SCORES ####
set <- pbt.generate(K4502, Conset, PBTlist219_AKIinduced)
set20 <- pbt.generate(K4502, Conset, PBTlist219_AKIinduced_diffpct20)
set30 <- pbt.generate(K4502, Conset, PBTlist219_AKIinduced_diffpct30)
set40 <- pbt.generate(K4502, Conset, PBTlist219_AKIinduced_diffpct40)


# WRANGLE THE PHENOTYPE DATA ####
set_pdata <- set %>%
    pData() %>%
    dplyr::select(injAA, contains("New")) %>%
    nest(.by = injAA) %>%
    mutate(means = map(data, summarise_all, mean)) %>%
    dplyr::select(-data) %>%
    unnest(means) %>%
    pivot_longer(cols = -injAA, names_to = "pbt") %>%
    pivot_wider(names_from = injAA, values_from = value) %>%
    dplyr::select(pbt, MildCKD, CKDAKI, AKI1, AKI2, Normal)

set20_pdata <- set20 %>%
    pData() %>%
    dplyr::select(injAA, contains("New")) %>%
    nest(.by = injAA) %>%
    mutate(means = map(data, summarise_all, mean)) %>%
    dplyr::select(-data) %>%
    unnest(means) %>%
    pivot_longer(cols = -injAA, names_to = "pbt") %>%
    pivot_wider(names_from = injAA, values_from = value) %>%
    dplyr::select(pbt, MildCKD, CKDAKI, AKI1, AKI2, Normal)

set30_pdata <- set30 %>%
    pData() %>%
    dplyr::select(injAA, contains("New")) %>%
    nest(.by = injAA) %>%
    mutate(means = map(data, summarise_all, mean)) %>%
    dplyr::select(-data) %>%
    unnest(means) %>%
    pivot_longer(cols = -injAA, names_to = "pbt") %>%
    pivot_wider(names_from = injAA, values_from = value) %>%
    dplyr::select(pbt, MildCKD, CKDAKI, AKI1, AKI2, Normal)

set40_pdata <- set40 %>%
    pData() %>%
    dplyr::select(injAA, contains("New")) %>%
    nest(.by = injAA) %>%
    mutate(means = map(data, summarise_all, mean)) %>%
    dplyr::select(-data) %>%
    unnest(means) %>%
    pivot_longer(cols = -injAA, names_to = "pbt") %>%
    pivot_wider(names_from = injAA, values_from = value) %>%
    dplyr::select(pbt, MildCKD, CKDAKI, AKI1, AKI2, Normal)


# SUMMARISE PBT SCORES ####
