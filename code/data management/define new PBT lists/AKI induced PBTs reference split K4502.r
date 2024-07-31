# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(FactoMineR) # install.packages('FactoMineR')
library(flextable) # install.packages('flextable')
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
# load datalock
if (!exists("K5086")) {
    load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/K5086_29jul2024.RData")
}
# load nephrectomies
load("Z:/Genome-Archive/RefData/KidneyReports2021REDCap/data/Conset.RData") # 4 Controls used for PBT calculations
# load new PBT list
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/PBTlist219_AKIinduced_26Jul2024_PTG.RData")


# DEFINE CORE SET ####
set <- K5086[, which(K5086$Cort > 0.1)]
pData(set) <- set %>%
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



# DEFINE SUBSETS FOR EACH PBT ####
df_set <- tibble(
    pbt = PBTlist219_AKIinduced %>% names(),
    datExpr = map(
        pbt,
        function(pbt) {
            pbt %>% print()
            set[featureNames(set) %in% c(PBTlist219_AKIinduced %>% purrr::pluck(pbt)), ] %>%
                Biobase::exprs() %>%
                t() %>%
                as.data.frame()
        }
    )
)
df_set$datExpr[[1]]



# ADD PBT SCORES TO NEPRECTOMIES SET ####
set_reference <- pbt.generate(Conset, Conset, PBTlist219_AKIinduced)




# FILTER PBTS USING REFERENCE MEAN METHOD ####
df_refsplit <- df_set %>%
    mutate(
        mean_population = map(
            datExpr,
            function(datExpr) {
                datExpr %>%
                    apply(., 2, FUN = mean) %>%
                    enframe(name = "AffyID", value = "mean_population")
            }
        ),
        mean_reference = map(
            pbt,
            function(pbt) {
                set_reference[featureNames(set_reference) %in% c(PBTlist219_AKIinduced %>% purrr::pluck(pbt)), ] %>%
                    exprs() %>%
                    t() %>%
                    as_tibble() %>%
                    apply(., 2, FUN = mean) %>%
                    enframe(name = "AffyID", value = "mean_reference")
            }
        ),
        mean_comparison = map2(
            mean_population, mean_reference,
            function(mean_population, mean_reference) {
                left_join(
                    mean_population,
                    mean_reference,
                    by = "AffyID"
                ) %>%
                    mutate(direction = ifelse(mean_population > mean_reference, "up", "down"))
            }
        ),
        up = map(
            mean_comparison,
            function(mean_comparison) {
                mean_comparison %>%
                    dplyr::filter(direction == "up") %>%
                    pull(AffyID)
            }
        ),
        dn = map(
            mean_comparison,
            function(mean_comparison) {
                mean_comparison %>%
                    dplyr::filter(direction == "down") %>%
                    pull(AffyID)
            }
        )
    )

# df_refsplit$mean_comparison[[1]]




# WRANGLE THE NEW PBTS ####
df_pbt_updn <- df_refsplit %>%
    dplyr::select(pbt, up, dn) %>%
    pivot_longer(-pbt, names_to = "direction", values_to = "pbtlist") %>%
    mutate(pbtmodule = map2_chr(pbt, direction, paste, sep = "_"))
names(df_pbt_updn$pbtlist) <- df_pbt_updn$pbtmodule


# COLLAPSE NEW PBTLIST ####
PBTlist219_AKIinduced_refsplit <- df_pbt_updn$pbtlist
PBT_names <- PBTlist219_AKIinduced_refsplit %>% names()


# ADD PBT SCORES TO REFERENCE SET ####
set1 <- pbt.generate(set, Conset, PBTlist219_AKIinduced_refsplit)


# SUMMARIZE N GENES IN SET ####
n_PBTlist219_AKIinduced_PCsplit <- PBTlist219_AKIinduced_refsplit %>%
    enframe(name = "pbt") %>%
    mutate("n genes in pbt" = map_dbl(value, length)) %>%
    dplyr::select(-value)


# WRANGLE THE DATA TO MAKE SUMMARY TABLES ####
set_pdata <- set1 %>%
    pData() %>%
    dplyr::select(injAA, contains("New")) %>%
    nest(.by = injAA) %>%
    mutate(means = map(data, summarise_all, mean)) %>%
    dplyr::select(-data) %>%
    unnest(means) %>%
    pivot_longer(cols = -injAA, names_to = "pbt") %>%
    pivot_wider(names_from = injAA, values_from = value) %>%
    mutate_if(is.numeric, ~ formatC(., digits = 3, format = "f")) %>%
    left_join(n_PBTlist219_AKIinduced_PCsplit, by = "pbt") %>%
    mutate(
        direction = pbt %>% str_extract("up|dn"),
        pbt = pbt %>% str_remove("_up|_dn")
    ) %>%
    dplyr::select(pbt, direction, "n genes in pbt", MildCKD, CKDAKI, AKI1, AKI2, Normal)


# SAVE THE PBTLIST ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/"
save(PBTlist219_AKIinduced_refsplit, file = paste(saveDir, "PBTlist219_AKIinduced_refsplit.RData", sep = ""))


# TABULATE PBT SCORES ####
set_pdata %>%
    flextable::flextable() %>%
    flextable::add_header_row(
        top = TRUE,
        values = rep("Mean PBT score of AKI markers filtered by mean expression in Nephrectomies", flextable::ncol_keys(.))
    ) %>%
    flextable::merge_h(part = "header") %>%
    flextable::merge_v(j = 1, part = "body") %>%
    flextable::border_remove() %>%
    flextable::border(part = "all", border = officer::fp_border()) %>%
    flextable::border(part = "footer", border.left = officer::fp_border(), border.right = officer::fp_border()) %>%
    flextable::border(i = 1, part = "footer", border.bottom = officer::fp_border()) %>%
    flextable::align(align = "center", part = "all") %>%
    flextable::font(fontname = "Arial", part = "all") %>%
    flextable::fontsize(size = 12, part = "all") %>%
    flextable::fontsize(i = 1, size = 15, part = "header") %>%
    flextable::bold(part = "header") %>%
    flextable::bg(bg = "white", part = "all") %>%
    flextable::padding(padding = 0, part = "all") %>%
    flextable::width(., width = dim(.)$widths * 30 / (flextable::flextable_dim(.)$widths), unit = "cm") %>%
    print(preview = "pptx")
