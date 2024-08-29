# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(ggbeeswarm) # install.packages("ggbeeswarm")
library(ggpubr) # install.packages("ggpubr")
library(patchwork) # install.packages("patchwork")
library(ggrepel) # install.packages("ggrepel")
library(gghalves) # pak::pak("erocoar/gghalves")
library(ggh4x) # install.packages("ggh4x")
library(broom) # install.packages("broom") #for tabular model object transformations
library(car) # install.packages("car") #for variable-II manovas
library(rstatix) # install.packages("rstatix") #for testing manova assumptions
library(vegan) # install.packages("vegan") #for permanova
library(pairwiseAdonis) # library(devtools) #install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(flextable) # install.packages("flextable")
library(officer) # install.packages("officer")
library(emmeans) # install.packages("emmeans") #for post-hoc testing and CLD
library(multcomp) # install.packages("multcomp") #for for CLD
library(ARTool) # install.packages("ARTool") #for non-parametric anova (aligned-rank test) and post-hoc
library(rcompanion) # install.packages("rcompanion") #for non-parametric anova (aligned-rank test) and post-hoc
# Bioconductor libraries
library(Biobase) # BiocManager::install("Biobase")
# Custom operators, functions, and datasets
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
se <- function(x) sd(x) / sqrt(length((x)))
log10zero <- scales::trans_new(
    name = "log10zero",
    transform = function(x) log10(x + 0.001),
    inverse = function(x) 10^x - 0.001
)
# Suppress pesky dplyr reframe info
options(dplyr.reframe.inform = FALSE)
# load data
loadDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/"
load(paste(loadDir, "data_scores_k1208.RData", sep = ""))
# load LM22 key
loadDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/A a New 3BMB report/data/deconvolution/signature matrix/"
load(paste(loadDir, "LM22_key.RData", sep = ""))



# DEFINE THE DATA ####
data <- data_scores_k1208 %>%
    dplyr::select(Center, Patient, Felzartamab, Group, Followup, Felzartamab_Group, Felzartamab_Followup, all_of(LM22_key$label_tsi)) %>%
    dplyr::filter(Group != "FU1b", Patient %nin% c(15, 18)) %>%
    left_join(., summarise(., sample_pairs = n(), .by = Felzartamab_Followup), by = "Felzartamab_Followup") %>%
    relocate(sample_pairs, .after = Felzartamab_Followup) %>%
    arrange(Felzartamab, Patient, Group)



# WRANGLE THE PHENOTYPE DATA ####
data_00 <- data %>%
    expand_grid(category = LM22_key %>% pull(lineage)) %>%
    nest(.by = category) %>%
    mutate(
        features = map(
            category,
            function(category) {
                LM22_key %>%
                    dplyr::filter(lineage == category) %>%
                    pull(label_tsi)
            }
        ),
        data = pmap(
            list(features, data),
            function(features, data) {
                data %>%
                    dplyr::select(
                        Center, Patient, Felzartamab, Group, Followup, Felzartamab_Group, Felzartamab_Followup, sample_pairs,
                        all_of(features)
                    )
            }
        )
    )
data_00$data[[2]]


# WRANGLE THE DATA FOR UNIVARIATE TESTS ####
data_01 <- data_00 %>%
    mutate(
        data_univariate = pmap(
            list(category, data, features),
            function(category, data, features) {
                data %>%
                    pivot_longer(
                        cols = all_of(features),
                        names_to = "variable",
                        values_to = "value"
                    ) %>%
                    distinct(Patient, Felzartamab_Followup, variable, .keep_all = TRUE) %>%
                    nest(.by = variable) %>%
                    left_join(
                        LM22_key %>%
                            dplyr::select(label_tsi, sublabel) %>%
                            dplyr::rename(variable = label_tsi),
                        by = "variable"
                    ) %>%
                    mutate(annotation = category) %>%
                    arrange(annotation, variable)
            }
        )
    ) %>%
    dplyr::select(category, data_univariate) %>%
    unnest(data_univariate) %>%
    mutate(
        score = sublabel,
        n_cat = score %>% unique() %>% length(), .by = category, .after = variable
    )

# data_01$data[[9]] %>%
#     print(n="all")


# UNIVARIATE MEANS AND MEDIANS ####
data_02 <- data_01 %>%
    mutate(
        means = map(
            data,
            function(data) {
                data %>%
                    reframe(
                        means = value %>% mean(),
                        se = value %>% se(),
                        sample_pairs = n(),
                        .by = c(Followup, Felzartamab, Felzartamab_Followup)
                    ) %>%
                    arrange(Felzartamab_Followup)
            }
        ),
        means_delta = map(
            means,
            function(means) {
                means %>%
                    reframe(
                        mean_delta = combn(means, 2, diff) %>% as.numeric(),
                        se_delta = combn(se, 2, mean) %>% as.numeric(),
                        sample_pairs = sample_pairs,
                        .by = c(Felzartamab)
                    ) %>%
                    mutate(
                        Followup_pairwise = rep(c("Baseline - Week24", "Baseline - Week52", "Week24 - Week52"), 2) %>%
                            factor(levels = c("Baseline - Week24", "Baseline - Week52", "Week24 - Week52")),
                        .before = sample_pairs
                    ) %>%
                    mutate(
                        mean_delta_delta = combn(mean_delta, 2, diff) %>% as.numeric(),
                        se_delta_delta = combn(se_delta, 2, mean) %>% as.numeric(),
                        .by = c(Followup_pairwise),
                        .before = sample_pairs
                    )
            }
        )
    )
data_02$data[[1]]
data_02$means[[1]]
data_02$means_delta[[1]]


# UNIVARIATE NONPARAMETRIC TESTS ####
data_03 <- data_02 %>%
    mutate(
        art = map(data, art, formula = value ~ Followup * Felzartamab + (1 | Patient)),
        art_lm = map(
            art,
            function(art) {
                art %>%
                    artlm(term = "Followup:Felzartamab", response = "aligned") %>%
                    summary() %>%
                    pluck("coefficients") %>%
                    as.data.frame() %>%
                    rownames_to_column("Term_lm") %>%
                    dplyr::filter(Term_lm != "(Intercept)")
            }
        ),
        art_aov = map(art, anova, type = "II"),
        art_aov_tidy = map(art_aov, tidy) %>% suppressWarnings(),
        art_con_interaction_default = map(
            art,
            art.con,
            formula = "Followup:Felzartamab", adjust = "fdr", method = "pairwise", interaction = TRUE, response = "art"
        ),
        art_con_interaction_default_tidy = map(art_con_interaction_default, tidy),
        art_con_cld = map(
            art_con_interaction_default_tidy,
            function(art_con_interaction_default_tidy) {
                art_con_interaction_default_tidy %>%
                    as.data.frame() %>%
                    cldList(adj.p.value ~ Followup:Felzartamab, data = .)
            }
        )
    )

names(data_03$art_aov_tidy) <- data_03$variable %>% as.character()
names(data_03$art_con_interaction_default_tidy) <- data_03$variable %>% as.character()
names(data_03$art_con_cld) <- data_03$variable %>% as.character()


# SAVE ART TEST RESULTS ####
felzartamab_ARTanova <- data_03
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/"
save(felzartamab_ARTanova, file = paste(saveDir, "felzartamab_ARTanova_cibersort.RData", sep = ""))


# CREATE FLEXTABLE SUMMARY OF ART MODELS ####
title_art <- paste("Table i. Non-parametric ANOVA (ART) of estimate cell type fractions in biopsies from treated vs untreated patients")
res_art_flextable <- data_03 %>%
    dplyr::select(annotation, n_cat, score, art_aov) %>%
    unnest(everything()) %>%
    dplyr::rename(p.value = `Pr(>F)`) %>%
    mutate(FDR = p.value %>% p.adjust(method = "fdr"), .by = annotation) %>%
    dplyr::select(annotation, score, Term, `F`, p.value, FDR) %>%
    flextable::flextable() %>%
    flextable::add_header_row(values = rep(title_art, ncol_keys(.))) %>%
    flextable::merge_h(part = "header") %>%
    flextable::merge_v(j = 1:2) %>%
    flextable::fontsize(size = 8, part = "all") %>%
    flextable::align(align = "center", part = "all") %>%
    flextable::bg(bg = "white", part = "all") %>%
    flextable::bg(i = ~ FDR < 0.05 & Term == "Followup:Felzartamab", j = 3:6, bg = "#fbff00") %>%
    flextable::colformat_double(j = 2:3, digits = 2) %>%
    flextable::colformat_double(j = 4:ncol_keys(.), digits = 3) %>%
    flextable::border_remove() %>%
    flextable::bold(part = "header") %>%
    flextable::padding(padding = 0, part = "all") %>%
    flextable::border(border = officer::fp_border(), part = "all") %>%
    flextable::autofit()


