# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(broom) # install.packages("broom") #for tabular model object transformations
library(rstatix) # install.packages("rstatix") #for testing manova assumptions
library(flextable) # install.packages("flextable")
library(officer) # install.packages("officer")
library(emmeans) # install.packages("emmeans") #for post-hoc testing and CLD
library(multcomp) # install.packages("multcomp") #for for CLD
library(ARTool) # install.packages("ARTool") #for non-parametric anova (aligned-rank test) and post-hoc
library(rcompanion) # install.packages("rcompanion") #for non-parametric anova (aligned-rank test) and post-hoc
# Bioconductor libraries
library(Biobase) # BiocManager::install("Biobase")
# Custom operators and functions
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
# Suppress pesky dplyr reframe info
options(dplyr.reframe.inform = FALSE)
# load reference set
load("data/cd38.RData")
# load affymap
load("data/affymap219.RData")


# DEFINE SEED ####
seed <- 42


# DEFINE CATEGORIES FOR FEATURES ####
vars_cfDNA <- c("cfDNA_cpml")
vars_abmr <- c("ABMRpm", "ggt0", "ptcgt0", "DSAST", "NKB")
vars_tcmr <- c("TCMRt", "tgt1", "igt1", "QCAT", "TCB")
vars_injury <- c("IRRAT30", "IRITD3", "IRITD5", "cigt1", "ctgt1")
vars_parenchyma <- c("KT1", "KT2")
vars_macrophage <- c("AMAT1", "QCMAT")


# DEFINE VARIABLES TO ASSESS ####
vars <- c(vars_cfDNA, vars_abmr, vars_tcmr, vars_macrophage, vars_injury, vars_parenchyma)


# DEFINE THE DATA ####
data <- set %>%
    pData()  %>% 
    dplyr::select(Center, Patient, Felzartamab, Group, Followup, Felzartamab_Group, Felzartamab_Followup, all_of(vars)) %>%
    dplyr::filter(Patient %nin% c(15, 18)) %>%
    left_join(., summarise(., sample_pairs = n(), .by = Felzartamab_Followup), by = "Felzartamab_Followup") %>%
    relocate(sample_pairs, .after = Felzartamab_Followup) %>%
    arrange(Felzartamab, Patient, Group)



# WRANGLE THE PHENOTYPE DATA ####
data_00 <- data %>%
    expand_grid(category = c("cfDNA", "ABMR", "TCMR", "macrophage", "injury", "parenchyma")) %>%
    nest(.by = category) %>%
    mutate(
        features = map(
            category,
            function(category) {
                if (category == "cfDNA") {
                    vars_cfDNA
                } else if (category == "ABMR") {
                    vars_abmr
                } else if (category == "TCMR") {
                    vars_tcmr
                } else if (category == "macrophage") {
                    vars_macrophage
                } else if (category == "injury") {
                    vars_injury
                } else if (category == "parenchyma") {
                    vars_parenchyma
                } 
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


# WRANGLE THE DATA FOR UNIVARIATE TESTS ####
data_01 <- data_00 %>%
    mutate(
        data_univariate = pmap(
            list(data, features),
            function(data, features) {
                data %>%
                    pivot_longer(
                        cols = all_of(features),
                        names_to = "variable",
                        values_to = "value"
                    ) %>%
                    nest(.by = variable) %>%
                    mutate(
                        variable = variable %>%
                            factor(
                                levels = vars
                            ),
                        annotation = case_when(
                            variable %in% vars_cfDNA ~ "cfDNA",
                            variable %in% vars_tcmr ~ "TCMR-related",
                            variable %in% vars_abmr ~ "ABMR-related",
                            variable %in% vars_macrophage ~ "macrophage-related",
                            variable %in% vars_injury ~ "injury-related",
                            variable %in% vars_parenchyma ~ "parenchyma-related",
                            TRUE ~ " "
                        ) %>%
                            factor(
                                levels = c(
                                    "cfDNA",
                                    "ABMR-related",
                                    "TCMR-related",
                                    "macrophage-related",
                                    "injury-related",
                                    "parenchyma-related",
                                    "archetypes",
                                    "rejection PC",
                                    "injury PC"
                                )
                            ),
                        score = case_when(
                            variable == "cfDNA_cpml" ~ "Donor-derived cell-free DNA (dd-cfDNA, cp/mL)",
                            variable == "TCMRt" ~ "TCMR classifier (TCMRProb)",
                            variable == "TCB" ~ "T cell burden (TCB)",
                            variable == "tgt1" ~ "Tubulitis classifier (t>1Prob)",
                            variable == "igt1" ~ "Interstitial infiltrate classifier (i>1Prob)",
                            variable == "QCAT" ~ "Cytotoxic T cell-associated (QCAT)",
                            variable == "ABMRpm" ~ "ABMR classifier (ABMRProb)",
                            variable == "DSAST" ~ "DSA-selective transcripts (DSAST)",
                            variable == "NKB" ~ "NK cell burden (NKB)",
                            variable == "ggt0" ~ "Glomerulitis classifier (g>0Prob)",
                            variable == "ptcgt0" ~ "Peritubular capillaritis classifier (ptc>0Prob)",
                            variable == "AMAT1" ~ "Alternatively activated macrophage (AMAT1)",
                            variable == "QCMAT" ~ "Constitutive macrophage (QCMAT)",
                            variable == "IRITD3" ~ "Injury-repair induced, day 3 (IRITD3)",
                            variable == "IRITD5" ~ "Injury-repair induced, day 5 (IRITD5)",
                            variable == "IRRAT30" ~ "Injury-repair associated (IRRAT30)",
                            variable == "cigt1" ~ "Fibrosis classifier (ci>1Prob)",
                            variable == "ctgt1" ~ "Atrophy classifier (ct>1Prob)",
                            variable == "KT1" ~ "Kidney parenchymal (KT1)",
                            variable == "KT2" ~ "Kidney parenchymal - no solute carriers (KT2)"
                        ), .before = 1
                    ) %>%
                    arrange(annotation, variable)
            }
        )
    ) %>%
    dplyr::select(category, data_univariate) %>%
    unnest(data_univariate) %>%
    mutate(n_cat = score %>% unique() %>% length(), .by = category, .after = variable)


# UNIVARIATE MEDIANS ####
data_02 <- data_01 %>%
    mutate(
        medians = map(
            data,
            function(data) {
                data %>%
                    reframe(
                        median = value %>% median(),
                        IQR = value %>% IQR(),
                        sample_pairs = n(),
                        .by = c(Followup, Felzartamab, Felzartamab_Followup)
                    ) %>%
                    arrange(Felzartamab_Followup)
            }
        ),
        medians_delta = map(
            medians,
            function(medians) {
                medians %>%
                    reframe(
                        median_delta = combn(median, 2, diff) %>% as.numeric(),
                        IQR_delta = combn(IQR, 2, mean) %>% as.numeric(),
                        sample_pairs = sample_pairs,
                        .by = c(Felzartamab)
                    ) %>%
                    mutate(
                        Followup_pairwise = rep(c("Baseline - Week24", "Baseline - Week52", "Week24 - Week52"), 2) %>%
                            factor(levels = c("Baseline - Week24", "Baseline - Week52", "Week24 - Week52")),
                        .before = sample_pairs
                    ) %>%
                    mutate(
                        median_delta_delta = combn(median_delta, 2, diff) %>% as.numeric(),
                        IQR_delta_delta = combn(IQR_delta, 2, mean) %>% as.numeric(),
                        .by = c(Followup_pairwise),
                        .before = sample_pairs
                    )
            }
        )
    )
data_02$data[[1]]
data_02 %>%
    dplyr::filter(variable == "ABMRpm") %>%
    pull(medians)
data_02 %>%
    dplyr::filter(variable == "ABMRpm") %>%
    pull(medians_delta)


# UNIVARIATE NONPARAMETRIC TESTS ####
artANOVA <- data_02 %>%
    mutate(
        art = map(data, art, formula = value ~ Followup * Felzartamab + (1 | Patient)),
        art_aov = map(art, anova, type = "II"),
        art_aov_tidy = map(art_aov, tidy) %>% suppressWarnings(),
        art_aov_contrast = map(
            art,
            art.con,
            formula = "Followup:Felzartamab", adjust = "fdr", method = "pairwise", interaction = TRUE, response = "art"
        ),
        art_aov_contrast_tidy = map(art_aov_contrast, tidy),
        art_aov_contrast_cld = map(
            art_aov_contrast_tidy,
            function(art_aov_contrast_tidy) {
                art_aov_contrast_tidy %>%
                    as.data.frame() %>%
                    cldList(adj.p.value ~ Followup:Felzartamab, data = .)
            }
        )
    )

names(artANOVA$art_aov_tidy) <- artANOVA$variable %>% as.character()
names(artANOVA$art_aov_contrast_tidy) <- artANOVA$variable %>% as.character()
names(artANOVA$art_aov_contrast_cld) <- artANOVA$variable %>% as.character()


# SAVE ART TEST RESULTS ####
saveDir <- "results/"
save(artANOVA, file = paste(saveDir, "artANOVA.RData", sep = ""))


# CREATE FLEXTABLE SUMMARY OF ART MODELS ####
title_art <- paste("Table i. Non-parametric ANOVA (ART) of molecular scores in biopsies from treated vs untreated patients")
artANOVA %>%
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



# END ####
