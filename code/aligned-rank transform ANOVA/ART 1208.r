# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(broom) # install.packages("broom") #for tabular model object transformations
library(rstatix) # install.packages("rstatix") #for testing manova assumptions
library(flextable) # install.packages("flextable") #for table outputs
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
# load reference data
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/data_felzartamab_k1208.RData")



# DEFINE SEED ####
seed <- 42


# DEFINE CATEGORIES FOR FEATURES ####
# Rejectionrelated <- c("GRIT3", "Rej-RAT", "RejAA_NR")
vars_cfDNA <- c("cfDNA")
vars_abmr <- c("ABMRpm", "ggt0", "ptcgt0", "DSAST", "NKB") # "cggt0", "RejAA_EABMR", "RejAA_FABMR", "RejAA_LABMR")
vars_tcmr <- c("TCMRt", "tgt1", "igt1", "QCAT", "TCB") # , "TCMR-RAT", )
vars_injury <- c("IRRAT30", "IRITD3", "IRITD5")
# Endothelium <- c("ENDAT")
# Parenchyma <- c("KT1", "KT2")
# Macrophage <- c("AMAT1", "QCMAT")
# Injuryrecent <- c("FICOL", "IRRAT30", "IRITD3", "IRITD5")
# Injurylate <- c("IGT", "MCAT", "BAT", "cigt1", "ctgt1")


# DEFINE VARIABLES TO ASSESS ####
# vars <- c("cfDNA", "ABMRpm", "ggt0", "ptcgt0", "DSAST", "NKB", "TCMRt", "tgt1", "igt1", "QCAT", "TCB")
vars <- c(vars_cfDNA, vars_abmr, vars_tcmr, vars_injury)


# DEFINE THE data ####
data <- data_felzartamab_k1208 %>%
    dplyr::select(Center, Patient, Felzartamab, Group, Followup, Felzartamab_Group, Felzartamab_Followup, all_of(vars)) %>%
    dplyr::filter(Group != "FU1b", Patient %nin% c(15, 18)) %>%
    left_join(., summarise(., sample_pairs = n(), .by = Felzartamab_Followup), by = "Felzartamab_Followup") %>%
    relocate(sample_pairs, .after = Felzartamab_Followup) %>%
    arrange(Felzartamab, Patient, Group)



# WRANGLE THE PHENOTYPE DATA ####
data_00 <- data %>%
    expand_grid(category = c("cfDNA", "ABMR", "TCMR", "injury")) %>%
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
                } else if (category == "injury") {
                    vars_injury
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
                            variable %in% vars_injury ~ "injury-related",
                            TRUE ~ " "
                        ) %>%
                            factor(
                                levels = c(
                                    "cfDNA",
                                    "ABMR-related",
                                    "TCMR-related",
                                    "injury-related"
                                )
                            ),
                        score = case_when(
                            variable == "cfDNA" ~ "Donor-derived cell-free DNA (dd-cfDNA, cp/mL)",
                            variable == "TCMRt" ~ "TCMR classifier (TCMRProb)",
                            variable == "TCMRt" ~ "TCMR classifier (TCMRProb)",
                            variable == "TCB" ~ "T-cell burden (TCB)",
                            variable == "TCMR-RAT" ~ "TCMR-associated RATs (TCMR-RAT)",
                            variable == "tgt1" ~ "Tubulitis classifier (t>1Prob)",
                            variable == "igt1" ~ "Interstitial infiltrate classifier (i>1Prob)",
                            # variable == "QCAT" ~ "Cytotoxic T cell-associated transcripts (QCAT)",
                            variable == "QCAT" ~ "Cytotoxic T cell-associated (QCAT)",
                            variable == "Rej-RAT" ~ "Rejection-associated transcripts (Rej-RAT)",
                            variable == "GRIT3" ~ "Interferon gamma-inducible transcripts (GRIT3)",
                            variable == "ABMR-RAT" ~ "ABMR-associated RATs (ABMR-RAT)",
                            variable == "ABMRpm" ~ "ABMR classifier (ABMRProb)",
                            variable == "DSAST" ~ "DSA-selective transcripts (DSAST)",
                            variable == "NKB" ~ "NK cell burden (NKB)",
                            variable == "ggt0" ~ "Glomerulitis classifier (g>0Prob)",
                            variable == "cggt0" ~ "Glomerular double contours classifier (cg>0Prob)",
                            variable == "ptcgt0" ~ "Peritubular capillaritis classifier (ptc>0Prob)",
                            variable == "AMAT1" ~ "Alternatively activated macrophage (AMAT1)",
                            variable == "QCMAT" ~ "Constitutive macrophage (QCMAT)",
                            variable == "FICOL" ~ "Fibrillar collagen (FICOL)",
                            variable == "DAMP" ~ "Damage-associated molecular pattern transcripts (DAMP)",
                            variable == "cIRIT" ~ "Cardiac injury and repair–induced transcripts (cIRIT)",
                            variable == "IRITD3" ~ "Injury-repair induced, day 3 (IRITD3)",
                            variable == "IRITD5" ~ "Injury-repair induced, day 5 (IRITD5)",
                            variable == "IRRAT30" ~ "Injury-repair associated (IRRAT30)",
                            variable == "cigt1" ~ "Fibrosis classifier (ci>1Prob)",
                            variable == "ctgt1" ~ "Atrophy classifier (ct>1Prob)",
                            variable == "txbxCorrCLAGdn_tbb" ~ "Transcripts downregulated in CLAD",
                            variable == "txbxCorrCLAGup_tbb" ~ "Transcripts upregulated in CLAD",
                            variable == "SFT" ~ "Surfactant-associated transcripts (SFT)",
                            variable == "ENDAT" ~ "Endothelial cell-associated transcripts (ENDAT)",
                            variable == "eDSAST" ~ "Endothelium-expressed DSA-selective transcripts (eDSAST)",
                            variable == "IGT" ~ "Immunoglobulin transcripts (IGT)",
                            variable == "BAT" ~ "B cell–associated transcripts (BAT)",
                            variable == "MCAT" ~ "Mast cell-associated transcripts (MCAT)",
                            variable == "KT1" ~ "Kidney parenchymal transcripts (KT1)",
                            variable == "KT2" ~ "Kindey parenchymal transcripts - no solute carriers (KT2)",
                            variable == "RejAA_NR" ~ "Archetypal No Rejection score (NR)",
                            variable == "RejAA_EABMR" ~ "Archetypal Early ABMR score (EABMR)",
                            variable == "RejAA_FABMR" ~ "Archetypal Full ABMR score (FABMR)",
                            variable == "RejAA_LABMR" ~ "Archetypal Late ABMR score (LABMR)",
                        ), .before = 1
                    ) %>%
                    arrange(annotation, variable)
            }
        )
    ) %>%
    dplyr::select(category, data_univariate) %>%
    unnest(data_univariate) %>%
    mutate(n_cat = score %>% unique() %>% length(), .by = category, .after = variable)


# UNIVARIATE MEANS AND MEDIANS ####
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
data_02$medians[[1]]
data_02$medians_delta[[1]]


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

# data_03$art_con_cld
# data_03$medians[[1]]
# data_03$medians_delta[[1]]
# data_03$art_aov_tidy[[1]]
# data_03$art_con_interaction_default[[1]]
# data_03$art_con_interaction_default_tidy[[1]]
# data_03$art_con_cld[[1]]


# SAVE ART TEST RESULTS ####
felzartamab_ARTanova <- data_03
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/"
save(felzartamab_ARTanova, file = paste(saveDir, "felzartamab_ARTanova.RData", sep = ""))



# CREATE FLEXTABLE SUMMARY OF ART MODELS ####
title_art <- paste("Table i. Non-parametric ANOVA (ART) of molecular scores in biopsies from treated vs untreated patients")
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
    flextable::border(border = fp_border(), part = "all") %>%
    flextable::autofit()



# END ####
