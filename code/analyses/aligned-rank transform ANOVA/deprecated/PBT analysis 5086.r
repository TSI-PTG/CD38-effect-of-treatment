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
library(car) # install.packages("car") #for type-II manovas
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
# source plot function
source("C:/R/CD38-effect-of-treatment/code/functions/plot.gg_violin_interaction.r")
# laod reference set
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0Georg Felz CD38 Vienna/G_Rstuff/data/vienna_5086_6Mar24.RData")


# DEFINE SEED ####
seed <- 42


# DEFINE CATEGORIES FOR FEATURES ####
Rejectionrelated <- c("GRIT3", "Rej-RAT", "RejAA_NR")
ABMRrelated <- c("DSAST", "NKB", "ABMRpm", "ggt0", "ptcgt0") # "cggt0", "RejAA_EABMR", "RejAA_FABMR", "RejAA_LABMR")
TCMRrelated <- c("QCAT", "TCB", "TCMRt", "tgt1", "igt1") # , "TCMR-RAT", )
Endothelium <- c("ENDAT")
Parenchyma <- c("KT1", "KT2")
Macrophage <- c("AMAT1", "QCMAT")
Injuryrecent <- c("FICOL", "IRRAT30", "IRITD3", "IRITD5")
Injurylate <- c("IGT", "MCAT", "BAT", "cigt1", "ctgt1")


# DEFINE FEATURES ####
features <- c(ABMRrelated, TCMRrelated, Macrophage, Injurylate)


# PATIENT SUMMARIES ####
patient_summary <- vienna_5086 %>%
    pData() %>%
    dplyr::select(Center, STUDY_EVALUATION_ID, Felzartamab_presumed, CEL, Group) %>%
    tibble() %>%
    mutate(
        Followup = case_when(Group == "Index" ~ "Index\n(Baseline)", Group == "FU1" ~ "FU1\n(week24)", Group == "FU2" ~ "FU2\n(week52)"),
        Felzartamab_presumed = Felzartamab_presumed %>% factor(),
        STUDY_EVALUATION_ID = STUDY_EVALUATION_ID %>% factor()
    ) %>%
    dplyr::filter(CEL %nin% c("FBN003_NBN010_B2_(PrimeView).CEL", "FVI022_FVI022_B2_(PrimeView).CEL")) %>%
    dplyr::select(-Group) %>%
    pivot_wider(names_from = Followup, values_from = CEL) %>%
    arrange(STUDY_EVALUATION_ID, Felzartamab_presumed)

patient_summary %>%
    mutate(across(c("Index\n(Baseline)", "FU1\n(week24)", "FU2\n(week52)"), ~ ifelse(. %>% is.na(), NA, "X"))) %>%
    flextable::flextable() %>%
    flextable::add_header_row(values = rep("Felzartamab study population", ncol_keys(.))) %>%
    flextable::merge_h(part = "header") %>%
    flextable::merge_v(j = 1:2) %>%
    flextable::fontsize(size = 8, part = "all") %>%
    flextable::align(align = "center", part = "all") %>%
    flextable::bg(bg = "white", part = "all") %>%
    flextable::colformat_double(j = 2:3, digits = 2) %>%
    flextable::border_remove() %>%
    flextable::bold(part = "header") %>%
    flextable::padding(padding = 0, part = "all") %>%
    flextable::border(border = officer::fp_border(), part = "all") %>%
    flextable::autofit()
# %>%
# print(preview = "pptx")

# table with cel id's
patient_summary %>%
    flextable::flextable() %>%
    flextable::add_header_row(values = rep("Felzartamab study population", ncol_keys(.))) %>%
    flextable::merge_h(part = "header") %>%
    flextable::merge_v(j = 1:2) %>%
    flextable::fontsize(size = 8, part = "all") %>%
    flextable::align(align = "center", part = "all") %>%
    flextable::bg(bg = "white", part = "all") %>%
    flextable::colformat_double(j = 2:3, digits = 2) %>%
    flextable::border_remove() %>%
    flextable::bold(part = "header") %>%
    flextable::padding(padding = 0, part = "all") %>%
    flextable::border(border = officer::fp_border(), part = "all") %>%
    flextable::autofit()

vienna_5086[, vienna_5086$STUDY_EVALUATION_ID %nin% c(15, 18)] %>%
    pData() %>%
    dplyr::select(Felzartamab_presumed, Group) %>%
    table() %>%
    as_tibble() %>%
    pivot_wider(names_from = Group, values_from = n)


# DEFINE THE SET ####
set <- vienna_5086[, vienna_5086$STUDY_EVALUATION_ID %nin% c(15, 18)]
# set <- vienna_5086[, vienna_5086$STUDY_EVALUATION_ID != 15 & vienna_5086$CEL != "FBN003_NBN010_B2_(PrimeView).CEL"]


# WRANGLE THE PHENOTYPE DATA ####
df00 <- set %>%
    pData() %>%
    tibble() %>%
    dplyr::rename(
        Patient = STUDY_EVALUATION_ID,
        Felz = Felzartamab_presumed
    ) %>%
    mutate(
        Patient = Patient %>% factor(),
        Group = Group %>% factor(levels = c("Index", "FU1", "FU2")),
        Felz = Felz %>% factor(labels = c("NoFelz", "Felz")),
        TxBx = TxBx %>% as.numeric(),
        Group_Felz = paste(Group, Felz, sep = ":") %>%
            factor(levels = c("Index_NoFelz", "FU1_NoFelz", "FU2_NoFelz", "Index_Felz", "FU1_Felz", "FU2_Felz"))
    ) %>%
    arrange(Patient, Group) %>%
    expand_grid(category = c("ABMR", "TCMR", "Macrophage", "Injurylate")) %>%
    nest(.by = category) %>%
    mutate(
        features = map(
            category,
            function(category) {
                if (category == "ABMR") {
                    ABMRrelated
                } else if (category == "TCMR") {
                    TCMRrelated
                } else if (category == "Macrophage") {
                    Macrophage
                } else if (category == "Injurylate") {
                    Injurylate
                }
            }
        ),
        data = pmap(
            list(features, data),
            function(features, data) {
                data %>%
                    dplyr::select(
                        CEL, Patient, Center, Group, Felz, Group_Felz, TxBx,
                        all_of(features)
                    )
            }
        )
    )


# PERMANOVA ####
df_permanova <- df00 %>%
    mutate(
        permanova = pmap(
            list(data, features),
            function(data, features) {
                set.seed(seed)
                features <- data %>% dplyr::select(all_of(features))
                Group <- data$Group
                Felz <- data$Felz
                adonis2(
                    features ~ Group * Felz + Patient,
                    data = data,
                    method = "euclidean",
                    by = "margin", # only specify margin if the sample sizes are unequal
                    permutations = 10000
                )
            }
        ),
        permanova_pairwise = pmap(
            list(data, features),
            function(data, features) {
                set.seed(seed)
                features <- data %>% dplyr::select(all_of(features))
                Group_Felz <- data$Group_Felz
                pairwise.adonis(
                    features,
                    Group_Felz,
                    reduce = "Group|Felz",
                    sim.method = "euclidean",
                    p.adjust.m = "fdr",
                    perm = 10000
                )
            }
        )
    )
df_permanova$permanova
# df_permanova$permanova_pairwise[[1]]


# WRANGLE THE DATA FOR UNIVARIATE TESTS ####
df_univariate_00 <- df_permanova %>%
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
                    group_by(variable) %>%
                    nest() %>%
                    tibble() %>%
                    mutate(
                        variable = variable %>%
                            factor(
                                levels = c(
                                    "Rej-RAT", "GRIT3", "RejAA_NR",
                                    "ABMRpm", "ggt0", "cggt0", "ptcgt0", "NKB", "DSAST",
                                    "RejAA_EABMR", "RejAA_FABMR", "RejAA_LABMR",
                                    "TCMRt", "tgt1", "igt1", "TCB", "TCMR-RAT", "QCAT",
                                    "AMAT1", "QCMAT", "BAT",
                                    "FICOL", "IRRAT30", "IRITD3", "IRITD5",
                                    "cigt1", "ctgt1", "IGT", "MCAT",
                                    "KT1", "KT2"
                                )
                            ),
                        annotation = case_when(
                            variable %in% Rejectionrelated ~ "Rejection-related",
                            variable %in% TCMRrelated ~ "TCMR-related",
                            variable %in% ABMRrelated ~ "ABMR-related",
                            variable %in% Endothelium ~ "Endothelium-related",
                            variable %in% Parenchyma ~ "Parenchyma-related",
                            variable %in% Macrophage ~ "Macrophage-related",
                            variable %in% Injurylate ~ "Atrophy-fibrosis-related",
                            variable %in% Injuryrecent ~ "Recent injury-related",
                            TRUE ~ " "
                        ) %>%
                            factor(
                                levels = c(
                                    "Rejection-related", "ABMR-related", "TCMR-related",
                                    "Macrophage-related", "Recent injury-related", "Atrophy-fibrosis-related",
                                    "Parenchyma-related"
                                )
                            ),
                        score = case_when(
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
                            variable == "KT2" ~ "Kidney parenchymal transcripts - no solute carriers (KT2)",
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
    unnest(data_univariate)


# UNIVARIATE MEANS AND MEDIANS ####
df_univariate_01 <- df_univariate_00 %>%
    mutate(
        n_cat = score %>% unique() %>% length(), .by = category, .after = variable
    ) %>%
    mutate(
        medians = map(
            data,
            function(data) {
                medians <- data %>%
                    reframe(
                        median = value %>% median(),
                        IQR = value %>% IQR(),
                        .by = c(Group, Felz)
                    ) %>%
                    mutate(
                        Group_Felz = paste(Felz, Group, sep = ":") %>%
                            factor(levels = c(
                                # "Index_NoFelz", "FU1_NoFelz", "FU2_NoFelz",
                                # "Index_Felz", "FU1_Felz", "FU2_Felz"
                                "NoFelz:Index", "NoFelz:FU1", "NoFelz:FU2",
                                "Felz:Index", "Felz:FU1", "Felz:FU2"
                            ))
                    ) %>%
                    arrange(Group_Felz)
            }
        ),
        medians_delta = map(
            medians,
            function(medians) {
                median <- medians %>%
                    reframe(
                        median_delta = combn(median, 2, diff) %>% as.numeric(),
                        IQR_delta = combn(IQR, 2, mean) %>% as.numeric(),
                        .by = c(Felz)
                    ) %>%
                    mutate(Group_pairwise = rep(c("Index - FU1", "Index - FU2", "FU1 - FU2"), 2)) %>%
                    mutate(
                        median_delta_delta = combn(median_delta, 2, diff) %>% as.numeric(),
                        IQR_delta_delta = combn(IQR_delta, 2, mean) %>% as.numeric(),
                        .by = c(Group_pairwise)
                    )
            }
        )
    )

df_univariate_01$data[[1]]
df_univariate_01$medians[[1]]
df_univariate_01$medians_delta[[1]]


# UNIVARIATE NONPARAMETRIC TESTS ####
df_univariate_02 <- df_univariate_01 %>%
    mutate(
        art = map(data, art, formula = value ~ Group * Felz + (1 | Patient)),
        art_lm = map(
            art,
            function(art) {
                art %>%
                    artlm(term = "Group:Felz", response = "aligned") %>%
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
            formula = "Group:Felz", adjust = "fdr", method = "pairwise", interaction = TRUE, response = "art"
        ),
        art_con_interaction_default_tidy = map(art_con_interaction_default, tidy),
        art_con_cld = map(
            art_con_interaction_default_tidy,
            function(art_con_interaction_default_tidy) {
                art_con_interaction_default_tidy %>%
                    as.data.frame() %>%
                    cldList(adj.p.value ~ Group:Felz, data = .) # %>%
                #         arrange(Group %>%
                #             factor(
                #                 levels = c(
                #                     "Index,NoFelz", "FU1,NoFelz",
                #                     "Index,Felz", "FU1,Felz"
                #                 )
                #             ))
            }
        )
    )

names(df_univariate_02$art_aov_tidy) <- df_univariate_02$variable %>% as.character()
names(df_univariate_02$art_con_interaction_default_tidy) <- df_univariate_02$variable %>% as.character()
names(df_univariate_02$art_con_cld) <- df_univariate_02$variable %>% as.character()

df_univariate_02$art_con_cld

df_univariate_02$medians[[1]]
df_univariate_02$medians_delta[[1]]
df_univariate_02$art_aov_tidy[[1]]
df_univariate_02$art_con_interaction_default[[1]]

df_univariate_02$art_con_interaction_default_tidy[[3]]
df_univariate_02$art_con_cld[[3]]


# CREATE FLEXTABLE OF ART MODELS ####
title_art <- paste("Table i. Non-parametric ANOVA (ART) of molecular scores in biopsies from treated vs untreated patients")
title_art_pairwise <- paste("Table i. Pairwise Non-parametric ANOVA (ART) of molecular scores in biopsies from treated vs untreated patients")

res_art_flextable <- df_univariate_02 %>%
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
    flextable::bg(i = ~ FDR < 0.05 & Term == "Group:Felz", j = 3:6, bg = "#fbff00") %>%
    flextable::colformat_double(j = 2:3, digits = 2) %>%
    flextable::colformat_double(j = 4:ncol_keys(.), digits = 3) %>%
    flextable::border_remove() %>%
    flextable::bold(part = "header") %>%
    flextable::padding(padding = 0, part = "all") %>%
    flextable::border(border = officer::fp_border(), part = "all") %>%
    flextable::autofit()

# res_art_flextable %>% print(preview = "pptx")


# FORMAT TABLES OF CONTRASTS ####
data_pairwise_formatted_delta <- df_univariate_02 %>%
    dplyr::select(annotation, score, art_con_interaction_default_tidy, art_con_cld, medians_delta) %>%
    unnest(c(art_con_interaction_default_tidy, art_con_cld), names_repair = tidyr_legacy) %>%
    unnest(medians_delta, names_repair = tidyr_legacy) %>%
    dplyr::rename(p_interaction = adj.p.value) %>%
    mutate(
        FDR_interaction = p_interaction %>% p.adjust(method = "fdr"), .by = annotation
    ) %>%
    mutate_at(
        vars(contains("p_i"), contains("FDR")),
        ~ ifelse(
            . < 0.01,
            formatC(., digits = 0, format = "e"),
            formatC(., digits = 3, format = "f")
        )
    ) %>%
    mutate(
        delta = paste(
            format(round(median_delta, 2), nsmall = 1),
            "\u00B1",
            round(IQR_delta, 2),
            # Letter,
            sep = " "
        ),
        deltadelta = paste(
            format(round(median_delta_delta, 2), nsmall = 1),
            "\u00B1",
            round(IQR_delta_delta, 2),
            # Letter,
            sep = " "
        )
    ) %>%
    dplyr::select(
        annotation, score, Group_pairwise, Group_pairwise1, FDR_interaction,
        Felz, delta, deltadelta
    ) %>%
    distinct(Group_pairwise1, score, Felz, .keep_all = TRUE) %>%
    nest(.by = c(annotation, score, Group_pairwise1)) %>%
    mutate(
        data = map(
            data,
            function(data) {
                data %>%
                    pivot_wider(names_from = Felz, values_from = delta) %>%
                    relocate(deltadelta, .after = last_col())
            }
        )
    )

data_delta_formatted <- df_univariate_02 %>%
    dplyr::select(annotation, score, art_con_interaction_default_tidy, art_con_cld, medians_delta) %>%
    unnest(medians_delta, names_repair = tidyr_legacy) %>%
    mutate_at(
        vars(contains("p_i"), contains("FDR")),
        ~ ifelse(
            . < 0.01,
            formatC(., digits = 0, format = "e"),
            formatC(., digits = 3, format = "f")
        )
    ) %>%
    mutate(
        delta = paste(
            format(round(median_delta, 2), nsmall = 1),
            "\u00B1",
            round(IQR_delta, 2),
            sep = " "
        ),
        deltadelta = paste(
            format(round(median_delta_delta, 2), nsmall = 1),
            "\u00B1",
            round(IQR_delta_delta, 2),
            sep = " "
        )
    ) %>%
    dplyr::select(annotation, score, Group_pairwise, Felz, delta, deltadelta) %>%
    distinct(Group_pairwise, score, Felz, .keep_all = TRUE) %>%
    nest(
        .by = c(annotation, score, Group_pairwise),
        medians = everything()
    ) %>%
    mutate(
        medians = map(
            medians,
            function(medians) {
                medians %>%
                    pivot_wider(names_from = Felz, values_from = delta) %>%
                    relocate(deltadelta, .after = last_col())
            }
        )
    )

data_res_art_formatted <- df_univariate_02 %>%
    dplyr::select(annotation, score, art_con_interaction_default_tidy, art_con_cld, medians_delta) %>%
    unnest(c(art_con_interaction_default_tidy, art_con_cld), names_repair = tidyr_legacy) %>%
    # mutate(FDR_interaction = adj.p.value %>% p.adjust(method = "fdr"), .by = annotation) %>%
    mutate(FDR_interaction = adj.p.value) %>%
    mutate_at(
        vars(contains("p_i"), contains("FDR")),
        ~ ifelse(
            . < 0.01,
            formatC(., digits = 0, format = "e"),
            formatC(., digits = 3, format = "f")
        )
    ) %>%
    dplyr::select(annotation, score, Group_pairwise, FDR_interaction, Letter) %>%
    nest(
        .by = c(annotation, score, Group_pairwise),
        contrasts = everything()
    )

data_pairwise_formatted <- data_delta_formatted %>%
    left_join(data_res_art_formatted, by = c("annotation", "score", "Group_pairwise")) %>%
    mutate(
        data = pmap(
            list(medians, contrasts),
            function(medians, contrasts) {
                medians %>%
                    left_join(contrasts, by = c("annotation", "score", "Group_pairwise")) %>%
                    # mutate(deltadelta = deltadelta %>% paste(Letter)) %>%
                    dplyr::select(-annotation, -score, -Group_pairwise, -Letter)
            }
        )
    ) %>%
    dplyr::select(-medians, -contrasts) %>%
    pivot_wider(names_from = Group_pairwise, values_from = data) %>%
    relocate(`FU1 - FU2`, .before = `Index - FU2`)

data_pairwise_formatted %>%
    dplyr::slice(1) %>%
    pull(4)

data_pairwise_formatted %>%
    dplyr::slice(1) %>%
    unnest(everything(), names_repair = tidyr_legacy)


# FORMAT FLEXTABLE ####
# define sample sizes
# df_n <- df_00 %>%
#     group_by(Group_Felz) %>%
#     tally()

# N <- df_n %>% pull(n)
title_art_pairwise <- paste("Table i. Median \u00B1 IQR molecular scores in biopsies from treated vs untreated patients")

footnoteText <- c(
    paste(
        "Grey shading denotes ANOVA interactive effect FDR < 0.05\n",
        "Bold denotes FDR < 0.05\n",
        "FDR correction was carried out within each annotation grouping",
        sep = ""
    ) # ,
    # "Scores are the mean fold change in expression for all Probes within the set vs the mean expression for all Probes in the NoPGD control set"
)

header1 <- c(
    "Annotation", "Score",
    rep("FU1 - Index", 4),
    rep("FU2 - FU1", 4),
    rep("FU2 - Index", 4)
)

header2 <- c(
    "Annotation", "Score",
    rep(c("\u394 NoFelz\n(N=11)", "\u394 Felz\n(N=11)", "\u394\u394", "\u394\u394 FDR"), 1),
    rep(c("\u394 NoFelz\n(N=10)", "\u394 Felz\n(N=10)", "\u394\u394", "\u394\u394 FDR"), 2)
)

cellWidths <- c(4, 11, rep(c(3, 3, 3, 2), 3))

category_vec1 <- str_remove(data_pairwise_formatted$score, "Prob\\)")
category_vec2 <- str_extract(data_pairwise_formatted$score, "Prob")
category_vec3 <- str_replace(category_vec2, "Prob", ")")


# CREATE FELXTABLE OF PAIRWISE ART MODELS ####
flextable_pairwise <- data_pairwise_formatted %>%
    unnest(everything(), names_repair = tidyr_legacy) %>%
    flextable::flextable() %>%
    flextable::compose(
        j = "score",
        value = as_paragraph(category_vec1, as_sub(category_vec2), category_vec3)
    ) %>%
    flextable::delete_part("header") %>%
    flextable::add_header_row(top = TRUE, values = header2) %>%
    flextable::add_header_row(top = TRUE, values = header1) %>%
    flextable::add_header_row(top = TRUE, values = rep(title_art_pairwise, ncol_keys(.))) %>%
    flextable::add_footer_row(values = footnoteText[[1]], colwidths = ncol_keys(.)) %>%
    flextable::merge_v(j = 1:2) %>%
    flextable::merge_v(part = "header") %>%
    flextable::merge_h(part = "header") %>%
    flextable::border_remove() %>%
    flextable::border(part = "header", border = officer::fp_border()) %>%
    flextable::border(part = "body", border = officer::fp_border()) %>%
    flextable::border(part = "footer", border.left = officer::fp_border(), border.right = officer::fp_border()) %>%
    flextable::border(i = 1, part = "footer", border.bottom = officer::fp_border()) %>%
    flextable::align(align = "center") %>%
    flextable::align(align = "center", part = "header") %>%
    # flextable::valign(i = 3, j = ncol_keys(.), valign = "bottom", part = "header") %>%
    flextable::font(fontname = "Arial", part = "all") %>%
    flextable::fontsize(size = 8, part = "all") %>%
    flextable::fontsize(size = 8, part = "footer") %>%
    flextable::fontsize(i = 1, size = 12, part = "header") %>%
    flextable::bold(part = "header") %>%
    flextable::bold(j = 1, part = "body") %>%
    flextable::bg(bg = "white", part = "all") %>%
    flextable::bg(i = ~ as.numeric(FDR_interaction) < 0.05, j = 2:ncol_keys(.), bg = "grey90", part = "body") %>%
    flextable::bold(i = ~ as.numeric(FDR_interaction) < 0.05, j = ncol_keys(.), part = "body") %>%
    flextable::padding(padding = 0, part = "all") %>%
    flextable::width(width = cellWidths, unit = "cm") %>%
    flextable::width(., width = dim(.)$widths * 30 / (flextable_dim(.)$widths), unit = "cm")


# PRINT THE FLEXTABLES ####
flextable_pairwise %>% print(preview = "pptx")


# PLOTTING GLOBALS ####
dodge <- 0.3


# MAKE BIOPSY PAIR PLOTS ####
plot_violin_pairs <- df_univariate_02 %>%
    mutate(
        plot_violin = pmap(
            list(data, variable, score, medians, medians_delta, art_con_interaction_default_tidy),
            gg_violin_interaction
        )
    )



# MAKE JOINT PATIENT PAIR PLOTS ####
plot_patient_pairs <- df_univariate_02 %>%
    mutate(
        plot_patient = pmap(
            list(data, variable, score),
            function(data, variable, score) {
                data <- data %>%
                    mutate(
                        variable = variable,
                        Felz = Felz %>% factor(labels = c("No Felzartamab", "Felzartamab"))
                    )
                data %>%
                    ggplot(
                        aes(
                            x = Group,
                            y = value,
                            col = Patient,
                            group = Patient
                        )
                    ) +
                    geom_point(size = 5, position = position_dodge(width = dodge)) +
                    geom_line(
                        linewidth = 0.75,
                        linetype = "dashed",
                        position = position_dodge(width = dodge),
                        show.legend = FALSE
                    ) +
                    geom_text(
                        aes(label = Patient),
                        size = 4,
                        col = "black",
                        show.legend = FALSE,
                        position = position_dodge(width = dodge)
                    ) +
                    labs(
                        x = NULL,
                        y = score %>% str_replace("\\(", "\n("),
                        col = "STUDY_EVALUATION_ID     ",
                        parse = TRUE
                    ) +
                    coord_cartesian(xlim = c(1.3, 2.7)) +
                    theme_bw() +
                    theme(
                        axis.title = element_text(size = 15),
                        axis.text = element_text(size = 12, colour = "black"),
                        panel.grid = element_blank(),
                        strip.text = element_text(size = 12, colour = "black"),
                        legend.position = "top",
                        legend.title = element_text(size = 15, hjust = 0.66, vjust = 1, face = "bold"),
                        legend.text = element_text(size = 15),
                    ) +
                    facet_wrap(data$Felz) +
                    guides(col = guide_legend(nrow = 1))
            }
        )
    )



# MAKE JOINT BIOPSY PAIR PLOTS ####
panel_pairs <- plot_violin_pairs %>%
    dplyr::filter(category %in% c("ABMR", "TCMR")) %>%
    pull(plot_violin) %>%
    ggarrange(
        plotlist = .,
        common.legend = TRUE,
        ncol = 5,
        nrow = 2,
        align = "hv",
        labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J"),
        font.label = list(size = 20, color = "black", face = "bold")
    )


panel_patient_pairs <- plot_patient_pairs %>%
    dplyr::filter(category %in% c("ABMR", "TCMR")) %>%
    pull(plot_patient) %>%
    ggarrange(
        plotlist = .,
        common.legend = TRUE,
        ncol = 5,

        nrow = 2,
        align = "hv",
        labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J"),
        font.label = list(size = 20, color = "black", face = "bold")
    )



# SAVE THE PLOTS ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0Georg Felz CD38 Vienna/G_Rstuff/output/"
# ggsave(
#     filename = paste(saveDir, "Felzartamab scores 5086.png"),
#     plot = panel_pairs,
#     dpi = 600,
#     width = 60,
#     height = 22,
#     units = "cm",
#     bg = "white"
# )
# ggsave(
#     filename = paste(saveDir, "Felzartamab patient pairs 5086.png"),
#     plot = panel_patient_pairs,
#     dpi = 600,
#     width = 60,
#     height = 22,
#     units = "cm",
#     bg = "white"
# )



# END ####
