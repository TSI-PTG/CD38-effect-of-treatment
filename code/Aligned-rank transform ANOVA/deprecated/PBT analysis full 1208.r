# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(ggbeeswarm) # install.packages("ggbeeswarm")
library(ggpubr) # install.packages("ggpubr")
library(ggrepel) # install.packages("ggrepel")
library(gghalves) # pak::pak("erocoar/gghalves")
library(ggh4x) # install.packages("ggh4x")
library(broom) # install.packages("broom") #for tabular model object transformations
library(car) # install.packages("car") #for type-II manovas
library(rstatix) # install.packages("rstatix") #for testing manova assumptions
library(vegan) # install.packages("vegan") #for permanova
library(pairwiseAdonis) # library(devtools) #install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
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
# laod reference set
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0Georg Felz CD38 Vienna/G_Rstuff/data/vienna_1208_16Feb24.RData")


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
features <- c(Rejectionrelated, ABMRrelated, Macrophage)


# DEFINE THE SET ####
set <- vienna_1208


# WRANGLE THE PHENOTYPE DATA ####
df00 <- set %>%
    pData() %>%
    dplyr::rename(
        Patient = STUDY_EVALUATION_ID,
        Felz = Felzartamab_presumed
    ) %>%
    mutate(
        Patient = Patient %>% factor(),
        Group = Group %>% factor(levels = c("Index", "FU1")),
        Felz = Felz %>% factor(labels = c("NoFelz", "Felz")),
        TxBx = TxBx %>% as.numeric(),
        Group_Felz = paste(Group, Felz, sep = ":") %>%
            factor(levels = c("Index:NoFelz", "FU1:NoFelz", "Index:Felz", "FU1:Felz"))
    ) %>%
    expand_grid(category = c("ABMR", "TCMR", "Macrophage", "Injurylate")) %>%
    group_by(category) %>%
    nest() %>%
    tibble() %>%
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
                } else if (category == "Injurylate") Injurylate
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
df01 <- df00 %>%
    mutate(
        permanova = pmap(
            list(data, features),
            function(data, features) {
                set.seed(seed)
                features <- data %>% dplyr::select(all_of(features))
                Group <- data$Group
                Felz <- data$Felz
                adonis2(
                    features ~ Group * Felz,
                    data = data,
                    method = "euclidean",
                    # by = "margin", only specific margin if the sample sizes are unequal
                    permutations = 100000
                )
            }
        )
    )
df01$permanova


# PAIRWISE PERMANOVA ####
df02 <- df01 %>%
    mutate(
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
df02$permanova_pairwise


# PRODUCE TABLE OF PERMANOVA RESULTS ####
df03 <- df02 %>%
    mutate(
        permanova_flextable = pmap(
            list(permanova),
            function(permanova) {
                title <- paste("Table i. PERMANOVA of molecular scores in biopsies from treated vs untreated patients")
                permanova %>%
                    tidy() %>%
                    suppressWarnings() %>%
                    dplyr::filter(term %nin% c("Residual", "Total")) %>%
                    dplyr::select(-df, -SumOfSqs) %>%
                    mutate(p.value = p.value %>% formatC(digits = 3, format = "fg")) %>%
                    flextable::flextable() %>%
                    flextable::add_header_row(values = rep(title, ncol_keys(.))) %>%
                    flextable::merge_h(part = "header") %>%
                    flextable::fontsize(size = 8, part = "all") %>%
                    flextable::align(align = "center", part = "all") %>%
                    flextable::bg(bg = "white", part = "all") %>%
                    flextable::colformat_double(j = 2:3, digits = 2) %>%
                    flextable::bg(i = ~ p.value %>% as.numeric() < 0.05, bg = "#fbff00", part = "body") %>%
                    flextable::border_remove() %>%
                    flextable::bold(part = "header") %>%
                    flextable::padding(padding = 0, part = "all") %>%
                    flextable::border(border = fp_border(), part = "all") %>%
                    flextable::autofit()
            }
        )
    )
df03$permanova_flextable %>% print(preview = "pptx")


# PRODUCE TABLE OF PAIRWISE ADONIS RESULTS ####
df04 <- df03 %>%
    mutate(
        permanova_pairwise_flextable = pmap(
            list(permanova_pairwise),
            function(permanova_pairwise) {
                title <- paste("Table i. Pairwise PERMANOVA of molecular scores in biopsies from treated vs untreated patients")
                permanova_pairwise %>%
                    dplyr::select(-Df, -sig, -SumsOfSqs) %>%
                    dplyr::rename("F-value" = F.Model, FDR = p.adjusted) %>%
                    mutate(
                        p.value = p.value %>% formatC(digits = 0, format = "e"),
                        FDR = FDR %>% formatC(digits = 0, format = "e"),
                    ) %>%
                    arrange(p.value) %>%
                    flextable::flextable() %>%
                    flextable::add_header_row(values = rep(title, ncol_keys(.))) %>%
                    flextable::merge_h(part = "header") %>%
                    flextable::fontsize(size = 8, part = "all") %>%
                    flextable::align(align = "center", part = "all") %>%
                    flextable::bg(bg = "white", part = "all") %>%
                    flextable::bg(i = ~ FDR %>% as.numeric() < 0.05, bg = "#fbff00") %>%
                    flextable::colformat_double(j = 2:3, digits = 2) %>%
                    flextable::border_remove() %>%
                    flextable::bold(part = "header") %>%
                    flextable::padding(padding = 0, part = "all") %>%
                    flextable::border(border = fp_border(), part = "all") %>%
                    flextable::autofit()
            }
        )
    )
df04$permanova_pairwise_flextable
# df04$permanova_pairwise_flextable  %>% print(preview = "docx")


# WRANGLE THE DATA FOR UNIVARIATE TESTS ####
df_univariate_00 <- df04 %>%
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
                            variable == "QCAT" ~ "Cytotoxic T cell-associated transcripts (QCAT)",
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
    unnest(data_univariate)


# UNIVARIATE MEANS AND MEDIANS ####
df_univariate_01 <- df_univariate_00 %>%
    mutate(
        n_cat = score %>% unique() %>% length(), .by = category, .after = variable
    ) %>%
    mutate(
        means = map(
            data,
            function(data) {
                data %>%
                    reframe(mean = mean(value), se = se(value), .by = Group_Felz) %>%
                    arrange(
                        Group_Felz %>%
                            factor(levels = c(
                                "Index:NoFelz", "FU1:NoFelz",
                                "Index:Felz", "FU1:Felz"
                            ))
                    )
            }
        ),
        means_delta = map(
            data,
            function(data) {
                data %>%
                    reframe(d = diff(value), .by = c(Patient, Felz)) %>%
                    reframe(mean_delta = mean(d), se_delta = se(d), .by = c(Felz)) %>%
                    arrange(Felz %>% factor(levels = c("NoFelz", "Felz")))
            }
        ),
        means_delta_delta = map(
            means_delta, reframe,
            mean_delta_delta = diff(mean_delta), se_delta_delta = mean(se_delta)
        ),
        medians = map(
            data,
            function(data) {
                data %>%
                    reframe(median = median(value), IQR = IQR(value), .by = c(Group, Felz)) %>%
                    mutate(
                        Group_Felz = paste(Group, Felz, sep = ":") %>%
                            factor(levels = c(
                                "Index:NoFelz", "FU1:NoFelz",
                                "Index:Felz", "FU1:Felz"
                            ))
                    ) %>%
                    arrange(Group_Felz) %>%
                    mutate(
                        median_delta = median %>% diff(),
                        IQR_delta = IQR %>% mean(),
                        .by = Felz
                    ) %>%
                    mutate(
                        median_delta_delta = median_delta %>% diff(),
                        IQR_delta_delta = IQR_delta %>% mean(),
                        .by = Group
                    )
            }
        ),
        medians_delta_inner = map(
            data,
            function(data) {
                data %>%
                    reframe(d = diff(value), .by = c(Patient, Felz)) %>%
                    reframe(median_delta = median(d), IQR_delta = IQR(d), .by = c(Felz)) %>%
                    arrange(Felz %>% factor(levels = c("NoFelz", "Felz")))
            }
        ),
        medians_delta_delta_inner = map(
            medians_delta_inner, reframe,
            median_delta_delta = diff(median_delta), IQR_delta_delta = mean(IQR_delta)
        )
    )

df_univariate_01$medians[[17]]
df_univariate_01$medians_delta_inner[[17]]
df_univariate_01$medians_delta_delta_inner[[17]]



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
            formula = "Group:Felz", adjust = "fdr", method = "pairwise", interaction = TRUE, response = "aligned"
        ),
        art_con_interaction = map(
            art,
            art.con,
            formula = "Group:Felz", adjust = "fdr", method = list("interaction" = c(1, -1, -1, 1)), response = "aligned"
        ),
        # art_con_groupwise = map(
        #     art,
        #     art.con,
        #     formula = "Group:Felz", adjust = "fdr", method = list("abs" = c(3, -1, -1, -1)), response = "aligned"
        # ),
        art_con_interaction_default_tidy = map(art_con_interaction_default, tidy),
        art_con_interaction_tidy = map(art_con_interaction, tidy),
        # art_con_groupwise_tidy = map(art_con_groupwise, tidy),
        # art_con_cld = map(
        #     art_con,
        #     function(art_con) {
        #         art_con %>%
        #             as.data.frame() %>%
        #             cldList(p.value ~ contrast, data = .) %>%
        #             arrange(Group %>%
        #                 factor(
        #                     levels = c(
        #                         "Index,NoFelz", "FU1,NoFelz",
        #                         "Index,Felz", "FU1,Felz"
        #                     )
        #                 ))
        #     }
        # )
    )

names(df_univariate_02$art_aov_tidy) <- df_univariate_02$variable %>% as.character()

names(df_univariate_02$art_con_interaction_default_tidy) <- df_univariate_02$variable %>% as.character()
names(df_univariate_02$art_con_interaction_tidy) <- df_univariate_02$variable %>% as.character()
names(df_univariate_02$art_con_groupwise_tidy) <- df_univariate_02$variable %>% as.character()
# names(df_univariate_02$art_con_cld) <- df_univariate_02$variable %>% as.character()

df_univariate_02$art_con_interaction_tidy[[1]]
df_univariate_02$art_con_interaction_default_tidy[[1]]
df_univariate_02$art_con_groupwise_tidy[[1]]
df_univariate_02$art_con_interaction[[1]]@misc


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
    flextable::border(border = fp_border(), part = "all") %>%
    flextable::autofit()

# res_art_flextable %>% print(preview = "pptx")


# FORMAT TABLES OF CONTRASTS ####
res_art_pairwise_formatted_delta <- df_univariate_02 %>%
    dplyr::select(annotation, score, art_aov_tidy, art_lm, medians) %>%
    unnest(c(art_aov_tidy, art_lm), names_repair = tidyr_legacy) %>%
    unnest(medians) %>%
    dplyr::filter(
        Term == "Group:Felz",
        Term_lm == "Group1:Felz1"
    ) %>%
    dplyr::rename(p_interaction = p.value) %>%
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
        Estimate = Estimate %>% round(3),
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
    distinct(Felz, .by = score, .keep_all = TRUE) %>%
    dplyr::select(annotation, score, Estimate, FDR_interaction, Felz, delta, deltadelta) %>%
    pivot_wider(names_from = c(Felz), values_from = c(delta)) %>%
    as.data.frame() %>%
    relocate(deltadelta, .after = last_col())

res_art_pairwise_formatted_medians <- df_univariate_02 %>%
    dplyr::select(annotation, score, art_con_groupwise_tidy, medians) %>%
    unnest(c(art_con_groupwise_tidy, medians), names_repair = tidyr_legacy) %>%
    dplyr::rename(p_groupwise = p.value) %>%
    mutate(
        FDR_groupwise = p_groupwise %>% p.adjust(method = "fdr"), .by = annotation
    ) %>%
    mutate_at(
        vars(contains("p_g"), contains("FDR")),
        ~ ifelse(
            . < 0.01,
            formatC(., digits = 0, format = "e"),
            formatC(., digits = 3, format = "f")
        )
    ) %>%
    mutate(
        medians =
            paste(
                format(round(median, 2), nsmall = 1),
                "\u00B1",
                round(IQR, 2),
                # Letter,
                sep = " "
            ),
    ) %>%
    dplyr::select(annotation, score, Group_Felz, medians, p_groupwise, FDR_groupwise) %>%
    pivot_wider(names_from = c(Group_Felz), values_from = c(medians)) %>%
    relocate(contains("groupwise"), .after = last_col()) %>%
    as.data.frame()

res_art_pairwise_formatted <- left_join(
    res_art_pairwise_formatted_delta,
    res_art_pairwise_formatted_medians,
    by = c("score", "annotation")
)


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
    rep("Aligned Rank-Transform ANOVA Interaction", 5)
)
header2 <- c(
    "Annotation", "Score",
    rep("Aligned Rank-Transform ANOVA Interaction", 5)
)
header3 <- c(
    "Annotation", "Score",
    "effect\n(aligned)", "FDR",
    "\u394\n(FU1 - Index)", "\u394\n(FU1 - Index)", "\u394\u394\n(N = 44)"
)
header4 <- c(
    "Annotation", "Score",
    "effect\n(aligned)", "FDR",
    "\u394\n(FU1 - Index)", "\u394\n(FU1 - Index)", "\u394\u394\n(N = 44)"
)

header5 <- c(
    "Annotation", "Score",
    "effect\n(aligned)", "FDR",
    "NoFelz\n(N = 22)", "Felz\n(N = 22)", "\u394\u394\n(N = 44)"
)

cellWidths <- c(4, 11, 2.2, rep(2, 1), rep(4, 3))

category_vec1 <- str_remove(res_art_pairwise_formatted$score, "Prob\\)")
category_vec2 <- str_extract(res_art_pairwise_formatted$score, "Prob")
category_vec3 <- str_replace(category_vec2, "Prob", ")")


# CREATE FELXTABLE OF PAIRWISE ART MODELS ####
res_art_pairwise_flextable <- res_art_pairwise_formatted_delta %>%
    flextable::flextable() %>%
    flextable::compose(
        j = "score",
        value = as_paragraph(category_vec1, as_sub(category_vec2), category_vec3)
    ) %>%
    flextable::delete_part("header") %>%
    flextable::add_header_row(top = TRUE, values = header5) %>%
    flextable::add_header_row(top = TRUE, values = header4) %>%
    flextable::add_header_row(top = TRUE, values = header3) %>%
    flextable::add_header_row(top = TRUE, values = header2) %>%
    flextable::add_header_row(top = TRUE, values = header1) %>%
    flextable::add_header_row(top = TRUE, values = rep(title_art_pairwise, ncol_keys(.))) %>%
    flextable::add_footer_row(values = footnoteText[[1]], colwidths = ncol_keys(.)) %>%
    flextable::merge_v(j = 1:2) %>%
    flextable::merge_v(part = "header") %>%
    flextable::merge_h(part = "header") %>%
    flextable::border_remove() %>%
    flextable::border(part = "header", border = fp_border()) %>%
    flextable::border(part = "body", border = fp_border()) %>%
    flextable::border(part = "footer", border.left = fp_border(), border.right = fp_border()) %>%
    flextable::border(i = 1, part = "footer", border.bottom = fp_border()) %>%
    flextable::align(align = "center") %>%
    flextable::align(align = "center", part = "header") %>%
    flextable::valign(i = 4, j = ncol_keys(.), valign = "bottom", part = "header") %>%
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
res_art_pairwise_flextable %>% print(preview = "pptx")



# PLOTTING GLOBALS ####
dodge <- 0.3



df_univariate_02$art_aov_tidy[[1]]

df_univariate_02$data[[1]] %>%
    left_join(
        df_univariate_02$data[[1]] %>%
            reframe(median = median(value), .by = c(Group, Felz, Patient)) %>%
            spread(Group, median) %>%
            mutate(delta = FU1 - Index) %>%
            pivot_longer(cols = c(Index, FU1), names_to = "Group", values_to = "value") %>%
            # mutate(delta = ifelse(delta < 0, "improved", "worsened")) %>%
            dplyr::select(Felz, Patient, Group, delta),
        by = c("Felz", "Patient", "Group")
    ) %>%
    left_join(df_univariate_02$medians[[1]] %>% dplyr::select(Group, Felz, median_delta_delta), by = c("Felz", "Group")) %>%
    mutate(
        Group = Group %>% factor(levels = c("Index", "FU1")),
        delta = delta %>% scale() %>% as.vector()
    )






# MAKE BIOPSY PAIR PLOTS ####
plot_violin_pairs <- df_univariate_02 %>%
    mutate(
        gg_line = pmap(
            list(data, variable, score, medians, art_aov_tidy),
            function(data, variable, score, medians, art_aov_tidy) {
                delta_delta <- medians %>%
                    mutate(
                        Felz = Felz %>% factor(labels = c("No Felzartamab", "Felzartamab"))
                    )
                delta_delta_fdr <- art_aov_tidy %>%
                    dplyr::filter(term == "Group:Felz") %>%
                    mutate(FDR = p.value %>% p.adjust(method = "fdr")) %>%
                    mutate_at(
                        vars(contains("p_i"), contains("FDR")),
                        ~ ifelse(
                            . < 0.01,
                            formatC(., digits = 0, format = "e"),
                            formatC(., digits = 3, format = "f")
                        )
                    ) %>%
                    pull(FDR)
                data <- data %>%
                    mutate(
                        variable = variable,
                        Felz = Felz %>% factor(labels = c("No Felzartamab", "Felzartamab"))
                    )
                medians <- data %>%
                    reframe(median = median(value), .by = c(Group, Felz)) %>%
                    spread(Group, median) %>%
                    mutate(delta = FU1 - Index)
                means <- data %>%
                    reframe(mean = median(value), .by = c(Group, Felz)) %>%
                    spread(Group, mean) %>%
                    mutate(delta = FU1 - Index)
                data <- data %>%
                    left_join(
                        data %>%
                            reframe(median = median(value), .by = c(Group, Felz, Patient)) %>%
                            spread(Group, median) %>%
                            mutate(delta = FU1 - Index) %>%
                            pivot_longer(cols = c(Index, FU1), names_to = "Group", values_to = "value") %>%
                            # mutate(delta = ifelse(delta < 0, "improved", "worsened")) %>%
                            dplyr::select(Felz, Patient, Group, delta),
                        by = c("Felz", "Patient", "Group")
                    ) %>%
                    left_join(delta_delta %>% dplyr::select(Group, Felz, median_delta_delta), by = c("Felz", "Group")) %>%
                    mutate(
                        Group = Group %>% factor(levels = c("Index", "FU1")),
                        # delta = delta %>% scale() %>% as.vector()
                    )
                data %>%
                    ggplot(
                        aes(
                            x = Group,
                            y = value,
                            # group = Patient,
                            # linetype = Felz
                        )
                    ) +
                    geom_half_violin(
                        inherit.aes = FALSE,
                        aes(
                            x = Group,
                            y = value
                        ),
                        # draw_quantiles = c(0.25, 0.75),
                        side = c("l", "r"),
                        fill = "grey95",
                        trim = FALSE,
                        scale = "width"
                    ) +
                    geom_point(
                        aes(col = delta),
                        position = position_nudge(x = ifelse(data$Group == "Index", 0.1, -0.1)),
                        size = 2, alpha = 0.75
                    ) +
                    geom_line(
                        aes(
                            # x = Group,
                            col = delta,
                            group = Patient
                        ),
                        position = position_nudge(x = ifelse(data$Group == "Index", 0.1, -0.1)),
                        linewidth = 0.5, alpha = 0.25,
                        show.legend = FALSE
                    ) +
                    # geom_text(
                    #     aes(label = Patient),
                    #     size = 4,
                    #     col = "#00000076",
                    #     show.legend = FALSE
                    # ) +
                    stat_summary(fun = median, geom = "point", size = 4) +
                    stat_summary(
                        aes(group = Felz),
                        fun = median,
                        geom = "line",
                        linewidth = 1,
                        linetype = "dashed"
                    ) +
                    geom_errorbar(
                        aes(x = 2.55, ymin = Index, ymax = FU1, y = NULL),
                        data = medians, width = 0.05, linewidth = 0.5, color = "black"
                    ) +
                    geom_text(
                        aes(
                            x = 2.7, y = (Index + FU1) / 2,
                            label = paste("\u394 =", round(delta, 2))
                        ),
                        data = medians, hjust = 0.5, vjust = 1, size = 5, color = "black", angle = 90
                    ) +
                    labs(
                        title = paste(
                            "\u394\u394 =",
                            data$median_delta_delta %>% round(2),
                            "| FDR =",
                            delta_delta_fdr
                        ),
                        x = NULL,
                        y = score %>% str_replace("\\(", "\n("),
                        col = "Individual patient response     ",
                        parse = TRUE
                    ) +
                    scale_color_gradient2(
                        low = "#00ff00bc",
                        mid = "grey60",
                        high = "red",
                        midpoint = 0,
                        breaks = c(min(data$delta), max(data$delta)),
                        labels = c("improved", "worsened"),
                        guide = guide_colorbar(
                            title.position = "top",
                            barwidth = 20,
                            ticks = FALSE,
                            label.hjust = c(1.1, -0.1),
                            label.vjust = 8,
                            reverse = TRUE
                        ),
                    ) +
                    scale_x_discrete(labels = c("Index\n(n = 11)", "Follow-up\n(n = 11)")) +
                    coord_cartesian(xlim = c(NA, 2.6)) +
                    theme_bw() +
                    theme(
                        axis.title = element_text(size = 15),
                        axis.text = element_text(size = 12, colour = "black"),
                        panel.grid = element_blank(),
                        strip.text = element_text(size = 12, colour = "black"),
                        legend.position = "top",
                        legend.title = element_text(size = 15, hjust = 0.5, vjust = 1, face = "bold"),
                        legend.text = element_text(size = 15),
                        plot.title = element_text(
                            colour = "black",
                            hjust = 0.5,
                            size = 12,
                            face = ifelse(delta_delta_fdr %>% as.numeric() < 0.05, "bold.italic", "italic")
                        )
                    ) +
                    facet_wrap(~Felz)
                # +
                # guides(col = guide_legend(nrow = 1, override.aes = list(size = 6)))
            }
        )
    )

# plot_violin_pairs$gg_line[[1]]


# MAKE JOINT BIOPSY PAIR PLOTS ####
panel_pairs <- plot_violin_pairs %>%
    dplyr::filter(category %in% c("ABMR", "TCMR")) %>%
    pull(gg_line) %>%
    ggarrange(
        plotlist = .,
        common.legend = TRUE, ncol = 5, nrow = 2,
        labels = c("A", "", "", "", "", "B"), font.label = list(size = 20, color = "black", face = "bold")
    )



# MAKE VIOLIN PLOTS ####
plot_violin <- df_univariate_02 %>%
    dplyr::select(category, annotation, score, variable, data, medians, means) %>%
    # dplyr::filter(category == "ABMR") %>%
    mutate(
        gg_violin = pmap(
            list(data, variable, medians, score),
            function(data, variable, medians, score) {
                data <- data %>%
                    mutate(
                        variable = variable,
                        Felz = Felz %>% factor(labels = c("No Felzartamab", "Felzartamab"))
                    )
                medians <- data %>%
                    reframe(median = median(value), .by = c(Group, Felz)) %>%
                    spread(Group, median) %>%
                    mutate(delta = FU1 - Index)
                means <- data %>%
                    reframe(mean = median(value), .by = c(Group, Felz)) %>%
                    spread(Group, mean) %>%
                    mutate(delta = FU1 - Index)
                data %>%
                    ggplot(aes(x = Group, y = value)) +
                    geom_violin(
                        scale = "width",
                        draw_quantiles = c(0.25, 0.75),
                        trim = FALSE, adjust = 1,
                        fill = "grey95"
                    ) +
                    geom_point(pch = "-", size = 4) +
                    stat_summary(fun = median, geom = "point", size = 4) +
                    stat_summary(
                        aes(group = Felz),
                        fun = median,
                        geom = "line",
                        linewidth = 1,
                        linetype = "dashed"
                    ) +
                    geom_errorbar(
                        aes(x = 2.55, ymin = Index, ymax = FU1, y = NULL),
                        data = medians, width = 0.05, linewidth = 0.5, color = "black"
                    ) +
                    geom_text(
                        aes(
                            x = 2.7, y = (Index + FU1) / 2,
                            label = paste("\u394 =", round(delta, 2))
                        ),
                        data = medians, hjust = 0.5, vjust = 1, size = 5, color = "black", angle = 90
                    ) +
                    scale_x_discrete(labels = c("Index\n(N = 11)", "FU1\n(N = 11)")) +
                    labs(
                        # title = annotation,
                        x = NULL,
                        y = score %>% str_replace("\\(", "\n("),
                        parse = TRUE
                    ) +
                    coord_cartesian(xlim = c(NA, 2.6)) +
                    theme_bw() +
                    theme(
                        axis.title = element_text(size = 15),
                        axis.text = element_text(size = 12, colour = "black"),
                        panel.grid = element_blank(),
                        legend.position = "top",
                        strip.text = element_text(size = 15, colour = "black"),
                        plot.title = element_text(colour = "red", size = 20, face = "bold")
                    ) +
                    facet_wrap(~Felz)
            }
        )
    )


# MAKE JOINT VIOLIN PLOTS ####
panel_ABMR <- plot_violin %>%
    dplyr::filter(category == "ABMR") %>%
    pull(gg_violin) %>%
    ggarrange(
        plotlist = ., common.legend = TRUE, nrow = 1
        # labels = "A", font.label = list(size = 20, color = "black", face = "bold")
    )

panel_TCMR <- plot_violin %>%
    dplyr::filter(category == "TCMR") %>%
    pull(gg_violin) %>%
    ggarrange(
        plotlist = ., common.legend = TRUE, nrow = 1
        # labels = "B", font.label = list(size = 20, color = "black", face = "bold")
    )

panel_violin <- ggarrange(
    panel_ABMR, panel_TCMR,
    common.legend = TRUE, nrow = 2,
    labels = c("A", "B"), font.label = list(size = 20, color = "black", face = "bold")
)




# SAVE THE PLOTS ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0Georg Felz CD38 Vienna/G_Rstuff/output/"
# ggsave(
#     filename = paste(saveDir, "Felzartamab effect of treatment patient data.png"),
#     plot = panel_pairs,
#     dpi = 300,
#     width = 60,
#     height = 22,
#     units = "cm",
#     bg = "white"
# )

# ggsave(
#     filename = paste(saveDir, "Felzartamab effect of treatment molecular scores.png"),
#     plot = panel_violin,
#     dpi = 300,
#     width = 60,
#     height = 20,
#     units = "cm",
#     bg = "white"
# )



# ggsave(
#     filename = paste(saveDir, "Felzartamab effect of treatment ABMR scores.png"),
#     plot = panel_ABMR,
#     dpi = 300,
#     width = 60,
#     height = 10,
#     units = "cm",
#     bg = "white"
# )
# ggsave(
#     filename = paste(saveDir, "Felzartamab effect of treatment TCMR scores.png"),
#     plot = panel_TCMR,
#     dpi = 300,
#     width = 60,
#     height = 10,
#     units = "cm",
#     bg = "white"
# )

# # END ####
