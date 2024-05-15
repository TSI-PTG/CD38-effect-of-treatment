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
library(flextable) # install.packages("flextable") #for table outputs
library(officer) # install.packages("officer")
library(emmeans) # install.packages("emmeans") #for post-hoc testing and CLD
library(multcomp) # install.packages("multcomp") #for for CLD
library(ARTool) # install.packages("ARTool") #for non-parametric anova (aligned-rank test) and post-hoc
library(rcompanion) # install.packages("rcompanion") #for non-parametric anova (aligned-rank test) and post-hoc
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
# load reference set
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0Georg Felz CD38 Vienna/G_Rstuff/data/data_K1208.RData")

# DEFINE SEED ####
seed <- 42


# DEFINE CATEGORIES FOR FEATURES ####
# Rejectionrelated <- c("GRIT3", "Rej-RAT", "RejAA_NR")
ABMRrelated <- c("DSAST", "NKB", "ABMRpm", "ggt0", "ptcgt0") # "cggt0", "RejAA_EABMR", "RejAA_FABMR", "RejAA_LABMR")
TCMRrelated <- c("QCAT", "TCB", "TCMRt", "tgt1", "igt1") # , "TCMR-RAT", )
# Endothelium <- c("ENDAT")
# Parenchyma <- c("KT1", "KT2")
# Macrophage <- c("AMAT1", "QCMAT")
# Injuryrecent <- c("FICOL", "IRRAT30", "IRITD3", "IRITD5")
# Injurylate <- c("IGT", "MCAT", "BAT", "cigt1", "ctgt1")


# DEFINE FEATURES ####
features <- c(ABMRrelated, TCMRrelated)


# DEFINE THE SET ####
set <- data_K1208 %>%
    dplyr::filter(STUDY_EVALUATION_ID %nin% c(15, 18))


data_K1208 %>% colnames()

set %>%
    ggplot(aes(x = Donor_Age, y = ABMRpm)) +
    geom_point() +
    geom_smooth(method = "gam") +
    facet_wrap(~Group)


# WRANGLE THE PHENOTYPE DATA ####
data_01 <- set %>%
    tibble() %>%
    relocate(Group, .after = Felzartamab) %>%
    dplyr::rename(
        Patient = STUDY_EVALUATION_ID,
        Felz = Felzartamab
    ) %>%
    mutate(
        Patient = Patient %>% factor(),
        Group = Group %>% factor(levels = c("Index", "FU1", "FU2")),
        Felz = Felz %>% factor(labels = c("NoFelz", "Felz")),
        # TxBx = TxBx %>% as.numeric(),
        Group_Felz = paste(Group, Felz, sep = ":") %>%
            factor(levels = c("Index:NoFelz", "FU1:NoFelz", "FU2:NoFelz", "Index:Felz", "FU1:Felz", "FU2:Felz")),
        .after = "Group"
    ) %>%
    arrange(Patient, Group) %>%
    expand_grid(category = c("ABMR", "TCMR")) %>%
    nest(.by = category) %>%
    mutate(
        features = map(
            category,
            function(category) {
                if (category == "ABMR") {
                    ABMRrelated
                } else if (category == "TCMR") {
                    TCMRrelated
                }
            }
        ),
        data_molecular = pmap(
            list(features, data),
            function(features, data) {
                data %>%
                    dplyr::select(
                        CEL, Patient, Trial_Center, Group, Felz, Group_Felz,
                        all_of(features)
                    )
            }
        ),
        data_clinical = pmap(
            list(features, data),
            function(features, data) {
                data %>%
                    dplyr::select(
                        -all_of(features)
                    )
            }
        )
    )


# WRANGLE THE DATA FOR UNIVARIATE TESTS ####
data_02 <- data_01 %>%
    mutate(
        data_univariate = pmap(
            list(data_molecular, features),
            function(data_molecular, features) {
                data_molecular %>%
                    pivot_longer(
                        cols = all_of(features),
                        names_to = "variable",
                        values_to = "value"
                    ) %>%
                    nest(.by = variable) %>%
                    mutate(
                        variable = variable %>%
                            factor(
                                levels = c(
                                    "ABMRpm", "ggt0", "cggt0", "ptcgt0", "NKB", "DSAST",
                                    "TCMRt", "tgt1", "igt1", "TCB", "TCMR-RAT", "QCAT"
                                )
                            ),
                        annotation = case_when(
                            variable %in% TCMRrelated ~ "TCMR-related",
                            variable %in% ABMRrelated ~ "ABMR-related",
                            TRUE ~ " "
                        ) %>%
                            factor(
                                levels = c(
                                    "ABMR-related",
                                    "TCMR-related"
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
    dplyr::select(category, data_univariate, data_clinical) %>%
    unnest(data_univariate)

data_02$data_clinical[[1]]


# CALCULATE CHANGE IN SCORES AMONG FOLLOW-UP ####
data_K1208_delta <- data_02 %>%
    mutate(
        delta = map(
            data,
            function(data) {
                data %>%
                    mutate(value_norm = value + min(value) + 0.01) %>%
                    mutate(
                        delta = value - lag(value),
                        delta = ifelse(delta %>% is.na(), last(value) - first(value), delta),
                        prop = value_norm / lag(value_norm),
                        prop = ifelse(prop %>% is.na(), last(value_norm) / first(value_norm), prop),
                        percent = ifelse(prop > 1, prop * 100, -(1 - prop) * 100),
                        contrast = case_when(
                            Group == "FU1" ~ "FU1 - Index",
                            Group == "FU2" ~ "FU2 - FU1",
                            TRUE ~ "FU2 - Index"
                        ) %>% factor(levels = c("FU1 - Index", "FU2 - FU1", "FU2 - Index")),
                        .by = Patient
                    )
            }
        )
    ) %>%
    mutate(data_delta = map2(data_clinical, delta, left_join)) %>%
    dplyr::select(-data, -data_clinical, -delta) %>%
    unnest(data_delta) %>%
    nest(.by = c("category", "annotation", "score", "variable", "contrast")) %>%
    arrange(variable, contrast)



# SAVE PLOTS ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0Georg Felz CD38 Vienna/G_Rstuff/data/"
save(data_K1208_delta, file = paste(saveDir, "data_K1208_patient_deltas.RData"))


# DEFINE CLINICAL FEATURES OF INTEREST ####
data_delta$data[[1]] %>% colnames()
vars_clinical_continuous <- c("Age_at_Tx", "Age_at_Screening", "Years_to_ScreeningVisit", "Screening_GFR", "Screening_Prot_Krea")
vars_clinical_discrete <- c("")

data_delta$data[[1]] %>%
    dplyr::select(Felz, DSA_HLA_class_I_Only_Screening) %>%
    arrange(Felz)


# MAKE SOME GENERAL PLOTS ####
data_plots <- data_delta %>%
    expand_grid(vars_clinical_continuous) %>%
    mutate(
        plot_scatter = pmap(
            list(variable, vars_clinical_continuous, contrast, data),
            function(variable, vars_clinical_continuous, contrast, data) {
                sym <- vars_clinical_continuous %>% sym()
                data %>%
                    mutate(
                        contrast = contrast,
                        variable = variable
                    ) %>%
                    ggplot(aes_string(x = sym, y = "percent", col = "Felz")) +
                    geom_point(
                        aes(fill = ABMRActivity_Banff19 %>% factor()),
                        pch = 21,
                        data = data %>% dplyr::filter(ABMRActivity_Banff19 == 1),
                        col = "#fa0d0d",
                        size = 2.5
                    ) +
                    geom_point() +
                    geom_smooth(
                        method = "gam",
                        formula = y ~ s(x, bs = "cs", k = 3),
                        na.rm = TRUE,
                        se = FALSE,
                        show.legend = FALSE
                    ) +
                    scale_colour_manual(values = c("black", "#005eff")) +
                    scale_fill_manual(values = c("#00000000")) +
                    labs(
                        fill = "non-responder",
                        col = NULL,
                        y = paste("%\u394", variable)
                    ) +
                    # facet_wrap(~contrast) +
                    theme_bw()
            }
        )
    )
data_plots$plot_scatter[[4]]












# MAKE JOINT PLOTS ####
plot_scatter_joint <- data_plots %>%
    dplyr::filter(category == "ABMR", contrast == "FU1 - Index") %>%
    pull(plot_scatter) %>%
    wrap_plots(.) +
    plot_layout(
        ncol = 5,
        guides = "collect",
        axes = "collect"
    ) &
    theme(
        legend.position = "top",
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15)
        ) &
        guides(col = guide_legend(
            override.aes = list(size = 8)
        ))




# SAVE PLOTS ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0Georg Felz CD38 Vienna/G_Rstuff/output/"
ggsave(
    filename = paste(saveDir, "Felzartamab scores vs time 1208.png"),
    plot = plot_scatter_joint,
    dpi = 600,
    width = 22,
    height = 25,
    units = "cm",
    bg = "white"
)
