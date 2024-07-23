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
source("C:/R/CD38-effect-of-treatment/code/functions/plot.gg_patient_pairs_interaction.r")
# load reference set
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/data_scores_k1208.RData")



# DEFINE SEED ####
seed <- 42



# DEFINE CATEGORIES FOR FEATURES ####
# Rejectionrelated <- c("GRIT3", "Rej-RAT", "RejAA_NR")
vars_cfDNA <- c("cfDNA_percent")
vars_abmr <- c("ABMRpm", "ggt0", "ptcgt0", "DSAST", "NKB") # "cggt0", "RejAA_EABMR", "RejAA_FABMR", "RejAA_LABMR")
vars_tcmr <- c("TCMRt", "tgt1", "igt1", "QCAT", "TCB") # , "TCMR-RAT", )

# Endothelium <- c("ENDAT")
vars_parenchyma <- c("KT1", "KT2")
vars_macrophage <- c("AMAT1", "QCMAT")
# Injuryrecent <- c("FICOL", "IRRAT30", "IRITD3", "IRITD5")
vars_atrophyfibrosis <- c("IGT", "MCAT", "BAT", "cigt1", "ctgt1")


# DEFINE VARIABLES TO ASSESS ####
# vars <- c("cfDNA", "ABMRpm", "ggt0", "ptcgt0", "DSAST", "NKB", "TCMRt", "tgt1", "igt1", "QCAT", "TCB")
vars <- c(vars_cfDNA, vars_abmr, vars_tcmr, vars_macrophage, vars_atrophyfibrosis, vars_parenchyma)


# DEFINE THE data ####
data <- data_scores_k1208 %>%
    dplyr::select(Center, Patient, Felzartamab, Group, Followup, Felzartamab_Group, Felzartamab_Followup, all_of(vars)) %>%
    dplyr::filter(Group != "FU1b", Patient %nin% c(15, 18)) %>%
    left_join(., summarise(., sample_pairs = n(), .by = Felzartamab_Followup), by = "Felzartamab_Followup") %>%
    relocate(sample_pairs, .after = Felzartamab_Followup) %>%
    arrange(Felzartamab, Patient, Group)


# WRANGLE THE PHENOTYPE DATA ####
df00 <- data %>%
    expand_grid(category = c("cfDNA", "ABMR", "TCMR", "macrophage", "atrophyfibrosis", "parenchyma")) %>%
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
                }else if (category == "macrophage") {
                    vars_macrophage
                }else if (category == "atrophyfibrosis") {
                    vars_atrophyfibrosis
                }else if (category == "parenchyma") {
                    vars_parenchyma
                }
            }
        ),
        data = pmap(
            list(features, data),
            function(features, data) {
                data %>%
                    dplyr::select(
                        Center, Patient, Felzartamab, Group, Followup, Felzartamab_Group, Felzartamab_Followup,
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
                features <- data %>% dplyr::select(dplyr::all_of(features))
                Followup <- data$Followup
                Felzartamab <- data$Felzartamab
                vegan::adonis2(
                    features ~ Followup * Felzartamab,
                    data = data,
                    method = "euclidean",
                    # by = "margin", # only specify margin if the sample sizes are unequal
                    permutations = 1000000
                )
            }
        ),
        permanova_pairwise = pmap(
            list(data, features),
            function(data, features) {
                set.seed(seed)
                features <- data %>% dplyr::select(dplyr::all_of(features))
                Felzartamab_Followup <- data$Felzartamab_Followup
                pairwiseAdonis::pairwise.adonis(
                    features,
                    Felzartamab_Followup,
                    reduce = "Followup|Felzartamab",
                    sim.method = "euclidean",
                    p.adjust.m = "fdr",
                    perm = 10000
                )
            }
        )
    )
names(df_permanova$permanova) <- df_permanova$category

df_permanova$permanova_pairwise[[2]]
df_permanova$permanova[[2]]

df_permanova$permanova
