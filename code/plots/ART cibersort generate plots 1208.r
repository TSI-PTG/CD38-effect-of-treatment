# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(ggbeeswarm) # install.packages("ggbeeswarm")
library(ggpubr) # install.packages("ggpubr")
library(patchwork) # install.packages("patchwork")
library(ggrepel) # install.packages("ggrepel")
library(gghalves) # pak::pak("erocoar/gghalves")
library(ggh4x) # install.packages("ggh4x")
library(ggprism) # install.packages("ggprism")
library(rstatix) # install.packages("rstatix")
# Custom operators, functions, and datasets
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
se <- function(x) sd(x) / sqrt(length((x)))
log10zero <- scales::trans_new(
    name = "log10zero",
    transform = function(x) log10(x + 0.001),
    inverse = function(x) 10^x - 0.001
)
# source plot function
source("C:/R/CD38-effect-of-treatment/code/functions/plot.gg_violin_interaction_means.r")
source("C:/R/CD38-effect-of-treatment/code/functions/plot.gg_patient_pairs_interaction.r")
# load reference data
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/felzartamab_ARTanova_cibersort.RData")



# MAKE INDIVIDUAL PLOTS ####
felzartamab_plots <- felzartamab_ARTanova %>%
    mutate(
        plot_violin = pmap(
            list(data, variable, score, means_delta, art_con_interaction_default_tidy),
            gg_violin_interaction, patient_label = c(4, 9, 13)
        ),
        plot_patient_pairs = pmap(
            list(data, variable, score),
            gg_patient_pairs_interaction
        )
    )

felzartamab_plots$plot_violin[[8]] # reference legend1
felzartamab_plots$plot_violin[[10]] # reference legend2

felzartamab_plots$plot_patient_pairs[[10]]


# SAVE THE PLOT DATA ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/"
save(felzartamab_plots, file = paste(saveDir, "Felzartamab_plots_cibersort.RData", sep = ""))
