# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(ggbeeswarm) # install.packages("ggbeeswarm")
library(ggpubr) # install.packages("ggpubr")
library(patchwork) # install.packages("patchwork")
library(ggrepel) # install.packages("ggrepel")
library(gghalves) # pak::pak("erocoar/gghalves") devtools::install_github("erocoar/gghalves")
library(ggh4x) # install.packages("ggh4x")
library(ggprism) # install.packages("ggprism")
library(rstatix) # install.packages("rstatix")
library(ggsci) # install.packages("ggsci")
# Custom operators and functions
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
# se <- function(x) sd(x) / sqrt(length((x)))
# log10zero <- scales::trans_new(
#     name = "log10zero",
#     transform = function(x) log10(x + 0.001),
#     inverse = function(x) 10^x - 0.001
# )
# source plot function
source("natmed/code/functions/plot.gg_violin_interaction.r")
source("natmed/code/functions/plot.gg_patient_pairs_interaction.r")
# load results
load("natmed/results/artANOVA.RData")


# UNIVERSAL PLOTTING PARAMETERS ####
# scales::show_col(pal_npg()(9))
col_low <- pal_npg()(9)[1]
col_high <- pal_npg()(9)[3]


# MAKE INDIVIDUAL PLOTS ####
plots_artANOVA <- artANOVA %>%
    mutate(
        plot_violin = pmap(
            list(data, variable, score, medians_delta, art_aov_contrast_tidy),
            gg_violin_interaction,
            patient_label = c(4, 9, 13),
            col_low = col_low, col_high = col_high
        )  %>% suppressWarnings(),
        plot_patient_pairs = pmap(
            list(data, variable, score),
            gg_patient_pairs_interaction
        )%>% suppressWarnings()
    ) %>%
    dplyr::select(category:variable, plot_violin, plot_patient_pairs)


# SAVE THE PLOT DATA ####
saveDir <- "natmed/results/"
save(plots_artANOVA, file = paste(saveDir, "plots_artANOVA.RData", sep = ""))
