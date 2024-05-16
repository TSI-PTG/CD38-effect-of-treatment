# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(ggpubr) # install.packages("ggpubr")
library(patchwork) # install.packages("patchwork")
library(ggrepel) # install.packages("ggrepel")
# Custom operators, functions, and datasets
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
# source plot function
source("C:/R/CD38-effect-of-treatment/code/functions/plot.gg_bland_altman.r")
# load reference data
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/felzartamab_ARTanova.RData")


# DEFINE THE DATA ####
data <- felzartamab_ARTanova %>%
    dplyr::select(category:data)


# MAKE INDIVIDUAL PLOTS ####
bland_altman_plots <- data %>%
    mutate(
        plot_bland_altman = pmap(
            list(data, variable, score),
            gg_bland_altman
        )
    )


# SAVE THE PLOT DATA ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/"
save(bland_altman_plots, file = paste(saveDir, "Felzartamab_bland_altman_plots.RData", sep = ""))

