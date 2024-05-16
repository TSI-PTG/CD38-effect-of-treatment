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
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/Felzartamab_bland_altman_plots.RData")



# MAKE PANEL OF BLAND ALTMAN PLOTS ####
panels_bland_altman <- bland_altman_plots %>%
    dplyr::filter(category %in% c("ABMR", "TCMR")) %>%
    pull(plot_bland_altman) %>%
    wrap_plots(nrow = 2, ncol = 5) +
    plot_annotation(
        title = "Week24 - Day0",
        tag_levels = list(c(LETTERS[1:15]))) &
    theme(
        legend.position = "none",
        plot.title = element_text(size = 20, face = "bold"),
        axis.text = element_text(size = 10, colour = "black"), plot.tag = element_text(size = 20, face = "bold", vjust = 1)
    )


# SAVE THE PLOTS ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/output/"
ggsave(
    filename = paste(saveDir, "Felzartamab Bland Altman plots ABMR and TCMR scores.png"),
    plot = panels_bland_altman,
    dpi = 300,
    width = 60,
    height = 22,
    units = "cm",
    bg = "white"
)
