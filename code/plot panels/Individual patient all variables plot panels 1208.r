# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(ggpubr) # install.packages("ggpubr")
library(patchwork) # install.packages("patchwork")
library(ggprism) # install.packages("ggprism")
# Custom operators, functions, and datasets
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
# load reference data
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/Felzartamab_plots.RData")



# MAKE PANEL OF PATIENT PAIR PLOTS ####
panel_patient <- felzartamab_plots %>%
    pull(plot_patient_pairs) %>%
    wrap_plots(nrow = 4, ncol = 6) &
    # plot_annotation(tag_levels = list(c(LETTERS[1:15]))) &
    theme(
        legend.position = "none",
        axis.text = element_text(size = 10, colour = "black"), plot.tag = element_text(size = 20, face = "bold", vjust = 1)
    )


# SAVE THE PLOTS ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/output/"
ggsave(
    filename = paste(saveDir, "Felzartamab patient plots all scores.png"),
    plot = panel_patient,
    dpi = 300,
    width = 80,
    height = 40,
    units = "cm",
    bg = "white"
)
# ggsave(
#     filename = paste(saveDir, "Felzartamab violin legend.png"),
#     plot = panel_legend,
#     dpi = 300,
#     width = 22,
#     height = 5,
#     units = "cm",
#     bg = "white"
# )



# END ####
