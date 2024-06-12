# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(patchwork) # install.packages("patchwork")
# Custom operators, functions, and datasets
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
# load reference data
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/gene_DEG_plots.RData")


# MAKE PLOTS PANELS ####
plot_panels_wide <- DEG_plots %>%
    dplyr::filter(geneset != "ABMR_activity") %>%
    pull(plot_wide) %>%
    wrap_plots() +
    plot_layout(axes = "collect") +

    plot_annotation(tag_levels = "A")


plot_panels_long <- DEG_plots %>%
    dplyr::filter(geneset != "ABMR_activity") %>%
    pull(plot_long) %>%
    wrap_plots() +
    # plot_layout(axes = "collect") +
    plot_annotation(tag_levels = "A")


# SAVE THE PLOTS ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/output/"
ggsave(
    filename = paste(saveDir, "Selective DEG panel wide.png"),
    plot = plot_panels_wide,
    dpi = 600,
    width = 60,
    height = 8,
    units = "cm",
    bg = "white"
)
ggsave(
    filename = paste(saveDir, "Selective DEG panel long.png"),
    plot = plot_panels_long,
    dpi = 600,
    width = 40,
    height = 8,
    units = "cm",
    bg = "white"
)
