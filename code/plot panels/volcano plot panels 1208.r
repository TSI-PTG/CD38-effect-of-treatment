# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(ggrepel) # install.packages("ggrepel")
library(ggpubr) # install.packages("ggpubr")
library(patchwork) # install.packages("patchwork")
library(readxl) # install.packages("readxl")
# Custom operators, functions, and datasets
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
# load plot data
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/Volcano plots IQR_filtered_probes_unique_genes_baseline_corrected_cortex_corrected 1208.RData")


# MAKE THE PLOT PANELS ####
plot_volcano_panel <- wrap_plots(
    plot_volcano$plot_volcano,
    nrow = 1, ncol = 3
) +
    plot_annotation(tag_levels = "A") +
    plot_layout(
        guides = "collect",
        axes = "collect", axis_titles = "collect"
    ) &
    theme(
        legend.position = "top",
        plot.background = element_rect(fill = "white"),
        plot.tag = element_text(vjust = 5)
    )


# SAVE THE PLOT ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/output/"
ggsave(
    plot_volcano_panel,
    file = paste(saveDir, "Volcano IQR_filtered_probes_unique_genes_baseline_corrected_cortex_corrected 1208.png", sep = ""),
    dpi = 300,
    width = 30,
    height = 12,
    units = "cm",
    bg = "white"
) %>% suppressWarnings()
