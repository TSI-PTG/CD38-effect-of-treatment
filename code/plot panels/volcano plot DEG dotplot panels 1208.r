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
# load dotplot data
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/gene_DEG_plots.RData")


# MAKE PLOTS PANELS ####
df_DEG_plot <- DEG_plots %>%
    dplyr::filter(geneset %in% c(
        "IFNG-inducible ABMR activity genes",
        "NK cell expressed ABMR activity genes",
        "ABMR-associated endothelial genes"
    )) 
    



# MAKE THE PLOT PANELS ####
plot_panel <- wrap_plots(
    c(plot_volcano$plot_volcano,
    df_DEG_plot$plot_long),
    nrow = 2, ncol = 3,
    heights = c(1, 0.6)
) +
    plot_annotation(tag_levels = "A") +
    plot_layout(
        guides = "collect",
        # axes = "collect", axis_titles = "collect"
    ) &
    theme(
        legend.position = "top",
        plot.background = element_rect(fill = "white"),
        plot.tag = element_text(vjust = 4)
    )


# plot_volcano_panel <- wrap_plots(
#     plot_volcano$plot_volcano,
#     nrow = 1, ncol = 3
# ) +
#     plot_annotation(tag_levels = "A") +
#     plot_layout(
#         guides = "collect",
#         axes = "collect", axis_titles = "collect"
#     ) &
#     theme(
#         legend.position = "top",
#         plot.background = element_rect(fill = "white"),
#         plot.tag = element_text(vjust = 6)
#     )


# SAVE THE PLOT ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/output/"
ggsave(
    plot_panel,
    file = paste(saveDir, "Volcano GEP IQR_filtered_probes_unique_genes_baseline_corrected_cortex_corrected 1208.png", sep = ""),
    dpi = 300,
    width = 32,
    height = 18,
    units = "cm",
    bg = "white"
) %>% suppressWarnings()
