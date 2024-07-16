# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(ggrepel) # install.packages("ggrepel")
library(ggpubr) # install.packages("ggpubr")
library(patchwork) # install.packages("patchwork")
library(readxl) # install.packages("readxl")
# Custom operators, functions, and datasets
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
# load volcano plot data
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/Volcano plots IQR_filtered_probes_unique_genes_baseline_corrected_cortex_corrected 1208.RData")
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/Volcano enrichment plots IQR_filtered_probes_unique_genes_baseline_corrected_cortex_corrected 1208.RData")
# load time series plot data
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/Volcano timeseries plots IQR_filtered_probes_unique_genes_baseline_corrected_cortex_corrected 1208.RData")
# load dotplot data
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/gene_DEG_plots.RData")



# UNIUVERSAL PLOTTING PARAMETERS ####
size_title <- 15


# MAKE PLOTS PANELS ####
df_DEG_plot <- DEG_plots %>%
    dplyr::filter(geneset %in% c(
        "ABMR activity genes",
        "IFNG-inducible ABMR activity genes",
        "NK cell-expressed ABMR activity genes",
        "ABMR-associated endothelial genes"
    ))



# MAKE THE DOTPLOT PANELS ####
panel_DEG2 <- wrap_plots(
    df_DEG_plot$plot_long,
    nrow = 1, ncol = 4
) +
    plot_annotation(tag_levels = list(c("D", "E", "F", "G"))) +
    plot_layout(
        guides = "collect",
        # axes = "collect", axis_titles = "collect"
    ) &
    theme(
        axis.text.y = element_text(size = 6),
        legend.position = "top",
        plot.background = element_rect(fill = "white"),
        plot.tag = element_text(vjust = 1)
    )

panel_DEG2 <- panel_DEG2 %>%
    ggarrange(., NULL, ncol = 2, nrow = 1, widths = c(1, 0))
    # %>%
    # ggpubr::annotate_figure(
    #     top = text_grob(
    #         "Effect felzartamab on ABMR activity genes",
    #         face = "bold.italic", size = size_title, hjust = 0.8
    #     )
    # )


# MAKE THE VOLCANO TIMESERIES PANELS ####
panel_volcano_timeseries <- plot_volcano_timeseries +
    plot_annotation(tag_levels = list(c("C"))) &
    theme(
        legend.position = "none",
        plot.background = element_rect(fill = "white"),
        plot.tag = element_text(vjust = 1)
    )

# panel_volcano_timeseries1 <- panel_volcano_timeseries %>%
#     ggarrange(NULL, ., NULL, ncol = 3, nrow = 1, widths = c(0.25, 1, 0.25)) %>%
#     ggpubr::annotate_figure(
#         top = text_grob(
#             "Longitudinal relapse of felzartamab effects",
#             face = "bold.italic", size = size_title, hjust = 0.63 # hjust = 0.95
#         )
#     )

panel_volcano_timeseries2 <- panel_volcano_timeseries %>%
    ggarrange(NULL, ., NULL, ncol = 3, nrow = 1, widths = c(0, 1, 0))
# %>%
# ggpubr::annotate_figure(
#     top = text_grob(
#         "Longitudinal relapse of felzartamab effects",
#         face = "bold.italic", size = size_title, hjust = 0.52 # hjust = 0.95
#     )
# )



# MAKE THE VOLCANO ENRICHMENT PANELS ####
panel_volcano_enrichment <- wrap_plots(
    A = plot_volcano_enrichment$plot_volcano_enrichment[[1]],
    B = plot_volcano_enrichment$plot_volcano_enrichment[[2]],
    # C = panel_volcano_timeseries,
    design = "AB"
) +
    plot_annotation(tag_levels = list(c("A", "B"))) +
    plot_layout(
        # guides = "collect"
        # axes = "collect", axis_titles = "collect"
    ) &
    theme(
        legend.position = "none",
        plot.background = element_rect(fill = "white", colour = "white"),
        plot.tag = element_text(vjust = 1, face = "plain")
    )

panel_volcano_enrichment <- panel_volcano_enrichment %>%
    ggarrange(., NULL, ncol = 2, nrow = 1, widths = c(1, 0)) 
    # %>%
    # ggpubr::annotate_figure(
    #     top = text_grob(
    #         "Genome-wide effect of felzartamab treatment",
    #         face = "bold.italic", size = size_title, hjust = 1.3
    #     )
    # )




# ARRANGE PANELS ####
panels_enrichment <- ggarrange(
    ggarrange(
        panel_volcano_enrichment,
        panel_volcano_timeseries2,
        nrow = 1, widths = c(1, 0.8)
    ),
    panel_DEG2,
    nrow = 2, ncol = 1,
    heights = c(1, 0.8)
)




# SAVE THE PLOT ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/output/"
ggsave(
    panels_enrichment,
    file = paste(saveDir, "Volcano enrichment timerseries GEP IQR_filtered_probes_unique_genes_baseline_corrected_cortex_corrected 1208.png", sep = ""),
    dpi = 300,
    width = 35,
    height = 20,
    units = "cm",
    bg = "white"
) %>% suppressWarnings()
