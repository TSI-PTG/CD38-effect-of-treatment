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
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/Volcano timeseries plots IQR_filtered_probes_unique_genes_baseline_corrected_cortex_corrected 1208.RData")
# load dotplot data
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/gene_DEG_plots.RData")



# UNIUVERSAL PLOTTING PARAMETERS ####
size_title <- 15



# MAKE PLOTS PANELS ####
df_DEG_plot <- DEG_plots %>%
    dplyr::filter(geneset %in% c(
        "IFNG-inducible ABMR activity genes",
        "NK cell expressed ABMR activity genes",
        "ABMR-associated endothelial genes"
    ))


# MAKE THE DOTPLOT PANELS ####
panel_DEG <- wrap_plots(
    c(
        df_DEG_plot$plot_long
    ),
    nrow = 3, ncol = 1
) +
    plot_annotation(tag_levels = list(c("E", "F", "G"))) +
    plot_layout(
        guides = "collect",
        # axes = "collect", axis_titles = "collect"
    ) &
    theme(
        axis.text.y = element_text(size = 5),
        legend.position = "top",
        plot.background = element_rect(fill = "white"),
        plot.tag = element_text(vjust = 4)
    )

panel_DEG <- panel_DEG %>%
    ggarrange(., NULL, ncol = 2, nrow = 1, widths = c(1, 0)) %>%
    ggpubr::annotate_figure(
        top = text_grob(
            "Effect Felzartamab on Select Genes",
            face = "bold.italic", size = size_title, hjust = 0.55
        )
    )
    

# MAKE THE VOLCANO PANELS ####
panel_volcano <- wrap_plots(
    A = plot_volcano$plot_volcano[[1]],
    B = plot_volcano$plot_volcano[[2]],
    C = plot_volcano$plot_volcano[[3]],
    design = "ABC"
) +
    plot_annotation(tag_levels = "A") +
    plot_layout(
        # guides = "collect"
        # axes = "collect", axis_titles = "collect"
    ) &
    theme(
        legend.position = "top",
        plot.background = element_rect(fill = "white"),
        plot.tag = element_text(vjust = 4)
    )

panel_volcano <- panel_volcano %>%
    ggarrange(., NULL, ncol = 2, nrow = 1, widths = c(1, 0)) %>%
    ggpubr::annotate_figure(
        top = text_grob(
            "Genome-wide Effect Felzartamab",
            face = "bold.italic", size = size_title, hjust = 1.28
        )
    )


# MAKE THE VOLCANO PANELS ####
panel_volcano_timeseries <- plot_volcano_timeseries +
    plot_annotation(tag_levels = list(c("D"))) &
    theme(
        legend.position = "top",
        plot.background = element_rect(fill = "grey95"),
        plot.tag = element_text(vjust = 4)
    )

panel_volcano_timeseries <- panel_volcano_timeseries %>%
    ggarrange(NULL, ., NULL, ncol = 3, nrow = 1, widths = c(0.25, 1, 0.25)) %>%
    ggpubr::annotate_figure(
        top = text_grob(
            "Longitiduinal Relapse of Effect Felzartamab",
            face = "bold.italic", size = size_title, hjust = 0.635 #hjust = 0.95
        )
    )


panels <- ggarrange(
    panel_volcano,
    panel_volcano_timeseries,
    nrow = 2, ncol = 1
) %>%ggarrange(
    panel_DEG, 
    nrow = 1, ncol = 2,
    widths = c(1, 0.5)
)
    





# SAVE THE PLOT ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/output/"
ggsave(
    panels,
    file = paste(saveDir, "Volcano timerseries GEP IQR_filtered_probes_unique_genes_baseline_corrected_cortex_corrected 1208.png", sep = ""),
    dpi = 300,
    width = 32,
    height = 18,
    units = "cm",
    bg = "white"
) %>% suppressWarnings()

# ggsave(
#     plot_panel,
#     file = paste(saveDir, "Volcano GEP IQR_filtered_probes_unique_genes_baseline_corrected_cortex_corrected 1208.png", sep = ""),
#     dpi = 300,
#     width = 32,
#     height = 18,
#     units = "cm",
#     bg = "white"
# ) %>% suppressWarnings()
