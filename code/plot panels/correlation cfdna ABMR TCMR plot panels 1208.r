# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(ggpubr) # install.packages("ggpubr")
library(patchwork) # install.packages("patchwork")
# Custom operators, functions, and datasets
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
# source plot function
# source("C:/R/CD38-effect-of-treatment/code/functions/plot.gg_violin_interaction.r")
# load reference data
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/Felzartamab_cfdna_percent_quantile_regression_plots_1208.RData")





# MAKE PANEL OF CORRELATION PLOTS ####
panel_cfdna_correlation_1 <- felzartamab_cfdna_qr_plots %>%
    dplyr::filter(
        Followup_pairwise %in% c("Baseline - Week24"),
        category %in% c("ABMR")
    ) %>%
    pull(plot_scatter) %>%
    wrap_plots(nrow = 1, ncol = 5) +
    plot_annotation(
        title = "Baseline - Week24",
        tag_levels = list(c("A", rep("", 9)))
    ) &
    theme(
        legend.position = "none",
        plot.title = element_text(size = 25, face = "bold.italic"),
        axis.text = element_text(size = 10, colour = "black"),
        plot.tag = element_text(size = 20, face = "bold", vjust = 1)
    )

panel_cfdna_correlation_2 <- felzartamab_cfdna_qr_plots %>%
    dplyr::filter(
        Followup_pairwise %in% c("Baseline - Week52"),
        category %in% c("ABMR")
    ) %>%
    pull(plot_scatter) %>%
    wrap_plots(nrow = 1, ncol = 5) +
    plot_annotation(
        title = "Baseline - Week52",
        tag_levels = list(c("B", rep("", 9)))
    ) &
    theme(
        legend.position = "none",
        plot.title = element_text(size = 25, face = "bold.italic"),
        axis.text = element_text(size = 10, colour = "black"),
        plot.tag = element_text(size = 20, face = "bold", vjust = 1)
    )

panel_cfdna_correlation_3 <- felzartamab_cfdna_qr_plots %>%
    dplyr::filter(
        Followup_pairwise %in% c("Week24 - Week52"),
        category %in% c("ABMR")
    ) %>%
    pull(plot_scatter) %>%
    wrap_plots(nrow = 1, ncol = 5) +
    plot_annotation(
        title = "Week24 - Week52",
        tag_levels = list(c("C", rep("", 9)))
    ) &
    theme(
        legend.position = "none",
        plot.title = element_text(size = 25, face = "bold.italic"),
        axis.text = element_text(size = 10, colour = "black"),
        plot.tag = element_text(size = 20, face = "bold", vjust = 1),
        plot.background = element_rect(fill = "grey95", colour = " white")
    )

panels <- ggarrange(
    panel_cfdna_correlation_1,
    panel_cfdna_correlation_3,
    panel_cfdna_correlation_2,
    nrow = 3
)


# SAVE THE PLOTS ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/output/"
ggsave(
    filename = paste(saveDir, "Felzartamab cfdna percent correlation plots cfDNA.png"),
    plot = panels,
    dpi = 300,
    width = 60,
    height = 30,
    units = "cm",
    bg = "white"
)
