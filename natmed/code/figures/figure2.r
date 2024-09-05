# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(ggpubr) # install.packages("ggpubr")
library(patchwork) # install.packages("patchwork")
library(ggprism) # install.packages("ggprism")
# Custom operators and functions
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
# load plot data
load("natmed/results/plots_artANOVA.RData")


# EXTRACT LEGEND FOR PLOTS ####
panel_legend <- plots_artANOVA %>%
    dplyr::filter(variable == "AMAT1") %>%
    pull(plot_violin) %>%
    ggpubr::get_legend() %>%
    ggpubr::as_ggplot() +
    theme(plot.margin = unit(c(0, 0, -1, 0), "cm"))


# ABMRpm PATIENT PANELS ####
panel_pairs_ABMRpm <- plots_artANOVA %>%
    dplyr::filter(variable == c("ABMRpm")) %>%
    pull(plot_patient_pairs) %>%
    wrap_plots() +
    theme(
        legend.position = "none",
        axis.text.x = ggplot2::element_text(size = 15, colour = "black"),
        plot.margin = unit(c(0, 0, -3, 0), "cm")
    )

panel_pairs_ABMRpm <- panel_pairs_ABMRpm %>%
    ggarrange(
        labels = "A",
        font.label = list(size = 25, face = "bold")
    ) %>%
    ggpubr::annotate_figure(
        top = text_grob("Tracking individual patient responses", face = "bold.italic", size = 25, hjust = 0.9)
    )


# MOLECULAR ABMR PANELS ####
panel_violin_abmr <- plots_artANOVA %>%
    dplyr::filter(category %in% c("ABMR")) %>%
    pull(plot_violin) %>%
    wrap_plots(nrow = 1, ncol = 5) &
    theme(
        axis.text = element_text(size = 10, colour = "black"),
        legend.position = "none",
        plot.background = element_rect(fill = "grey95", colour = " white")
    )

panel_violin_abmr <- panel_violin_abmr %>%
    ggarrange(
        labels = "B",
        font.label = list(size = 25, face = "bold"),
        legend = "none"
    ) %>%
    ggpubr::annotate_figure(
        top = text_grob("Effect of felzartamab on molecular ABMR activity scores",
            face = "bold.italic",
            size = 25,
            hjust = 1.27
        )
    )


# MOLECULAR TCMR PANELS ####
panel_violin_tcmr <- plots_artANOVA %>%
    dplyr::filter(category %in% c("TCMR")) %>%
    pull(plot_violin) %>%
    wrap_plots(nrow = 1, ncol = 5) &
    theme(
        axis.text = element_text(size = 10, colour = "black"),
        legend.position = "none",
    )

panel_violin_tcmr <- panel_violin_tcmr %>%
    ggarrange(
        labels = "C",
        font.label = list(size = 25, face = "bold"),
        legend = "none"
    ) %>%
    ggpubr::annotate_figure(
        top = text_grob("Effect of felzartamab on molecular TCMR activity scores",
            face = "bold.italic",
            size = 25,
            hjust = 1.27
        )
    )


# COMBINED VIOLIN PANELS ####
panels_violin <- ggarrange(
    panel_violin_abmr,
    panel_violin_tcmr,
    nrow = 2,
    heights = c(1, 1)
)

panels_violin_with_patient_scores <- ggarrange(
    panel_pairs_ABMRpm,
    panels_violin,
    ncol = 2,
    widths = c(5, 10)
)

figure2 <- ggarrange(
    panel_legend,
    panels_violin_with_patient_scores,
    nrow = 2,
    heights = c(0.125, 1)
)


# SAVE THE PLOTS ####
saveDir <- "natmed/png/"
ggsave(
    filename = paste0(saveDir, "figure2.png"),
    plot = figure2,
    dpi = 300,
    width = 90,
    height = 25,
    units = "cm",
    bg = "white"
)

# END ####
