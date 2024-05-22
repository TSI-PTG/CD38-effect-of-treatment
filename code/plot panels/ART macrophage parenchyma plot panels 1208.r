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


# EXTRACT LEGEND FOR PLOTS ####
panel_legend <- felzartamab_plots %>%
    dplyr::filter(variable == "AMAT1") %>%
    pull(plot_violin) %>%
    ggpubr::get_legend() %>%
    ggpubr::as_ggplot() +
    theme(plot.margin = unit(c(0, 0, -1, 0), "cm"))


# MAKE PANEL OF ARCHETYPE PLOTS ####
panel_violin_macrophages <- felzartamab_plots %>%
    dplyr::filter(category %in% c("macrophage")) %>%
    pull(plot_violin) %>%
    wrap_plots(nrow = 1, ncol = 2) &
    theme(
        legend.position = "none",
        axis.text = element_text(size = 10, colour = "black"),
        plot.tag = element_text(size = 20, face = "bold", vjust = 1),
    )

panel_violin_macrophages <- panel_violin_macrophages %>%
    ggarrange(
        labels = "A",
        font.label = list(size = 25, face = "bold"),
        legend = "none"
    ) %>%
    ggpubr::annotate_figure(
        top = text_grob(
            "Effect of Felzartamab on Macrophage Scores",
            face = "bold.italic", size = 25, hjust = 0.912
        )
    )


# MAKE PANEL OF REJECTION PC PLOTS ####
panel_violin_parenchyma <- felzartamab_plots %>%
    dplyr::filter(category %in% c("parenchyma")) %>%
    pull(plot_violin) %>%
    wrap_plots(nrow = 1, ncol = 2) &
    theme(
        legend.position = "none",
        axis.text = element_text(size = 10, colour = "black"),
        plot.tag = element_text(size = 20, face = "bold", vjust = 1),
        plot.background = element_rect(fill = "grey95", colour = " white")
    )

panel_violin_parenchyma <- panel_violin_parenchyma %>%
    ggarrange(
        labels = "B",
        font.label = list(size = 25, face = "bold"),
        legend = "none"
    ) %>%
    ggpubr::annotate_figure(
        #         fig.lab = "fuckyou",   fig.lab.size = 25,
        #   fig.lab.face = "bold.italic"
        top = text_grob(
            "Effect of Felzartamab on Kidney Parenchyma Scores",
            face = "bold.italic", size = 25, hjust = 0.775
        )
    )

panels_violin_legend <- ggarrange(
    panel_legend,
    panel_violin_macrophages,
    panel_violin_parenchyma,
    nrow = 3,
    heights = c(0.125, 0.5, 0.5)
)


# SAVE THE PLOTS ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/output/"
ggsave(
    filename = paste(saveDir, "Felzartamab violin plots macrophage paranchymal scores.png"),
    plot = panels_violin_legend,
    dpi = 300,
    width = 35,
    height = 28,
    units = "cm",
    bg = "white"
)
# ggsave(
#     filename = paste(saveDir, "Felzartamab patient plots ABMR and TCMR scores.png"),
#     plot = panel_patient,
#     dpi = 300,
#     width = 60,
#     height = 22,
#     units = "cm",
#     bg = "white"
# )


# END ####
