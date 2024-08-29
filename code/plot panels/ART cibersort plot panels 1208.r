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
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/Felzartamab_plots_cibersort.RData")


felzartamab_plots$category
felzartamab_plots <- felzartamab_plots %>%
    arrange(category, score)


# MOLECULAR ABMR PANELS ####
panel_violin_nk <- felzartamab_plots %>%
    dplyr::filter(category %in% c("NK cells")) %>%
    pull(plot_violin) %>%
    wrap_plots(nrow = 1, ncol = 2) &
    # plot_annotation(
    #     tag_level = list(c(
    #         "B", rep("", 5)
    #     ))
    # ) &
    theme(
        axis.text = element_text(size = 10, colour = "black"),
        legend.position = "none",
        # plot.tag = element_text(size = 25, face = "bold", vjust = 1, hjust = 1),
        plot.background = element_rect(fill = "grey95", colour = " white")
    )

panel_violin_nk <- panel_violin_nk %>%
    ggarrange(
        labels = "A",
        font.label = list(size = 25, face = "bold"),
        legend = "none"
    ) %>%
    ggpubr::annotate_figure(
        top = text_grob("NK cell fractions",
            face = "bold.italic",
            size = 25,
            hjust = 1.5
        )
    )


# MOLECULAR ABMR PANELS ####
panel_violin_monomacro <- felzartamab_plots %>%
    dplyr::filter(category %in% c("Monocytes", "Macrophages")) %>%
    pull(plot_violin) %>%
    wrap_plots(nrow = 1, ncol = 4) &
    # plot_annotation(
    #     tag_level = list(c(
    #         "B", rep("", 5)
    #     ))
    # ) &
    theme(
        axis.text = element_text(size = 10, colour = "black"),
        legend.position = "none",
        # plot.tag = element_text(size = 25, face = "bold", vjust = 1, hjust = 1),
        plot.background = element_rect(fill = "grey95", colour = " white")
    )

panel_violin_monomacro <- panel_violin_monomacro %>%
    ggarrange(
        labels = "B",
        font.label = list(size = 25, face = "bold"),
        legend = "none"
    ) %>%
    ggpubr::annotate_figure(
        top = text_grob("Monocyte and Macrophage cell fractions",
            face = "bold.italic",
            size = 25,
            hjust = 1.3
        )
    )


# MOLECULAR TCMR PANELS ####
panel_violin_tcell <- felzartamab_plots %>%
    dplyr::filter(category %in% c("T cells")) %>%
    pull(plot_violin) %>%
    wrap_plots(nrow = 1, ncol = 7) &
    # plot_annotation(
    #     tag_level = list(c(
    #         "B", rep("", 5)
    #     ))
    # ) &
    theme(
        axis.text = element_text(size = 10, colour = "black"),
        legend.position = "none",
        # plot.tag = element_text(size = 25, face = "bold", vjust = 1, hjust = 1),
    )

panel_violin_tcell <- panel_violin_tcell %>%
    ggarrange(
        labels = "C",
        font.label = list(size = 25, face = "bold"),
        legend = "none"
    ) %>%
    ggpubr::annotate_figure(
        top = text_grob("T cell fractions",
            face = "bold.italic",
            size = 25,
            hjust = 5.5
        )
    )




# COMBINED VIOLIN PANELS ####
panels_violin1 <- ggarrange(
    panel_violin_nk,
    panel_violin_monomacro,
    nrow = 1,
    widths = c(0.5, 1)
)

panels_violin2 <- ggarrange(
    panels_violin1,
    panel_violin_tcell,
    nrow = 2,
    heights = c(1, 1)
)

panels_violin_wlegend <- ggarrange(
    panel_legend,
    panels_violin2,
    nrow = 2,
    heights = c(0.125, 1)
)


# # COMBINED FINAL PANELS ####
# panels <- wrap_plots(
#     # A = plot_spacer(),
#     B = panel_pairs_ABMRpm,
#     C = panels_violin_wlegend,
#     design = c(
#         "A
#         B
#         C
#         C"
#     ),
#     heights = c(4, 3, 2)
# )


# SAVE THE PLOTS ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/output/"
ggsave(
    filename = paste(saveDir, "Felzartamab violin plots cibersort.png"),
    plot = panels_violin_wlegend,
    dpi = 300,
    width = 70,
    height = 25,
    units = "cm",
    bg = "white"
)

# END ####
