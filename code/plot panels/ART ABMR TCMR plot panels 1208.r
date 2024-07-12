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



# ABMRpm PATIENT PANELS ####
panel_pairs_ABMRpm <- felzartamab_plots %>%
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
        labels = "B",
        font.label = list(size = 25, face = "bold")
    )  %>%
    ggpubr::annotate_figure(
        top = text_grob("Tracking individual patient responses", face = "bold.italic", size = 25, hjust = 0.73)
    )


# MOLECULAR ABMR PANELS ####
panel_violin_abmr <- felzartamab_plots %>%
    dplyr::filter(category %in% c("ABMR")) %>%
    pull(plot_violin) %>%
    wrap_plots(nrow = 1, ncol = 5) &
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

panel_violin_abmr <- panel_violin_abmr %>%
    ggarrange(
        labels = "C",
        font.label = list(size = 25, face = "bold"),
        legend = "none"
    ) %>%
    ggpubr::annotate_figure(
        top = text_grob("Effect of felzartamab on molecular ABMR activity scores",
            face = "bold.italic",
            size = 25,
            hjust = 1.26
        )
    )


# MOLECULAR TCMR PANELS ####
panel_violin_tcmr <- felzartamab_plots %>%
    dplyr::filter(category %in% c("TCMR")) %>%
    pull(plot_violin) %>%
    wrap_plots(nrow = 1, ncol = 5) &
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

panel_violin_tcmr <- panel_violin_tcmr %>%
    ggarrange(
        labels = "D",
        font.label = list(size = 25, face = "bold"),
        legend = "none"
    ) %>%
    ggpubr::annotate_figure(
        top = text_grob("Effect of felzartamab on molecular TCMR activity scores",
            face = "bold.italic",
            size = 25,
            hjust = 1.26
        )
    )



# COMBINED VIOLIN PANELS ####
panels_violin <- ggarrange(
    panel_violin_abmr,
    panel_violin_tcmr,
    nrow = 2,
    heights = c(1, 1)
)

panels_violin_wlegend <- ggarrange(
    panel_legend,
    panels_violin,
    nrow = 2,
    heights = c(0.125, 1)
)


# COMBINED FINAL PANELS ####
panels <- wrap_plots(
    A = plot_spacer(),
    B = panel_pairs_ABMRpm,
    C = panels_violin_wlegend,
    design = c(
        "AAABB
        CCCCC
        CCCCC"
    ),
    nrow = 2
    # heights =
)



# SAVE THE PLOTS ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/output/"
ggsave(
    filename = paste(saveDir, "Felzartamab violin plots cfdna ABMR and TCMR scores.png"),
    plot = panels,
    dpi = 300,
    width = 60,
    height = 36,
    units = "cm",
    bg = "white"
)

# END ####
