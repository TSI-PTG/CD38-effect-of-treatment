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
panel_violin_archetypes <- felzartamab_plots %>%
    dplyr::filter(
        variable %>% str_detect("RejAA_NR|RejAA_FABMR|RejAA_TCMR1")
    ) %>%
    pull(plot_violin) %>%
    wrap_plots(nrow = 1, ncol = 3) &
    # plot_annotation(tag_levels = list(c("A", "", ""))) &
    theme(
        legend.position = "none",
        axis.text = element_text(size = 10, colour = "black"),
        plot.tag = element_text(size = 20, face = "bold", vjust = 1),
    )

panel_violin_archetypes <- panel_violin_archetypes %>%
    ggarrange(
        labels = "A",
        font.label = list(size = 25, face = "bold"),
        legend = "none"
    ) %>%
    ggpubr::annotate_figure(
        top = text_grob(
            "Effect of felzartamab on rejection archetype scores",
            face = "bold.italic", size = 25, hjust = 0.919
        )
    )


# MAKE PANEL OF REJECTION PC PLOTS ####
panel_violin_RejPC <- felzartamab_plots %>%
    dplyr::filter(
        variable %>% str_detect("RejPC")
    ) %>%
    dplyr::filter() %>%
    pull(plot_violin) %>%
    wrap_plots(nrow = 1, ncol = 3) &
    # plot_annotation(tag_levels = list(c("B", "", ""))) &
    theme(
        legend.position = "none",
        axis.text = element_text(size = 10, colour = "black"),
        plot.tag = element_text(size = 20, face = "bold", vjust = 1),
        plot.background = element_rect(fill = "grey95", colour = " white")
    )

panel_violin_RejPC <- panel_violin_RejPC %>%
    ggarrange(
        labels = "B",
        font.label = list(size = 25, face = "bold"),
        legend = "none"
    ) %>%
    ggpubr::annotate_figure(
        #         fig.lab = "fuckyou",   fig.lab.size = 25,
        #   fig.lab.face = "bold.italic"
        top = text_grob(
            "Effect of felzartamab on rejection principal component scores",
            face = "bold.italic", size = 25, hjust = 0.762
        )
    )


# MAKE PANEL OF INJURY PC PLOTS ####
panel_violin_InjPC <- felzartamab_plots %>%
    dplyr::filter(
        variable %>% str_detect("InjPC")
    ) %>%
    dplyr::filter() %>%
    pull(plot_violin) %>%
    wrap_plots(nrow = 1, ncol = 3) &
    # plot_annotation(
    #     title = "InjPC",
    #     tag_levels = list(c("C", "", ""))
    # ) &
    theme(
        legend.position = "none",
        # plot.title = element_text(size = 25, face = "bold.italic"),
        axis.text = element_text(size = 10, colour = "black"),
        plot.tag = element_text(size = 20, face = "bold", vjust = 1),
    )

panel_violin_InjPC <- panel_violin_InjPC %>%
    ggarrange(
        labels = "C",
        font.label = list(size = 25, face = "bold"),
        legend = "none"
    ) %>%
    ggpubr::annotate_figure(
        top = text_grob(
            "Effect of felzartamab on injury principal component scores",
            face = "bold.italic", size = 25, hjust = 0.80
        )
    )




panels_violin_legend <- ggarrange(
    panel_legend,
    panel_violin_archetypes,
    panel_violin_RejPC,
    panel_violin_InjPC,
    nrow = 4,
    heights = c(0.125, 0.5, 0.5, 0.5)
)


# SAVE THE PLOTS ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/output/"
ggsave(
    filename = paste(saveDir, "Felzartamab violin plots archetype PC scores.png"),
    plot = panels_violin_legend,
    dpi = 300,
    width = 40,
    height = 33,
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
