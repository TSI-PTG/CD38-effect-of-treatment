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



# CFDNA SCORE PANELS ####
panel_violin_cfdna <- felzartamab_plots %>%
    dplyr::filter(category %in% c("cfDNA")) %>%
    pull(plot_violin) %>%
    wrap_plots() +
    theme(
        legend.position = "none",
        axis.text.x = ggplot2::element_text(size = 10, colour = "black")
    )

panel_violin_cfdna <- panel_violin_cfdna %>%
    ggarrange(
        labels = "A",
        font.label = list(size = 25, face = "bold")
        )  %>% 
    ggpubr::annotate_figure(
        top = text_grob("Effect of Felzartamab on dd-cfDNA", face = "bold.italic", size = 25, hjust = 0.55)
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
        labels = "B",
        font.label = list(size = 25, face = "bold"), 
        legend = "none"
    ) %>%
    ggpubr::annotate_figure(
        top = text_grob("Effect of Felzartamab on molecular ABMR activity", face = "bold.italic", size = 25, hjust = 1.55)
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
        labels = "C",
        font.label = list(size = 25, face = "bold"),
        legend = "none"
    ) %>%
    ggpubr::annotate_figure(
        top = text_grob("Effect of Felzartamab on molecular TCMR activity", face = "bold.italic", size = 25, hjust = 1.55)
    )

panels_violin_scores <- ggarrange(
    panel_violin_abmr,
    panel_violin_tcmr,
    nrow = 2,
    heights = c(1, 1)
)

panels_violin <- ggarrange(
    panel_violin_cfdna,
    panels_violin_scores,
    nrow = 1,
    widths = c(0.25, 1)
)

panels_violin_wlegend <- ggarrange(
    panel_legend,
    panels_violin,
    nrow = 2,
    heights = c(0.125, 1)
)




# panels_violin <- wrap_plots(
#     panel_violin_cfdna,
#     panel_violin_scores
# ) +
#     plot_layout(
#         widths = c(0.25, 1)
#     ) +
#     plot_annotation(
#         tag_level = list(c(
#             "A",
#             "B", rep("", 4),
#             "C", rep("", 4)
#         ))
#     ) &
#     theme(
#         legend.position = "none",
#         axis.text = element_text(size = 10, colour = "black"),
#         plot.title = element_text(size = 12, face = "bold.italic"),
#         plot.tag = element_text(size = 25, face = "bold", vjust = 1)
#     )

# panels_violin_wlegend <- ggarrange(
#     panel_legend,
#     panels_violin,
#     nrow = 2,
#     heights = c(0.1, 1)
# )





# # MAKE PANEL OF PATIENT PAIR PLOTS ####
# panel_patient <- felzartamab_plots %>%
#     dplyr::filter(category %in% c("ABMR", "TCMR")) %>%
#     pull(plot_patient_pairs) %>%
#     wrap_plots(nrow = 2, ncol = 5) +
#     plot_annotation(tag_levels = list(c(LETTERS[1:15]))) &
#     theme(
#         legend.position = "none",
#         axis.text = element_text(size = 10, colour = "black"), plot.tag = element_text(size = 20, face = "bold", vjust = 1)
#     )


# SAVE THE PLOTS ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/output/"
ggsave(
    filename = paste(saveDir, "Felzartamab violin plots cfdna and score.png"),
    plot = panels_violin_wlegend,
    dpi = 300,
    width = 80,
    height = 22,
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

# ggsave(
#     filename = paste(saveDir, "Felzartamab violin legend.png"),
#     plot = panel_legend,
#     dpi = 300,
#     width = 22,
#     height = 5,
#     units = "cm",
#     bg = "white"
# )



# END ####
