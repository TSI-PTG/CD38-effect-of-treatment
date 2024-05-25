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
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/Felzartamab_cfdna_quantile_regression_plots_1208.RData")




# EXTRACT LEGEND FOR PLOTS ####
legend <- felzartamab_plots %>%
    dplyr::filter(variable == "AMAT1") %>%
    pull(plot_violin) %>%
    pluck(1) +
    guides(
        fill = guide_colorbar(
            title.position = "top",
            barwidth = 12,
            draw.ulim = FALSE,
            draw.llim = FALSE,
            label.hjust = c(-0.075, 1.5), # first value is improved, second value worsened (i.e., reverse = TRUE)
            label.vjust = 7.75,
            reverse = TRUE
        )
    ) +
    theme(legend.title = ggplot2::element_text(size = 15, hjust = 0.2, vjust = 1, face = "bold"))

panel_legend <- legend %>% 
    ggpubr::get_legend() %>%
    ggpubr::as_ggplot()


# MAKE PANEL OF VIOLIN PLOTS ####
plot_cfdna <- felzartamab_plots %>%
    dplyr::filter(category %in% c("cfDNA")) %>%
    pull(plot_violin) %>%
    pluck(1) +
    facet_wrap(~Felzartamab, nrow = 2, ncol = 1, scales = "free_x") +
    theme(
        legend.position = "none",
        axis.text = element_text(size = 12, colour = "black"),
        plot.tag = element_text(size = 20, face = "bold.italic", vjust = 1)
    )

panel_violin <- plot_cfdna + inset_element(
    panel_legend, 
    ignore_tag = TRUE,
    left = 0.085, top = 0.25, right = 1, bottom = 1
    ) +
    plot_layout(ncol = 1) +
    plot_annotation(
        tag_levels = list(c(LETTERS[1:15]))
    ) &
    theme(
        # legend.position = "none",
        axis.text = element_text(size = 12, colour = "black"),
        plot.tag = element_text(size = 20, face = "bold.italic", vjust = 1),
        strip.text = element_text(size = 20)

    )


panel_violin_legend <- panel_violin  %>%
    ggarrange(.,NULL, ncol = 2, nrow = 1, widths = c(1,0)) %>% 
    ggpubr::annotate_figure(
        top = text_grob(
            "Treatment Effect on dd-cfDNA",
            face = "bold.italic", size = 25, hjust = 0.575
        )
    )



# MAKE PANEL OF QUANTILE REGRESSION PLOTS ####
panel_cfdna_correlation_1 <- felzartamab_cfdna_qr_plots %>%
    dplyr::filter(
        Followup_pairwise %in% c("Baseline - Week24"),
        category %in% c("ABMR")
    ) %>%
    pull(plot_scatter) %>%
    wrap_plots(nrow = 1, ncol = 5) +
    plot_annotation(
        title = "B\nBaseline - Week24",
        tag_levels = list(c("", rep("", 9)))
    ) &
    theme(
        legend.position = "none",
        plot.title = element_text(size = 20, face = "bold.italic"),
        axis.text = element_text(size = 10, colour = "black"),
        plot.tag = element_text(size = 20, face = "bold", vjust = 1),
        plot.background = element_rect(fill = "grey95", colour = " white")
    )

panel_cfdna_correlation_2 <- felzartamab_cfdna_qr_plots %>%
    dplyr::filter(
        Followup_pairwise %in% c("Week24 - Week52"),
        category %in% c("ABMR")
    ) %>%
    pull(plot_scatter) %>%
    wrap_plots(nrow = 1, ncol = 5) +
    plot_annotation(
        title = "C\nWeek24 - Week52",
        tag_levels = list(c("", rep("", 9)))
    ) &
    theme(
        legend.position = "none",
        plot.title = element_text(size = 20, face = "bold.italic"),
        axis.text = element_text(size = 10, colour = "black"),
        plot.tag = element_text(size = 20, face = "bold", vjust = 1),
        plot.background = element_rect(fill = "white", colour = " white")
    )


panel_cfdna_correlation_3 <- felzartamab_cfdna_qr_plots %>%
    dplyr::filter(
        Followup_pairwise %in% c("Baseline - Week52"),
        category %in% c("ABMR")
    ) %>%
    pull(plot_scatter) %>%
    wrap_plots(nrow = 1, ncol = 5) +
    plot_annotation(
        title = "D\nBaseline - Week52",
        tag_levels = list(c("", rep("", 9)))
    ) &
    theme(
        legend.position = "none",
        plot.title = element_text(size = 20, face = "bold.italic"),
        axis.text = element_text(size = 10, colour = "black"),
        plot.tag = element_text(size = 20, face = "bold", vjust = 1),
        plot.background = element_rect(fill = "grey95", colour = " white")
    )


panels_quantile_regression <- ggarrange(
    panel_cfdna_correlation_1,
    panel_cfdna_correlation_2,
    panel_cfdna_correlation_3,
    nrow = 3
) %>%
    ggpubr::annotate_figure(
        top = text_grob(
            "Relationship of \u394 dd-cfDNA to \u394 ABMR Activity Scores",
            face = "bold.italic", size = 25, hjust = 1.3085
        )
    )

panels <- ggarrange(
    panel_violin_legend,
    panels_quantile_regression,
    nrow = 1, ncol = 2,
    widths = c(0.25, 1)
)

# panels_legend <- ggarrange(
#     panel_legend,
#     panels,
#     nrow = 2, ncol = 1,
#     heights = c(0.075, 1)
# )


# SAVE THE PLOTS ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/output/"
ggsave(
    filename = paste(saveDir, "Felzartamab all cfDNA panel.png"),
    plot = panels,
    dpi = 600,
    width = 75,
    height = 31,
    units = "cm",
    bg = "white"
)


# END ####
