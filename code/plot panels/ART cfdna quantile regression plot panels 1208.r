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
panel_legend <- felzartamab_plots %>%
    dplyr::filter(variable == "AMAT1") %>%
    pull(plot_violin) %>%
    get_legend() %>%
    as_ggplot()


# MAKE PANEL OF VIOLIN PLOTS ####
plot_cfdna <- felzartamab_plots %>%
    dplyr::filter(category %in% c("cfDNA")) %>%
    pull(plot_violin) %>%
    pluck(1) +
    facet_wrap(~Felzartamab, nrow = 2, ncol = 1, scales = "free_x")

panel_violin <- plot_cfdna %>%
    wrap_plots(., nrow = 1, ncol = 1) +
    plot_annotation(
        tag_levels = list(c(LETTERS[1:15]))
    ) &
    theme(
        legend.position = "none",
        axis.text = element_text(size = 12, colour = "black"),
        plot.tag = element_text(size = 20, face = "bold", vjust = 1)
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
        title = "Baseline - Week24",
        tag_levels = list(c("B", rep("", 9)))
    ) &
    theme(
        legend.position = "none",
        plot.title = element_text(size = 20, face = "bold.italic"),
        axis.text = element_text(size = 10, colour = "black"),
        plot.tag = element_text(size = 20, face = "bold", vjust = 1)
    )

panel_cfdna_correlation_2 <- felzartamab_cfdna_qr_plots %>%
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
        plot.title = element_text(size = 20, face = "bold.italic"),
        axis.text = element_text(size = 10, colour = "black"),
        plot.tag = element_text(size = 20, face = "bold", vjust = 1),
        plot.background = element_rect(fill = "grey95", colour = " white")
    )


panel_cfdna_correlation_3 <- felzartamab_cfdna_qr_plots %>%
    dplyr::filter(
        Followup_pairwise %in% c("Baseline - Week52"),
        category %in% c("ABMR")
    ) %>%
    pull(plot_scatter) %>%
    wrap_plots(nrow = 1, ncol = 5) +
    plot_annotation(
        title = "Baseline - Week52",
        tag_levels = list(c("D", rep("", 9)))
    ) &
    theme(
        legend.position = "none",
        plot.title = element_text(size = 20, face = "bold.italic"),
        axis.text = element_text(size = 10, colour = "black"),
        plot.tag = element_text(size = 20, face = "bold", vjust = 1)
    )


panels_quantile_regression <- ggarrange(
    panel_cfdna_correlation_1,
    panel_cfdna_correlation_2,
    panel_cfdna_correlation_3,
    nrow = 3
) %>%
    ggpubr::annotate_figure(
        top = text_grob(
            "Relationship of dd-cfDNA to ABMR Activity Scores",
            face = "bold.italic", size = 25, hjust = 1.4
        )
    )


panel <- ggarrange(
    panel_violin,
    panels_quantile_regression,
    nrow = 1, ncol = 2,
    widths = c(0.25, 1)
)



# SAVE THE PLOTS ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/output/"
ggsave(
    filename = paste(saveDir, "Felzartamab all cfDNA panel.png"),
    plot = panel,
    dpi = 300,
    width = 75,
    height = 30,
    units = "cm",
    bg = "white"
)


# END ####
