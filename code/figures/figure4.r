# HOUSEKEEPING ####
# CRAN packages
library(tidyverse) # install.packages("tidyverse")
library(flextable) # install.packages("flextable")
library(ggpubr) # install.packages("ggpubr")
library(patchwork) # install.packages("patchwork")
library(ggprism) # install.packages("ggprism")
# Custom operators and functions
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
# load regression plots
load("results/plots_regression.RData")
# load violin plots
load("results/plots_artANOVA.RData")
# load volcano enrichment plots
load("results/plots_volcano.RData")


# DEFINE LEGEND FOR VIOLIN PLOTS ####
panel_legend_violin <- plots_artANOVA %>%
    dplyr::filter(variable == "AMAT1") %>%
    pull(plot_violin) %>%
    ggpubr::get_legend() %>%
    ggpubr::as_ggplot() +
    theme(plot.margin = unit(c(0, 0, -1, 0), "cm"))


# DEFINE LEGEND FOR REGRESSION PLOTS ####
panel_legend_regression <- plots_regression %>%
    dplyr::filter(variable == "IRRAT30") %>%
    pull(plot) %>%
    ggpubr::get_legend() %>%
    ggpubr::as_ggplot()


# DEFIN JOINT LEGEND FOR REGRESSION AND VIOLIN PLOTS ####
panel_legend_all_injury <- panel_legend_violin + panel_legend_regression


# DEFINE INJURY VIOLIN PLOTS ####
plots_violin <- plots_artANOVA %>%
    dplyr::filter(category %in% c("injury")) %>%
    pull(plot_violin)


# DEFINE VOLCANO ENRICHMENT PLOT ####
plot_volcano <- plot_volcano %>%
    dplyr::filter(design == "Baseline_vs_Week52") %>%
    pull(plot_volcano_enrichment) %>%
    pluck(1) %>% 
    wrap_plots(
        A = plot_spacer(),
        B = .,
        C = plot_spacer(),
        design = c(
            "ABBBBC"
        )
    ) +
    patchwork::plot_annotation(tag_levels = list(c("G"))) &
    theme(
        legend.position = "none",
        plot.title = element_blank(),
        # plot.title = element_text(size = 15),
        axis.text = element_text(size = 10, colour = "black"),
        axis.title = element_text(size = 12, colour = "black"),
        plot.tag = element_text(size = 20, face = "bold", vjust = 1)
    )



# MAKE PANEL OF SLOPE AND VIOLIN PLOTS ####
plot_panels_all_injury <- patchwork::wrap_plots(
    plots_violin[[1]], plots_violin[[2]], plots_violin[[3]],
    plots_regression$plot[[1]], plots_regression$plot[[2]], plots_regression$plot[[3]],
    plots_regression$grob_table[[1]], plots_regression$grob_table[[2]], plots_regression$grob_table[[3]],
    nrow = 3,
    ncol = 3,
    guides = "collect",
    heights = c(1, 1, 0.5)
) +
    patchwork::plot_annotation(tag_levels = list(c(LETTERS[1:6], rep("", 3)))) &
    theme(
        legend.position = "none",
        # plot.title = element_blank(),
        axis.text = element_text(size = 10, colour = "black"),
        axis.title = element_text(size = 12, colour = "black"),
        plot.tag = element_text(size = 20, face = "bold", vjust = 1)
    )


figure4 <- ggarrange(
    panel_legend_all_injury,
    plot_panels_all_injury,
    plot_volcano,
    nrow = 3,
    heights = c(0.15, 1.5, 1.5)
)


# SAVE THE PLOTS ####
saveDir <- "png/"
ggsave(
    filename = paste0(saveDir, "figure4.png"),
    plot = figure4,
    dpi = 300,
    width = 40,
    height = 50,
    units = "cm",
    bg = "white"
)
