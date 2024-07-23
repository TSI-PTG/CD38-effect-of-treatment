# HOUSEKEEPING ####
# CRAN packages
library(tidyverse) # install.packages("tidyverse")
library(flextable) # install.packages("flextable")
library(ggpubr) # install.packages("ggpubr")
library(patchwork) # install.packages("patchwork")
library(ggprism) # install.packages("ggprism")
# Custom operators, functions, and datasets
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
# load injury results
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/results_injury_k1208.RData")
# load violin plots
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/Felzartamab_plots.RData")
# load volcano enrichment plots
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/Volcano enrichment plots baseline to week 52 IQR_filtered_probes_unique_genes_baseline_corrected_cortex_corrected 1208.RData")



# DEFINE THE REFERENCE INJURY DATA ###
m <- resuls_injury$reference_data[[1]]


# MAKE PLOTS OF IRRAT30 MODEL RESULTS ####
# TODO: make iterative workflow for these plots
plot_IRRAT30 <- resuls_injury %>%
    dplyr::filter(variable == "IRRAT30") %>%
    pull(plot_data) %>%
    pluck(1) %>%
    ggplot(aes(x = x, y = predicted, color = group)) +
    geom_line(aes(group = group), alpha = 1, linewidth = 1.2) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.4, color = NA) +
    geom_point(
        data = m,
        aes(x = time, y = IRRAT30, color = factor(Felzartamab), group = Patient),
        alpha = 0,
        size = 0.1
    ) +
    geom_line(
        data = m,
        aes(x = time, y = IRRAT30, color = factor(Felzartamab), group = Patient, linetype = linetype),
        stat = "smooth",
        method = "lm",
        alpha = 0.3,
        show.legend = FALSE
    ) +
    # annotate(
    #   geom = "text",
    #   label = "Slope difference: -0.60 (95%CI -1.10 to -0.10), p=0.055",
    #   x = 0.5, y = 1.8, size = 4
    # ) +
    scale_color_manual(
        name = "Treatment",
        values = c("Placebo" = "black", "Felzartamab" = "#bc3c29"),
        labels = c("Placebo", "Felzartamab")
    ) +
    scale_fill_manual(
        name = "Treatment",
        values = c("Placebo" = "#b3b3b3", "Felzartamab" = "#bc3c29"),
        labels = c("Placebo", "Felzartamab")
    ) +
    scale_x_continuous(
        breaks = c(0, 0.46, 0.99),
        labels = c("0", "24", "52")
    ) +
    scale_y_continuous(limits = c(-1.5, 2)) +
    labs(
        # title = "D",
        x = "Time post-treatment (weeks)",
        y = "Injury-repair associated\n(IRRAT30)"
    ) +
    theme_classic() +
    theme(
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        legend.key.width = unit(4, "cm"),
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 16, face = "bold")
    )


# MAKE PLOTS OF IRITD3 MODEL RESULTS ####
# TODO: make iterative workflow for these plots

plot_IRITD3 <- resuls_injury %>%
    dplyr::filter(variable == "IRITD3") %>%
    pull(plot_data) %>%
    pluck(1) %>%
    ggplot(aes(x = x, y = predicted, color = group)) +
    geom_line(aes(group = group), alpha = 1, linewidth = 1.2) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.4, color = NA) +
    geom_point(
        data = m,
        aes(x = time, y = IRITD3, color = factor(Felzartamab), group = Patient),
        alpha = 0,
        size = 0.1
    ) +
    geom_line(
        data = m,
        aes(x = time, y = IRITD3, color = factor(Felzartamab), group = Patient, linetype = linetype),
        stat = "smooth",
        method = "lm",
        alpha = 0.3,
        show.legend = FALSE
    ) +
    # annotate(
    #   geom = "text",
    #   label = "Slope difference: -0.20 (95%CI -0.37 to -0.02), p=0.041",
    #   x = 0.5, y = 0.65, size = 4
    # ) +
    labs(
        # title = "E",
        x = "Time post-treatment (weeks)",
        y = "Injury-repair induced, day 3\n(IRITD3)"
    ) +
    scale_color_manual(
        name = "Treatment",
        values = c("Placebo" = "black", "Felzartamab" = "#bc3c29"),
        labels = c("Placebo", "Felzartamab")
    ) +
    scale_fill_manual(
        name = "Treatment",
        values = c("Placebo" = "#b3b3b3", "Felzartamab" = "#bc3c29"),
        labels = c("Placebo", "Felzartamab")
    ) +
    scale_x_continuous(
        breaks = c(0, 0.46, 1),
        labels = c("0", "24", "52")
    ) +
    scale_y_continuous(limits = c(-0.7, 0.7)) +
    theme_classic() +
    theme(
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 16, face = "bold")
    )


# MAKE PLOTS OF IRITD5 MODEL RESULTS ####
# TODO: make iterative workflow for these plots

plot_IRITD5 <- resuls_injury %>%
    dplyr::filter(variable == "IRITD5") %>%
    pull(plot_data) %>%
    pluck(1) %>%
    ggplot(aes(x = x, y = predicted, color = group)) +
    geom_line(aes(group = group), alpha = 1, linewidth = 1.2) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.4, color = NA) +
    geom_point(
        data = m,
        aes(x = time, y = IRITD5, color = factor(Felzartamab), group = Patient),
        alpha = 0,
        size = 0.1
    ) +
    geom_line(
        data = m,
        aes(x = time, y = IRITD5, color = factor(Felzartamab), group = Patient, linetype = linetype),
        stat = "smooth",
        method = "lm",
        alpha = 0.3,
        show.legend = FALSE
    ) +
    # annotate(
    #   geom = "text",
    #   label = "Slope difference: -0.23 (95%CI -0.41 to -0.05), p=0.042",
    #   x = 0.5, y = 0.75, size = 4
    # ) +
    labs(
        # title = "F",
        x = "Time post-treatment (weeks)",
        y = "Injury-repair induced, day 5\n(IRITD5)"
    ) +
    scale_color_manual(
        name = "Treatment",
        values = c("Placebo" = "black", "Felzartamab" = "#bc3c29"),
        labels = c("Placebo", "Felzartamab")
    ) +
    scale_fill_manual(
        name = "Treatment",
        values = c("Placebo" = "#b3b3b3", "Felzartamab" = "#bc3c29"),
        labels = c("Placebo", "Felzartamab")
    ) +
    scale_x_continuous(
        breaks = c(0, 0.46, 0.99),
        labels = c("0", "24", "52")
    ) +
    scale_y_continuous(limits = c(-0.20, 0.8)) +
    theme_classic() +
    theme(
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 16, face = "bold")
    )


# GENERATE GROBS FOR SLOPE TABLES ####
table_slope_IRRAT30_grob <- resuls_injury %>%
    dplyr::filter(variable == "IRRAT30") %>%
    pull(slope_tables) %>%
    pluck(1) %>%
    gen_grob(fit = "fixed", just = "centre") %>%
    as_ggplot()

table_slope_IRITD3_grob <- resuls_injury %>%
    dplyr::filter(variable == "IRITD3") %>%
    pull(slope_tables) %>%
    pluck(1) %>%
    gen_grob(fit = "fixed", just = "centre") %>%
    as_ggplot()

table_slope_IRITD5_grob <- resuls_injury %>%
    dplyr::filter(variable == "IRITD5") %>%
    pull(slope_tables) %>%
    pluck(1) %>%
    gen_grob(fit = "fixed", just = "centre") %>%
    as_ggplot()


# EXTRACT LEGEND FOR VIOLIN PLOTS ####
panel_legend_violin <- felzartamab_plots %>%
    dplyr::filter(variable == "AMAT1") %>%
    pull(plot_violin) %>%
    ggpubr::get_legend() %>%
    ggpubr::as_ggplot() +
    theme(plot.margin = unit(c(0, 0, -1, 0), "cm"))


# EXTRACT LEGEND FOR REGRESSION PLOTS ####
panel_legend_regression <- plot_IRRAT30 %>%
    ggpubr::get_legend() %>%
    ggpubr::as_ggplot()


# EXTRACT JOINT LEGEND FOR REGRESSION AND VIOLIN PLOTS ####
panel_legend_all_injury <- panel_legend_violin + panel_legend_regression



# EXTRACT INJURY VIOLIN PLOTS ####
plots_violin <- felzartamab_plots %>%
    dplyr::filter(category %in% c("injury")) %>%
    pull(plot_violin)


# MAKE PANEL OF SLOPE PLOTS ####
plot_panels_regression <- patchwork::wrap_plots(
    plot_IRRAT30, plot_IRITD3, plot_IRITD5,
    table_slope_IRRAT30_grob, table_slope_IRITD3_grob, table_slope_IRITD5_grob,
    nrow = 2,
    ncol = 3,
    guides = "collect",
    heights = c(1, 0.5)
) +
    patchwork::plot_annotation(tag_levels = list(c("D", "E", "F", rep("", 3)))) &
    theme(
        legend.position = "none",
        plot.title = element_blank(),
        axis.text = element_text(size = 10, colour = "black"),
        plot.tag = element_text(size = 20, face = "bold", vjust = 1)
    )

plot_panels_regression_legend <- ggarrange(
    panel_legend_regression,
    plot_panels_regression,
    nrow = 2,
    heights = c(0.15, 1.5)
)


# MAKE PANEL OF VIOLIN PLOTS ####
# plot_panel_violin <- plots_violin %>%
#   wrap_plots(nrow = 1, ncol = 3) +
#   plot_annotation(tag_levels = list(c(LETTERS[1:15]))) &
#   theme(
#     legend.position = "none",
#     axis.text = element_text(size = 10, colour = "black"), plot.tag = element_text(size = 20, face = "bold", vjust = 1)
#   )

# panels_violin_legend <- ggarrange(
#   panel_legend_violin,
#   plot_panel_violin,
#   nrow = 2,
#   heights = c(0.25, 1.5)
# )


# DEFINE VOLCANO ENRICHMENT PLOT ####
plot_volcano <- plot_volcano_enrichment %>%
    pull() %>%
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
# plot_panels_all_injury <- patchwork::wrap_plots(
#     A = plots_violin[[1]],
#     B = plots_violin[[2]],
#     C = plots_violin[[3]],
#     D = plot_IRRAT30,
#     E = plot_IRITD3,
#     F = plot_IRITD5,
#     I = table_slope_IRRAT30_grob,
#     J = table_slope_IRITD3_grob,
#     K = table_slope_IRITD5_grob,
#     # nrow = 3,
#     # ncol = 3,
#     design = "ABC
#     DEF
#     IJK
#     GGG",
#     guides = "collect",
#     heights = c(1, 1, 0.5, 2.5)
# ) +
#     patchwork::plot_annotation(tag_levels = list(c(LETTERS[1:6], rep("", ), "G"))) &
#     theme(
#         legend.position = "none",
#         # plot.title = element_blank(),
#         axis.text = element_text(size = 10, colour = "black"),
#         axis.title = element_text(size = 12, colour = "black"),
#         plot.tag = element_text(size = 20, face = "bold", vjust = 1)
#     )

plot_panels_all_injury <- patchwork::wrap_plots(
    plots_violin[[1]], plots_violin[[2]], plots_violin[[3]],
    plot_IRRAT30, plot_IRITD3, plot_IRITD5,
    table_slope_IRRAT30_grob, table_slope_IRITD3_grob, table_slope_IRITD5_grob,
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


panels_all_injury <- ggarrange(
    panel_legend_all_injury,
    plot_panels_all_injury,
    plot_volcano,
    nrow = 3,
    heights = c(0.15, 1.5, 1.5)
)


# SAVE THE PLOTS ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/output/"
ggsave(
    filename = paste(saveDir, "Felzartamab injury ART, regression, volcano.png"),
    plot = panels_all_injury,
    dpi = 300,
    width = 40,
    height = 50,
    units = "cm",
    bg = "white"
)
