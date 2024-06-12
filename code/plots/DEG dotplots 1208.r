# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(patchwork) # install.packages("patchwork")
# load DEG results
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/gene_tables_limma_1208.RData")


# WRANGLE DATA FOR PLOTTING ####
df_plot00 <- gene_tables %>%
    mutate(
        design = design %>%
            factor(
                levels = c(
                    "Baseline_vs_Week24",
                    "Week24_vs_Week52",
                    "Baseline_vs_Week52"
                ),
                labels = c(
                    "Baseline - Week24",
                    "Week24 - Week52",
                    "Baseline - Week52"
                )
            )
    ) %>%
    dplyr::select(design, geneset, gene_tables) %>%
    nest(.by = "geneset") %>%
    mutate(
        geneset = geneset %>% factor(levels = c("ABMR_activity", "IFNG", "NK", "Endothelial")),
        data = map(data, unnest, everything())
    ) %>%
    arrange(geneset)


# MAKE DOT PLOTS ####
df_plot <- df_plot00 %>%
    mutate(
        plot_long = pmap(
            list(geneset, data),
            function(geneset, data) {
                data <- data %>%
                    dplyr::rename(logFC = `<U+0394><U+0394> logFC`, p = `<U+0394><U+0394> p`) %>%
                    dplyr::mutate(
                        col = dplyr::case_when(
                            p < 0.05 & logFC < 0 ~ "green",
                            p < 0.05 & logFC > 0 ~ "red",
                            TRUE ~ "grey20"
                        )
                    )
                order_symb <- data %>%
                    dplyr::filter(design == "Baseline - Week24") %>%
                    dplyr::arrange(logFC %>% desc()) %>%
                    dplyr::pull(Symb)
                data %>%
                    dplyr::mutate(Symb = Symb %>% factor(levels = order_symb)) %>%
                    ggplot2::ggplot(mapping = ggplot2::aes(y = Symb, x = logFC)) +
                    ggplot2::geom_vline(xintercept = 0, linetype = "solid", col = "grey10", size = 0.1) +
                    ggplot2::geom_point(col = data$col) +
                    ggplot2::geom_segment(
                        mapping = ggplot2::aes(x = 0, xend = logFC, y = Symb),
                        col = data$col, size = 0.25
                    ) +
                    ggplot2::labs(
                        x = "\u0394\u0394 logFC",
                        y = NULL,
                        title = paste(geneset %>% str_replace("_", " "), "genes")
                    ) +
                    ggplot2::theme_bw() +
                    ggplot2::theme(
                        panel.grid.minor = element_blank(),
                        panel.grid.major = element_line(size = 0.25),
                        axis.text.x = element_text(colour = "black", size = 6),
                        axis.text.y = element_text(colour = "black", size = 7)
                    ) +
                    ggplot2::facet_wrap(~design, scales = "free_x", nrow = 1)
            }
        ),
        plot_wide = pmap(
            list(geneset, data),
            function(geneset, data) {
                data <- data %>%
                    dplyr::rename(logFC = `<U+0394><U+0394> logFC`, p = `<U+0394><U+0394> p`) %>%
                    dplyr::mutate(
                        col = dplyr::case_when(
                            p < 0.05 & logFC < 0 ~ "green",
                            p < 0.05 & logFC > 0 ~ "red",
                            TRUE ~ "grey20"
                        )
                    )
                order_symb <- data %>%
                    dplyr::filter(design == "Baseline - Week24") %>%
                    dplyr::arrange(logFC) %>%
                    dplyr::pull(Symb)
                data %>%
                    dplyr::mutate(Symb = Symb %>% factor(levels = order_symb)) %>%
                    ggplot2::ggplot(mapping = ggplot2::aes(x = Symb, y = logFC)) +
                    ggplot2::geom_hline(yintercept = 0, linetype = "solid", col = "grey10", size = 0.1) +
                    ggplot2::geom_point(col = data$col) +
                    ggplot2::geom_segment(
                        mapping = ggplot2::aes(y = 0, yend = logFC, x = Symb),
                        col = data$col, size = 0.25
                    ) +
                    ggplot2::labs(
                        y = "\u0394\u0394 logFC",
                        x = NULL,
                        title = paste(geneset %>% str_replace("_", " "), "genes")
                    ) +
                    ggplot2::theme_bw() +
                    ggplot2::theme(
                        panel.grid.minor = element_blank(),
                        panel.grid.major = element_line(size = 0.25),
                        axis.text.y = element_text(colour = "black"),
                        axis.text.x = element_text(colour = "black", size = 7, angle = 90, hjust = 1, vjust = 0.5)
                    ) +
                    ggplot2::facet_wrap(~design, scales = "fixed", nrow = 1)
            }
        )
    )
df_plot$plot_long[[3]]


# MAKE PLOTS PANELS ####
plot_panels_wide <- df_plot %>%
    dplyr::filter(geneset != "ABMR_activity") %>%
    pull(plot_wide) %>%
    wrap_plots() +
    plot_layout(axes = "collect") +

    plot_annotation(tag_levels = "A")


plot_panels_long <- df_plot %>%
    dplyr::filter(geneset != "ABMR_activity") %>%
    pull(plot_long) %>%
    wrap_plots() +
    # plot_layout(axes = "collect") +
    plot_annotation(tag_levels = "A")


# SAVE THE PLOTS ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/output/"
ggsave(
    filename = paste(saveDir, "Selective DEG panel wide.png"),
    plot = plot_panels_wide,
    dpi = 600,
    width = 60,
    height = 8,
    units = "cm",
    bg = "white"
)
ggsave(
    filename = paste(saveDir, "Selective DEG panel long.png"),
    plot = plot_panels_long,
    dpi = 600,
    width = 40,
    height = 8,
    units = "cm",
    bg = "white"
)
