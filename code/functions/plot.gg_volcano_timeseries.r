gg_volcano_timeseries <- function(
    data, data_trimmed, xbreak1, xbreak2, alpha_lines,
    point_size = 2.5, point_size_null = 1.25, col_null = "grey20",
    label_n = 10, label_size = 2,
    x_label1_x, x_label2_x) {
    require(tidyverse)
    top_dn <- data %>%
        dplyr::filter(
            design == "Baseline_vs_Week24",
            logFC < 0
        ) %>%
        dplyr::slice_min(`<U+0394><U+0394> p`, n = label_n) %>%
        pull(AffyID)
    top_up <- data %>%
        dplyr::filter(
            design == "Baseline_vs_Week24",
            logFC > 0
        ) %>%
        dplyr::slice_min(`<U+0394><U+0394> p`, n = label_n) %>%
        pull(AffyID)
    data_labels <- data %>% dplyr::filter(AffyID %in% c(top_dn, top_up))
    data %>%
        ggplot2::ggplot(mapping = ggplot2::aes(x = p, y = logFC)) +
        ggplot2::geom_vline(xintercept = xbreak1, linetype = "solid", linewidth = 0.2) +
        ggplot2::geom_vline(xintercept = xbreak2, linetype = "solid", linewidth = 0.2) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
        ggplot2::geom_vline(xintercept = -log10(0.05), linetype = "dashed") +
        ggplot2::geom_vline(xintercept = -log10(0.05) + min_p, linetype = "dashed") +
        ggplot2::geom_segment(
            inherit.aes = FALSE,
            mapping = ggplot2::aes(x = -Inf, xend = xbreak1, y = -Inf),
            linewidth = 0.2
        ) +
        ggplot2::geom_segment(
            inherit.aes = FALSE,
            mapping = ggplot2::aes(x = xbreak2, xend = Inf, y = -Inf),
            linewidth = 0.2
        ) +
        ggplot2::geom_segment(
            inherit.aes = FALSE,
            mapping = ggplot2::aes(x = -Inf, xend = xbreak1, y = Inf),
            linewidth = 0.2
        ) +
        ggplot2::geom_segment(
            inherit.aes = FALSE,
            mapping = ggplot2::aes(x = xbreak2, xend = Inf, y = Inf),
            linewidth = 0.2
        ) +
        ggplot2::geom_segment(
            inherit.aes = FALSE,
            mapping = ggplot2::aes(x = Inf, y = -Inf, yend = Inf),
            linewidth = 0.2
        ) +
        ggplot2::geom_point(
            data = data %>% dplyr::filter(AffyID %nin% c(probes_supressed, probes_increased)),
            col = col_null,
            size = point_size_null,
            alpha = 0.25
        ) +
        ggplot2::geom_line(
            data = data_trimmed,
            mapping = ggplot2::aes(x = p, y = logFC, group = AffyID),
            col = data_trimmed$col,
            alpha = alpha_lines,
            linewidth = 0.5
        ) +
        ggrepel::geom_label_repel(
            data = data_labels,
            mapping = ggplot2::aes(label = Symb),
            size = label_size,
            fontface = "bold",
            min.segment.length = 0,
            segment.color = "grey20",
            segment.size = 0.25,
            col = "white",
            fill = data_labels %>% dplyr::pull(col),
            label.padding = 0.1,
            box.padding = 0.1,
            label.size = 0.1,
            direction = "y",
            # nudge_x = 0.01,
            nudge_y = ifelse(data_labels$logFC > 0, 0.25, -0.25),
            show.legend = FALSE,
            seed = 42
        ) +
        ggplot2::geom_point(
            data = data_trimmed,
            shape = 21,
            fill = data_trimmed$col,
            alpha = data_trimmed$alpha,
            stroke = 0.125,
            size = point_size
        ) +
        ggplot2::geom_text(
            x = -Inf,
            y = Inf,
            vjust = -0.25,
            hjust = 0,
            label = "Baseline - Week24",
            fontface = "bold.italic"
        ) +
        ggplot2::geom_text(
            x = xbreak2,
            y = Inf,
            vjust = -0.25,
            hjust = 0,
            label = "Week24 - Week52",
            fontface = "bold.italic"
        ) +
        ggplot2::geom_text(
            x = x_label1_x,
            y = -Inf,
            vjust = 1.75,
            hjust = 0.5,
            label = "Effect of treatment\n(\u394\u394 p-value)",
            fontface = "plain"
        ) +
        ggplot2::geom_text(
            x = mean(c(min_p, 1 + min_p, 2 + min_p, 3 + min_p)),
            y = -Inf,
            vjust = 1.75,
            hjust = 0.5,
            label = "Refractory response post-treatment\n(\u394\u394 p-value)"
        ) +
        ggplot2::scale_x_continuous(
            breaks = c(0, 1, 2, 3, min_p, 1 + min_p, 2 + min_p, 3 + min_p),
            labels = 10^-c(0, 1, 2, 3, 0, 1, 2, 3)
        ) +
        ggplot2::scale_y_continuous(
            sec.axis = ggplot2::sec_axis(~., name = "Refractory response post-treatment\n(\u394\u394 log2FC)"),
            expand = c(0, 0.05)
        ) +
        ggplot2::coord_cartesian(clip = "off") +
        ggplot2::labs(
            y = "Effect of treatment\n(\u394\u394 log2FC)",
            x = NULL,
            col = NULL
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            panel.border = ggplot2::element_blank(),
            axis.line.x = ggplot2::element_blank(),
            axis.line.y = ggplot2::element_line(linewidth = 0.2),
            panel.grid = ggplot2::element_blank(),
            legend.position = "none",
            axis.title = ggplot2::element_text(size = 12, face = "plain"),
            axis.text = ggplot2::element_text(colour = "black"),
            plot.margin = ggplot2::unit(c(0.5, 0.1, 1.25, 0.1), "cm"),
            plot.background = ggplot2::element_rect(fill = "grey95", colour = " white")
        )
}
