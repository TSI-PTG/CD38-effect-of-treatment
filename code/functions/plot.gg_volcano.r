gg_volcano <- function(
    data, design = NULL, x_break = 1.5, 
    point_size = 2.5, point_size_null = 1.25,
    alpha = 0.8, alpha_null = 0.5, alpha_lines = 0.25,
    col_dn = "#005EFF", col_up = "#ff0040", col_null = "grey30",
    labels_probes = NULL, labels_probes_n = 10, labels_probes_size = 2) {
    require(tidyverse)
    probes_sig_up <- data %>%
        dplyr::filter(
            `<U+0394><U+0394> p` < 0.05,
            `<U+0394><U+0394> logFC` > 0
        ) %>%
        dplyr::pull(AffyID)
    probes_sig_dn <- data %>%
        dplyr::filter(
            `<U+0394><U+0394> p` < 0.05,
            `<U+0394><U+0394> logFC` < 0
        ) %>%
        dplyr::pull(AffyID)
    data <- data %>%
        dplyr::mutate(
            col = dplyr::case_when(
                AffyID %in% probes_sig_dn ~ col_dn,
                AffyID %in% probes_sig_up ~ col_up,
                TRUE ~ col_null
            ),
            order = col %>% factor(levels = c(col_null, col_dn, col_up))
        ) %>%
        dplyr::arrange(order)
    data_sig <- data %>%
        dplyr::filter(AffyID %in% c(probes_sig_up, probes_sig_dn))
    data_labels <- data %>%
        dplyr::filter(AffyID %in% labels_probes)
    min_p <- data %>%
        dplyr::slice_min(`<U+0394><U+0394> p`) %>%
        dplyr::pull(`<U+0394><U+0394> p`) %>%
        log10() * -1
    data %>%
        ggplot2::ggplot(mapping = ggplot2::aes(x = p, y = logFC)) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
        ggplot2::geom_vline(xintercept = -log10(0.05), linetype = "dashed") +
        ggplot2::geom_point(
            data = data %>% dplyr::filter(AffyID %nin% c(probes_sig_up, probes_sig_dn)),
            col = col_null,
            size = point_size_null,
            alpha = 0.25
        ) +
        ggplot2::geom_point(
            data = data_sig,
            shape = 21,
            fill = data_sig$col,
            alpha = alpha,
            stroke = 0.125,
            size = point_size
        ) +
        ggrepel::geom_label_repel(
            data = data_labels,
            mapping = ggplot2::aes(label = Symb),
            size = labels_probes_size,
            fontface = "bold",
            min.segment.length = 0,
            segment.color = "grey20",
            segment.size = 0.25,
            col = "white",
            fill = data_labels %>% dplyr::pull(col),
            label.padding = 0.1,
            box.padding = 0.1,
            point.padding = 0.5,
            label.size = 0.1,
            direction = "y",
            # nudge_x = 0.01,
            nudge_y = ifelse(data_labels$logFC > 0, 0.25, -0.25),
            show.legend = FALSE,
            seed = 42
        ) +
        ggplot2::geom_text(
            x = -Inf,
            y = Inf,
            vjust = -0.25,
            hjust = 0,
            label = design  %>% str_replace("_vs_", " - "),
            fontface = "bold.italic"
        ) +
        ggplot2::geom_text(
            x = x_break,
            y = -Inf,
            vjust = 1.75,
            hjust = 0.5,
            label = "Effect of treatment\n(\u394\u394 p-value)",
            fontface = "plain"
        ) +
        ggplot2::scale_colour_manual(values = c(col_dn, col_up)) +
        ggplot2::scale_x_continuous(
            breaks = c(0, 1, 2, 3, 4, 5),
            labels = 10^-c(0, 1, 2, 3, 4, 5)
        ) +
        ggplot2::coord_cartesian(clip = "off") +
        ggplot2::labs(
            y = "Effect of treatment\n(\u394\u394 log2FC)",
            # x = "Effect of treatment\n(\u394\u394 p-value)",
            x = NULL,
            col = NULL
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            # panel.border = ggplot2::element_blank(),
            axis.line = ggplot2::element_line(linewidth = 0.2),
            panel.grid = ggplot2::element_blank(),
            legend.position = "none",
            axis.title = ggplot2::element_text(size = 12, face = "plain"),
            axis.text = ggplot2::element_text(colour = "black"),
            plot.margin = ggplot2::unit(c(0.5, 0.1, 1.25, 0.1), "cm"),
            plot.background = ggplot2::element_rect(fill = "grey95", colour = " white")
        )
}
