gg_volcano_enrichment <- function(
    data_plot, data_annotation, GO_lines, design = NULL, x_break = 1.5,
    point_size = 2.5, point_size_null = 1.25,
    xlim = NULL, ylim = NULL,
    alpha = 0.8, alpha_null = 0.5, alpha_lines = 0.25,
    col_dn = "#005EFF", col_up = "#ff0040", col_null = "grey30",
    labels_probes = NULL, labels_probes_n = 10, labels_probes_size = 2,
    curvature = NULL, angle = 30, ties = FALSE) {
    require(tidyverse)
    "%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
    ylim <- data_plot$logFC %>% range()
    xlim <- c(0, data_plot$p %>% max() * 1.4)
    probes_sig_up <- data_plot %>%
        dplyr::filter(
            `<U+0394><U+0394> p` < 0.05,
            `<U+0394><U+0394> logFC` > 0
        ) %>%
        dplyr::pull(AffyID)
    probes_sig_dn <- data_plot %>%
        dplyr::filter(
            `<U+0394><U+0394> p` < 0.05,
            `<U+0394><U+0394> logFC` < 0
        ) %>%
        dplyr::pull(AffyID)
    data_plot <- data_plot %>%
        dplyr::mutate(
            col = dplyr::case_when(
                AffyID %in% probes_sig_dn ~ col_dn,
                AffyID %in% probes_sig_up ~ col_up,
                TRUE ~ col_null
            ),
            order = col %>% factor(levels = c(col_null, col_dn, col_up))
        ) %>%
        dplyr::arrange(order)
    data_sig <- data_plot %>%
        dplyr::filter(AffyID %in% c(probes_sig_up, probes_sig_dn))
    data_insig <- data_plot %>%
        dplyr::filter(AffyID %nin% c(probes_sig_up, probes_sig_dn))
    # y_min <- data_plot %>%
    #     dplyr::pull(logFC) %>%
    #     min()
    # y_max <- data_plot %>%
    #     dplyr::pull(logFC) %>%
    #     max()
    # y_diff <- c(y_min, y_max) %>% diff()
    # x_end <- data_plot %>%
    #     dplyr::pull(p) %>%
    #     max() * 1.1
    # n_group <- data_annotation %>%
    #     dplyr::pull(n_group) %>%
    #     unique()
    # interval <- (y_max - (y_min)) / (n_group + 1)
    # GO_lines <- data_annotation %>%
    #             dplyr::arrange(count %>% dplyr::desc()) %>%
    #             dplyr::distinct(Symb, group, .keep_all = TRUE) %>%
    #             dplyr::mutate(
    #                 x_end = x_end,
    #                 y_end = y_min + (groupID * interval),
    #                 curvature = dplyr::case_when(
    #                     groupID == 1 ~ -0.25,
    #                     groupID == 2 ~ 0.25,
    #                     groupID == 3 ~ -0.25,
    #                     groupID == 4 ~ 0.25,
    #                     TRUE ~ 0.25
    #                 ),
    #                 conflict = dplyr::case_when(
    #                     logFC < 0 & NES > 0 ~ FALSE,
    #                     logFC > 0 & NES < 0 ~ FALSE,
    #                     TRUE ~ TRUE
    #                 )
    #             ) %>%
                # dplyr::filter(conflict)
    GO_labels <- data_annotation %>%
        dplyr::distinct(Symb, group, .keep_all = TRUE) %>%
        dplyr::slice_max(prop, n = 5, by = c("group"), with_ties = FALSE) %>%
        dplyr::arrange(group, prop %>% dplyr::desc()) %>%
        dplyr::left_join(GO_lines %>% dplyr::distinct(Symb, group, .keep_all = TRUE))
    data_plot %>%
        ggplot2::ggplot(mapping = ggplot2::aes(x = p, y = logFC)) +
        # purrr::map(
        #     GO_lines$curvature,
        #     function(i) {
        #         ggplot2::geom_curve(
        #             data = GO_lines,
        #             mapping = ggplot2::aes(
        #                 x = p,
        #                 y = logFC,
        #                 xend = x_end,
        #                 yend = y_end,
        #                 col = group,
        #                 group = group,
        #             ),
        #             curvature = i,
        #             angle = 30,
        #             linewidth = 0.1
        #         )
        #     }
        # ) +
        ggplot2::geom_curve(
            data = GO_lines,
            mapping = ggplot2::aes(
                x = p,
                y = logFC,
                xend = x_end,
                yend = y_end,
                # col = group,
                group = group,
            ),
            col = GO_lines$col_group,
            # curvature = curvature,
            angle = 30,
            linewidth = 0.1
        ) +
        ggnewscale::new_scale_colour() +
        ggplot2::geom_segment(
            inherit.aes = FALSE,
            data = dplyr::tibble(
                x0 = seq(
                    -log10(0.05),
                    GO_lines$x_end %>% max(),
                    length.out = 300
                ),
                xend = x0,
                y0 = -Inf,
                yend = Inf,
                col = x0
            ),
            mapping = ggplot2::aes(
                x = x0,
                xend = xend,
                y = y0,
                yend = yend,
                col = x0
            ),
            show.legend = FALSE
        ) +
        ggplot2::scale_colour_gradient(low = "#ffffff95", high = "#ffffff00") +
        ggnewscale::new_scale_colour() +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
        ggplot2::geom_vline(xintercept = -log10(0.05), linetype = "dashed") +
        ggplot2::geom_point(
            data = data_insig,
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
        #     # ggrepel::geom_label_repel(
        #     #     data_plot = GO_labels,
        #     #     mapping = ggplot2::aes(label = Symb),
        #     #     size = labels_probes_size,
        #     #     fontface = "bold",
        #     #     min.segment.length = 0,
        #     #     segment.color = "grey20",
        #     #     segment.size = 0.25,
        #     #     col = "white",
        #     #     fill = GO_labels %>% dplyr::pull(col),
        #     #     label.padding = 0.1,
        #     #     box.padding = 0.1,
        #     #     point.padding = 0.5,
        #     #     label.size = 0.1,
        #     #     direction = "y",
        #     #     # nudge_x = 0.01,
        #     #     nudge_y = ifelse(GO_labels$logFC > 0, 0.25, -0.25),
        #     #     show.legend = FALSE,
        #     #     seed = 42
        #     # ) +
        ggplot2::geom_text(
            x = -Inf,
            y = Inf,
            vjust = -0.25,
            hjust = 0,
            label = design %>% stringr::str_replace("_vs_", " - "),
            fontface = "bold.italic"
        ) +
        ggplot2::scale_colour_manual(values = c(col_dn, col_up)) +
        ggplot2::scale_x_continuous(
            breaks = c(0, 1, 2, 3, 4, 5),
            labels = 10^-c(0, 1, 2, 3, 4, 5)
        ) +
        coord_cartesian(
            xlim = xlim,
            ylim = ylim,
            clip = "off"
        ) +
        ggplot2::labs(
            y = "\u394\u394 logFC",
            x = "\u394\u394 p-value",
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
            plot.margin = ggplot2::unit(c(0.5, 0.1, 0.1, 0.1), "cm"),
            # plot.background = ggplot2::element_rect(fill = "grey95", colour = " white")
        ) %>%
        suppressWarnings()
}
