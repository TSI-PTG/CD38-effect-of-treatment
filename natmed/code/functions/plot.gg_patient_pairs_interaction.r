gg_patient_pairs_interaction <- function(data, variable, score) {
    require(tidyverse)
    require(gghalves)
    require(ggprism)
    require(ggsci)
    log10zero <- scales::trans_new(
        name = "log10zero",
        transform = function(x) log10(x + 0.15),
        inverse = function(x) 10^x - 0.15
    )
    if (score %>% stringr::str_detect("Prob")) {
        ylims <- c(0, 1)
        point_size <- 2
    } else if (score %>% stringr::str_detect("cfDNA")) {
        ylims <- c(0, 1000)
        point_size <- 2
    } else {
        ylims <- c(NA, NA)
        point_size <- 2
    }
    if (variable == "cfDNA") logy <- TRUE else logy <- FALSE
    log10_ticks <- c(
        seq(0, 1, length.out = 11),
        seq(1, 10, length.out = 10),
        seq(10, 100, length.out = 10),
        seq(100, 1000, length.out = 10)
    )
    log10_labels <- c(0, 1, 10, 100, 1000)
    xlabels <- data %>%
        dplyr::distinct(Followup, .keep_all = TRUE) %>%
        dplyr::select(Followup, sample_pairs) %>%
        dplyr::mutate(xlabels = paste(Followup, "\n(n = ", sample_pairs, ")", sep = "")) %>%
        dplyr::pull(xlabels)
    dodge <- 0.5
    data <- data %>%
        dplyr::mutate(
            variable = variable,
            Felzartamab = Felzartamab %>% factor(labels = c("Placebo", "Felzartamab"))
        )
    n <- data %>%
        dplyr::distinct(Patient) %>%
        nrow()
    cols <- c(ggsci::pal_npg()(10), ggsci::pal_npg(alpha = 0.7)(10), ggsci::pal_npg(alpha = 0.4)(10))[1:n]

    plot <- data %>%
        ggplot2::ggplot(
            mapping = ggplot2::aes(
                x = Followup,
                y = value,
                col = Patient,
                group = Patient
            )
        ) +
        ggplot2::geom_point(
            size = 5,
            alpha = 0.5,
        ) +
        ggplot2::geom_line(
            linewidth = 0.5,
            linetype = "solid",
            show.legend = FALSE
        ) +
        ggplot2::geom_point(
            size = 2,
            alpha = 1,
            col = "black"
        ) +
        ggrepel::geom_label_repel(
            seed = 42,
            mapping = ggplot2::aes(label = Patient),
            col = "black",
            fill = alpha("white", 0.75),
            box.padding = 0.35,
            label.padding = 0.15,
            label.r = 0.5,
            label.size = 0,
            show.legend = FALSE,
            force_pull = 5,
            size = 8
        ) +
        ggplot2::labs(
            x = NULL,
            y = score %>% stringr::str_replace("\\(", "\n("),
            col = "patient",
            parse = TRUE
        ) +
        ggplot2::scale_color_manual(values = cols) +
        ggplot2::scale_x_discrete(labels = xlabels) +
        ggplot2::coord_cartesian(xlim = c(1.3, 2.7)) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            axis.title = ggplot2::element_text(size = 15),
            axis.text = ggplot2::element_text(size = 15, colour = "black"),
            axis.text.x = ggplot2::element_text(size = 15, colour = "black"),
            panel.grid = ggplot2::element_blank(),
            strip.text = ggplot2::element_text(size = 15, colour = "black"),
            legend.position = "top",
            legend.title = ggplot2::element_text(size = 15, face = "bold"),
            legend.text = ggplot2::element_text(size = 15),
        ) +
        ggplot2::facet_wrap(~Felzartamab) +
        ggplot2::guides(col = ggplot2::guide_legend(nrow = 1))
    if (logy) {
        plot <- plot +
            ggplot2::scale_y_continuous(
                trans = log10zero,
                breaks = log10_labels,
                labels = log10_labels,
                minor_breaks = log10_ticks,
                limits = c(0, 1000),
                expand = c(0.05, 0.05)
            ) +
            ggplot2::theme(axis.ticks.length.y = unit(0.25, "cm")) +
            ggplot2::guides(y = ggprism::guide_prism_minor())
    } else {
        plot
    }
    plot
}
