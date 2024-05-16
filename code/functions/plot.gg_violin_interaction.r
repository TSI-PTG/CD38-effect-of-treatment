gg_violin_interaction <- function(data, variable, score, medians_delta, art_con_interaction_default_tidy) {
    require(tidyverse)
    require(gghalves)
    require(ggprism)
    log10zero <- scales::trans_new(
        name = "log10zero",
        transform = function(x) log10(x + 0.15),
        inverse = function(x) 10^x - 0.15
    )
    if (score %>% stringr::str_detect("Prob")) {
        ylims <- c(0, 1)
        point_size <- 2
    } else if (variable == "cfDNA") {
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
    dodge <- 0.3
    delta <- medians_delta %>%
        dplyr::mutate(Felzartamab = Felzartamab %>% factor(labels = c("Placebo", "Felzartamab")))
    delta_delta <- medians_delta %>%
        dplyr::distinct(Followup_pairwise, .keep_all = TRUE) %>%
        dplyr::mutate(Followup_pairwise = c("Week24 - Day0", "Week52 - Day0", "Week52 - Week24"))
    delta_delta_p <- art_con_interaction_default_tidy %>%
        dplyr::select(Followup_pairwise, adj.p.value) %>%
        dplyr::mutate(
            FDR = ifelse(
                adj.p.value < 0.01,
                formatC(adj.p.value, digits = 0, format = "e"),
                formatC(adj.p.value, digits = 3, format = "f")
            )
        )
    data <- data %>%
        dplyr::mutate(
            variable = variable,
            Felzartamab = Felzartamab %>% factor(labels = c("Placebo", "Felzartamab"))
        )
    data <- data %>%
        dplyr::left_join(
            data %>%
                dplyr::reframe(median = median(value), .by = c(Followup, Felzartamab, Patient)) %>%
                rstatix::spread(Followup, median) %>%
                dplyr::mutate(
                    delta = dplyr::case_when(
                        variable == "cfDNA" ~ log2(Week24 / Day0),
                        TRUE ~ Week24 - Day0
                    ),
                    delta = dplyr::case_when(
                        delta == -Inf ~ min(delta[delta != -Inf]),
                        delta == Inf ~ max(delta[delta != Inf]),
                        TRUE ~ delta
                    ),
                    delta2 = dplyr::case_when(
                        variable == "cfDNA" ~ log2(Week52 / Week24),
                        TRUE ~ Week52 - Week24
                    ),
                    delta2 = dplyr::case_when(
                        delta2 == -Inf ~ min(delta2[delta2 != -Inf]),
                        delta2 == Inf ~ max(delta2[delta2 != Inf]),
                        TRUE ~ delta2
                    )
                ) %>%
                tidyr::pivot_longer(cols = c(Day0, Week24, Week52), names_to = "Followup", values_to = "value") %>%
                dplyr::select(Felzartamab, Patient, Followup, delta, delta2),
            by = c("Felzartamab", "Patient", "Followup")
        ) %>%
        dplyr::mutate(Followup = Followup %>% factor(levels = c("Day0", "Week24", "Week52")))
    bw <- ifelse(data$variable[[1]] %in% c("ABMRpm", "ggt0", "ptcgt0", "TCMRt", "tgt1", "igt1"), 0.05, 0.1)
    midpoint <- 0
    min_delta <- data %>%
        dplyr::select(delta, delta2) %>%
        flatten_dbl() %>%
        min()
    max_delta <- data %>%
        dplyr::select(delta, delta2) %>%
        flatten_dbl() %>%
        max()
    plot <- data %>%
        ggplot2::ggplot(aes(x = Followup, y = value)) +
        gghalves::geom_half_violin(
            inherit.aes = FALSE,
            data = data %>% dplyr::filter(Followup %in% c("Day0")),
            mapping = ggplot2::aes(
                x = Followup,
                y = value
            ),
            bw = bw,
            side = c("l"),
            fill = "grey95",
            col = "grey60",
            trim = FALSE,
            scale = "width"
        ) +
        ggplot2::geom_violin(
            inherit.aes = FALSE,
            data = data %>% dplyr::filter(Followup == c("Week24")),
            mapping = ggplot2::aes(
                x = Followup,
                y = value
            ),
            bw = bw,
            fill = "#f2f2f2ec",
            col = "#cdcdcdb8",
            trim = FALSE,
            scale = "width"
        ) +
        gghalves::geom_half_violin(
            inherit.aes = FALSE,
            data = data %>% dplyr::filter(Followup %in% c("Week52")),
            mapping = ggplot2::aes(
                x = Followup,
                y = value
            ),
            bw = bw,
            side = c("r"),
            fill = "grey95",
            col = "grey60",
            trim = FALSE,
            scale = "width"
        ) +
        ggplot2::geom_point(
            inherit.aes = FALSE,
            data = data %>% dplyr::filter(Followup %in% c("Day0")),
            mapping = ggplot2::aes(
                x = Followup,
                y = value,
                col = delta
            ),
            position = ggplot2::position_nudge(x = 0.1),
            size = 2,
            alpha = 0.75
        ) +
        ggplot2::geom_point(
            inherit.aes = FALSE,
            data = data %>% dplyr::filter(Followup %in% c("Week24")),
            mapping = ggplot2::aes(
                x = Followup,
                y = value,
                col = delta
            ),
            size = 2,
            alpha = 0.75
        ) +
        ggplot2::geom_point(
            inherit.aes = FALSE,
            data = data %>% dplyr::filter(Followup %in% c("Week52")),
            mapping = ggplot2::aes(
                x = Followup,
                y = value,
                col = delta2
            ),
            position = ggplot2::position_nudge(x = -0.1),
            size = 2,
            alpha = 0.75
        ) +
        ggplot2::geom_line(
            data = data %>% dplyr::filter(Followup %in% c("Day0", "Week24")),
            mapping = aes(
                x = Followup,
                col = delta,
                group = Patient
            ),
            position = ggplot2::position_nudge(
                x = ifelse(data %>%
                    dplyr::filter(Followup %in% c("Day0", "Week24")) %>%
                    dplyr::pull(Followup) == "Day0", 0.1, 0)
            ),
            linewidth = 0.5,
            alpha = 0.25,
            show.legend = FALSE
        ) +
        ggplot2::geom_line(
            data = data %>% dplyr::filter(Followup %in% c("Week24", "Week52")),
            mapping = ggplot2::aes(
                x = Followup,
                col = delta2,
                group = Patient
            ),
            position = ggplot2::position_nudge(
                x = ifelse(data %>%
                    dplyr::filter(Followup %in% c("Week24", "Week52")) %>%
                    dplyr::pull(Followup) == "Week24", 0, -0.1)
            ),
            linewidth = 0.5,
            alpha = 0.25,
            show.legend = FALSE
        ) +
        ggplot2::stat_summary(fun = median, geom = "point", size = 4) +
        ggplot2::stat_summary(
            ggplot2::aes(group = Felzartamab),
            fun = median,
            geom = "line",
            linewidth = 1,
            linetype = "dashed"
        ) +
        ggplot2::labs(
            title = paste(
                paste(
                    "\u394\u394", delta_delta$Followup_pairwise[[1]], "=",
                    delta_delta$median_delta_delta[[1]] %>% round(2),
                    "| FDR =",
                    delta_delta_p$FDR[[1]]
                ),
                paste(
                    "\u394\u394", delta_delta$Followup_pairwise[[3]], "=",
                    delta_delta$median_delta_delta[[3]] %>% round(2),
                    "| FDR =",
                    delta_delta_p$FDR[[3]]
                ),
                paste(
                    "\u394\u394", delta_delta$Followup_pairwise[[2]], "=",
                    delta_delta$median_delta_delta[[2]] %>% round(2),
                    "| FDR =",
                    delta_delta_p$FDR[[2]]
                ),
                sep = "\n"
            ),
            x = NULL,
            y = score %>% stringr::str_replace("\\(", "\n("),
            col = "Individual patient response     ",
            parse = TRUE
        ) +
        ggplot2::scale_color_gradient2(
            low = "#00ff00bc",
            mid = "grey60",
            high = "red",
            midpoint = midpoint,
            breaks = c(min(data$delta), max(data$delta)),
            labels = c("improved", "worsened"),
            guide = ggplot2::guide_colorbar(
                title.position = "top",
                barwidth = 20,
                ticks = FALSE,
                label.hjust = c(5.35, min_delta), # first value is improved, second value worsened (i.e., reverse = TRUE)
                label.vjust = 7.75,
                reverse = TRUE
            )
        ) +
        ggplot2::scale_x_discrete(labels = xlabels) +
        ggplot2::coord_cartesian(
            xlim = c(NA, NA),
            ylim = ylims
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            axis.title = ggplot2::element_text(size = 15),
            axis.text = ggplot2::element_text(size = 12, colour = "black"),
            panel.grid = ggplot2::element_blank(),
            strip.text = ggplot2::element_text(size = 12, colour = "black"),
            legend.position = "top",
            legend.title = ggplot2::element_text(size = 15, hjust = 0.66, vjust = 1, face = "bold"),
            legend.text = ggplot2::element_text(size = 15),
            plot.title = ggplot2::element_text(
                colour = "black",
                hjust = 0,
                size = 12,
                face = dplyr::case_when(
                    delta_delta_p$adj.p.value[[1]] < 0.05 ~ "bold.italic",
                    delta_delta_p$adj.p.value[[2]] < 0.05 ~ "bold.italic",
                    TRUE ~ "italic"
                )
            )
        ) +
        ggplot2::facet_wrap(~Felzartamab)
    if (logy) {
        plot <- plot +
            ggplot2::scale_y_continuous(
                trans = log10zero,
                breaks = log10_labels,
                labels = log10_labels,
                minor_breaks = log10_ticks,
                expand = c(0.05, 0.05)
            ) +
            ggplot2::theme(axis.ticks.length.y = unit(0.25, "cm")) +
            ggplot2::guides(y = ggprism::guide_prism_minor())
    } else {
        plot <- plot + ggplot2::scale_y_continuous(expand = c(0, 0), trans = "identity")
    }
    plot
}