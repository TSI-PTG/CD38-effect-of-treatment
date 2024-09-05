gg_violin_interaction <- function(data, variable, score, medians_delta, art_con_interaction_default_tidy, patient_label) {
    require(tidyverse)
    require(gghalves)
    require(ggprism)
    require(ggrepel)
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
    if (variable == "cfDNA_cpml") logy <- TRUE else logy <- FALSE
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
        # dplyr::mutate(Followup_pairwise = c("Week24 - Baseline", "Week52 - Baseline", "Week52 - Week24"))
    dplyr::mutate(Followup_pairwise = c("Baseline - Week 24", "Baseline - Week 52", "Week 24 - Week 52"))

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
                        variable == "cfDNA" ~ log2(Week24 / Baseline),
                        TRUE ~ Week24 - Baseline
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
                tidyr::pivot_longer(cols = c(Baseline, Week24, Week52), names_to = "Followup", values_to = "value") %>%
                dplyr::select(Felzartamab, Patient, Followup, delta, delta2),
            by = c("Felzartamab", "Patient", "Followup")
        ) %>%
        dplyr::mutate(
            Followup = Followup %>% factor(levels = c("Baseline", "Week24", "Week52")),
            shape = 21,
            # stroke = 0,
            # shape = ifelse(Patient == 9, 24, 21),
            stroke = ifelse(Patient %in% patient_label, 0.5, 0)
        )
    if (variable %>% stringr::str_detect("cfDNA")) {
        data_patient_labels <- data %>% dplyr::filter(Patient %in% c(patient_label, 9, 19))
    } else {
        data_patient_labels <- data %>% dplyr::filter(Patient %in% patient_label)
    }
    bw <- dplyr::case_when(
        data$variable[[1]] %>% stringr::str_detect("ABMRpm|ggt0|ptcgt0|TCMRt|tgt1|igt1") ~ 0.05,
        data$variable[[1]] %>% stringr::str_detect("PC") ~ 0.2,
        TRUE ~ 0.1
    )
    size_point <- ifelse(variable == "cfDNA", 4, 2.5)
    midpoint <- 0
    if (variable %in% c("KT1", "KT2", "InjPC2_5086Set")) {
        gradient_labels <- c("Worsened", "Improved")
        gradient_labels_hjust <- c(-0.075, 1.55)
        col_low <- "red"
        col_high <- "#00ff00bc"
    } else {
        gradient_labels <- c("Improved", "Worsened")
        col_low <- "#00ff00bc"
        col_high <- "red"
        gradient_labels_hjust <- c(-0.075, 1.85)
    }
    plot <- data %>%
        ggplot2::ggplot(ggplot2::aes(x = Followup, y = value)) +
        gghalves::geom_half_violin(
            inherit.aes = FALSE,
            data = data %>% dplyr::filter(Followup %in% c("Baseline")),
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
            data = data %>% dplyr::filter(Followup %in% c("Baseline")),
            mapping = ggplot2::aes(
                x = Followup,
                y = value,
                fill = delta,
                stroke = stroke
            ),
            col = "black",
            shape = data %>% dplyr::filter(Followup %in% c("Baseline")) %>% dplyr::pull(shape),
            position = ggplot2::position_nudge(x = 0.1),
            size = size_point,
            alpha = 0.75
        ) +
        ggplot2::geom_point(
            inherit.aes = FALSE,
            data = data %>% dplyr::filter(Followup %in% c("Week24")),
            mapping = ggplot2::aes(
                x = Followup,
                y = value,
                fill = delta,
                stroke = stroke
            ),
            col = "black",
            shape = data %>% dplyr::filter(Followup %in% c("Week24")) %>% dplyr::pull(shape),
            size = size_point,
            alpha = 0.75
        ) +
        ggplot2::geom_point(
            inherit.aes = FALSE,
            data = data %>% dplyr::filter(Followup %in% c("Week52")),
            mapping = ggplot2::aes(
                x = Followup,
                y = value,
                fill = delta2,
                stroke = stroke
            ),
            col = "black",
            shape = data %>% dplyr::filter(Followup %in% c("Week52")) %>% dplyr::pull(shape),
            position = ggplot2::position_nudge(x = -0.1),
            size = size_point,
            alpha = 0.75
        ) +
        ggplot2::geom_line(
            data = data %>% dplyr::filter(Followup %in% c("Baseline", "Week24")),
            mapping = aes(
                x = Followup,
                col = delta,
                group = Patient
            ),
            position = ggplot2::position_nudge(
                x = ifelse(data %>%
                    dplyr::filter(Followup %in% c("Baseline", "Week24")) %>%
                    dplyr::pull(Followup) == "Baseline", 0.1, 0)
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
        ggrepel::geom_text_repel(
            seed = 42,
            data = data_patient_labels,
            mapping = ggplot2::aes(label = Patient),
            position = ggplot2::position_nudge(
                x = dplyr::case_when(
                    data_patient_labels$Followup == "Baseline" ~ 0.1,
                    data_patient_labels$Followup == "Week24" ~ 0,
                    data_patient_labels$Followup == "Week52" ~ -0.1
                )
            ),
            size = 4,
            col = "#000000",
            segment.color = "black",
            min.segment.length = 0.01,
            segment.size = 0.25,
            # direction = "y",
            max.overlaps = Inf,
            box.padding = 0.5,
            point.padding = 0.25,
            # nudge_x = dplyr::case_when(
            #     data_patient_labels$Followup == "Baseline" ~ 0.1,
            #     data_patient_labels$Followup == "Week24" ~ 0,
            #     data_patient_labels$Followup == "Week52" ~ -0.1
            # ),
            show.legend = FALSE
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
            fill = "Individual patient response     ",
            parse = TRUE
        ) +
        ggplot2::scale_fill_gradient2(
            low = col_low,
            mid = "grey60",
            high = col_high,
            midpoint = midpoint,
            breaks = c(min(data$delta), max(data$delta)),
            labels = gradient_labels,
            # labels = c("fuck", "you"),
            guide = ggplot2::guide_colorbar(
                title.position = "top",
                barwidth = 20,
                draw.ulim = FALSE,
                draw.llim = FALSE,
                label.hjust = gradient_labels_hjust, # first value is improved, second value worsened (i.e., reverse = TRUE)
                label.vjust = 7.75,
                reverse = TRUE
            )
        ) +
        ggplot2::scale_color_gradient2(
            low = col_low,
            mid = "grey60",
            high = col_high,
            midpoint = midpoint,
            breaks = c(min(data$delta), max(data$delta)),
            labels = gradient_labels,
            # labels = c("fuck", "you"),
            guide = "none"
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
