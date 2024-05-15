gg_violin_interaction <- function(data, variable, score, medians, medians_delta, art_con_interaction_default_tidy) {
    require(tidyverse)
    require(gghalves)
    require(ggprism)




    delta <- medians_delta %>%
        dplyr::mutate(Felzartamab = Felzartamab %>% factor(labels = c("Placebo", "Felzartamab")))

    # delta %>% print()
    delta_delta <- medians_delta %>%
        # dplyr::filter(Followup_pairwise != "Index - FU2") %>%
        dplyr::distinct(Followup_pairwise, .keep_all = TRUE) %>%
        dplyr::mutate(Followup_pairwise = c("FU1 - Index", "FU2 - Index", "FU2 - FU1"))

    delta_delta_p <- art_con_interaction_default_tidy %>%
        # dplyr::filter(Followup_pairwise != "Index - FU2") %>%
        dplyr::select(Followup_pairwise, adj.p.value) %>%
        dplyr::mutate(FDR = ifelse(
            adj.p.value < 0.01,
            formatC(adj.p.value, digits = 0, format = "e"),
            formatC(adj.p.value, digits = 3, format = "f")
        ))

    data <- data %>%
        dplyr::mutate(
            variable = variable,
            Felzartamab = Felzartamab %>% factor(labels = c("Placebo", "Felzartamab"))
        )

    ymax <- data$value %>% max()
    # medians <- data %>%
    #     reframe(median = median(value), .by = c(Group, Felzartamab)) %>%
    #     spread(Group, median) %>%
    #     mutate(delta = FU1 - Index, delta2 = FU2 - FU1)
    data <- data %>%
        dplyr::left_join(
            data %>%
                dplyr::reframe(median = median(value), .by = c(Group, Felzartamab, Patient)) %>%
                rstatix::spread(Group, median) %>%
                dplyr::mutate(delta = FU1 - Index, delta2 = FU2 - FU1) %>%
                tidyr::pivot_longer(cols = c(Index, FU1, FU2), names_to = "Group", values_to = "value") %>%
                dplyr::select(Felzartamab, Patient, Group, delta, delta2),
            by = c("Felzartamab", "Patient", "Group")
        ) %>%
        # left_join(
        #     delta_delta %>%
        #         dplyr::select(
        #             Group,
        #             Felzartamab,
        #             #   median_delta_delta
        #         ),
        #     by = c("Felzartamab", "Group")
        # ) %>%
        dplyr::mutate(
            Group = Group %>% factor(levels = c("Index", "FU1", "FU2")),
            # delta = delta %>% scale() %>% as.vector()
        )
    # bw <- 0.075
    # plot_proto <-
    data %>%
        ggplot2::ggplot(aes(x = Group, y = value)) +
        gghalves::geom_half_violin(
            inherit.aes = FALSE,
            data = data %>% dplyr::filter(Group %in% c("Index")),
            mapping = ggplot2::aes(
                x = Group,
                y = value
            ),
            bw = ifelse(data$variable[[1]] %in% c("ABMRpm", "ggt0", "ptcgt0", "TCMRt", "tgt1", "igt1"), 0.05, 0.1),
            side = c("l"),
            fill = "grey95",
            col = "grey60",
            trim = FALSE,
            scale = "width"
        ) +
        dplyr::geom_violin(
            inherit.aes = FALSE,
            data = data %>% dplyr::filter(Group == c("FU1")),
            mapping = ggplot2::aes(
                x = Group,
                y = value
            ),
            bw = ifelse(data$variable[[1]] %in% c("ABMRpm", "ggt0", "ptcgt0", "TCMRt", "tgt1", "igt1"), 0.05, 0.1),
            fill = "#f2f2f2ec",
            col = "#cdcdcdb8",
            trim = FALSE,
            scale = "width"
        ) +
        gghalves::geom_half_violin(
            inherit.aes = FALSE,
            data = data %>% dplyr::filter(Group %in% c("FU2")),
            mapping = ggplot2::aes(
                x = Group,
                y = value
            ),
            bw = ifelse(data$variable[[1]] %in% c("ABMRpm", "ggt0", "ptcgt0", "TCMRt", "tgt1", "igt1"), 0.05, 0.1),
            side = c("r"),
            fill = "grey95",
            col = "grey60",
            trim = FALSE,
            scale = "width"
        ) +
        ggplot2::geom_point(
            inherit.aes = FALSE,
            data = data %>% dplyr::filter(Group %in% c("Index")),
            mapping = ggplot2::aes(
                x = Group,
                y = value,
                col = delta
            ),
            position = ggplot2::position_nudge(x = 0.1),
            size = 2, alpha = 0.75
        ) +
        ggplot2::geom_point(
            inherit.aes = FALSE,
            data = data %>% dplyr::filter(Group %in% c("FU1")),
            mapping = ggplot2::aes(
                x = Group,
                y = value,
                col = delta
            ),
            size = 2, alpha = 0.75
        ) +
        ggplot2::geom_point(
            inherit.aes = FALSE,
            data = data %>% dplyr::filter(Group %in% c("FU2")),
            mapping = ggplot2::aes(
                x = Group,
                y = value,
                col = delta2
            ),
            position = ggplot2::position_nudge(x = -0.1),
            size = 2, alpha = 0.75
        ) +
        ggplot2::geom_line(
            data = data %>% dplyr::filter(Group %in% c("Index", "FU1")),
            mapping = aes(
                x = Group,
                col = delta,
                group = Patient
            ),
            position = ggplot2::position_nudge(
                x = ifelse(data %>%
                    dplyr::filter(Group %in% c("Index", "FU1")) %>%
                    dplyr::pull(Group) == "Index", 0.1, 0)
            ),
            linewidth = 0.5, alpha = 0.25,
            show.legend = FALSE
        ) +
        ggplot2::geom_line(
            data = data %>% dplyr::filter(Group %in% c("FU1", "FU2")),
            mapping = ggplot2::aes(
                x = Group,
                col = delta2,
                group = Patient
            ),
            position = ggplot2::position_nudge(
                x = ifelse(data %>%
                    dplyr::filter(Group %in% c("FU1", "FU2")) %>%
                    dplyr::pull(Group) == "FU1", 0, -0.1)
            ),
            linewidth = 0.5, alpha = 0.25,
            show.legend = FALSE
        ) +
        ggplot2::stat_summary(
            fun = median, geom = "point", size = 4
        ) +
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
            midpoint = 0,
            breaks = c(min(data$delta), max(data$delta)),
            labels = c("improved", "worsened"),
            guide = ggplot2::guide_colorbar(
                title.position = "top",
                barwidth = 20,
                ticks = FALSE,
                label.hjust = c(2, -0.05),
                label.vjust = 8,
                reverse = TRUE
            ),
        ) +
        # ggplot2::scale_x_discrete(labels = c("Index\n(n = 11)", "FU1\n(n = 11)", "FU2\n(n = 11)")) +
        ggplot2::scale_x_discrete(labels = c("Day0", "Week24", "Week52")) +
        ggplot2::scale_y_continuous(expand = c(0, 0), trans = "identity") +
        ggplot2::coord_cartesian(
            xlim = c(NA, NA),
            ylim = c(
                ifelse(data$variable[[1]] %in% c("ABMRpm", "ggt0", "ptcgt0", "TCMRt", "tgt1", "igt1"), 0, NA),
                ifelse(data$variable[[1]] %in% c("ABMRpm", "ggt0", "ptcgt0", "TCMRt", "tgt1", "igt1"), 1, NA)
            )
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
}
