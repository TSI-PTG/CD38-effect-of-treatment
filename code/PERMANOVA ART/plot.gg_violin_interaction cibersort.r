gg_violin_interaction <- function(data, annotation, variable, score, medians, medians_delta, art_con_interaction_default_tidy) {
    delta <- medians_delta %>%
        mutate(
            Felz = Felz %>% factor(labels = c("No Felzartamab", "Felzartamab"))
        )

    # delta %>% print()
    delta_delta <- medians_delta %>%
        # dplyr::filter(Group_pairwise != "Index - FU2") %>%
        distinct(Group_pairwise, .keep_all = TRUE) %>%
        mutate(Group_pairwise = c("FU1 - Index", "FU2 - Index", "FU2 - FU1"))

    delta_delta_p <- art_con_interaction_default_tidy %>%
        # dplyr::filter(Group_pairwise != "Index - FU2") %>%
        dplyr::select(Group_pairwise, adj.p.value) %>%
        mutate(FDR = ifelse(
            adj.p.value < 0.01,
            formatC(adj.p.value, digits = 0, format = "e"),
            formatC(adj.p.value, digits = 3, format = "f")
        ))

    data <- data %>%
        mutate(
            variable = variable,
            Felz = Felz %>% factor(labels = c("No Felzartamab", "Felzartamab"))
        )

    ymax <- data$value %>% max()
    # medians <- data %>%
    #     reframe(median = median(value), .by = c(Group, Felz)) %>%
    #     spread(Group, median) %>%
    #     mutate(delta = FU1 - Index, delta2 = FU2 - FU1)
    data <- data %>%
        left_join(
            data %>%
                reframe(median = median(value), .by = c(Group, Felz, Patient)) %>%
                spread(Group, median) %>%
                mutate(delta = FU1 - Index, delta2 = FU2 - FU1) %>%
                pivot_longer(cols = c(Index, FU1, FU2), names_to = "Group", values_to = "value") %>%
                dplyr::select(Felz, Patient, Group, delta, delta2),
            by = c("Felz", "Patient", "Group")
        ) %>%
        # left_join(
        #     delta_delta %>%
        #         dplyr::select(
        #             Group,
        #             Felz,
        #             #   median_delta_delta
        #         ),
        #     by = c("Felz", "Group")
        # ) %>%
        mutate(
            Group = Group %>% factor(levels = c("Index", "FU1", "FU2")),
            # delta = delta %>% scale() %>% as.vector()
        )
    data %>%
        ggplot(aes(x = Group, y = value)) +
        geom_half_violin(
            inherit.aes = FALSE,
            data = data %>% dplyr::filter(Group %in% c("Index")),
            mapping = aes(
                x = Group,
                y = value
            ),
            # bw = ifelse(data$variable[[1]] %in% c("ABMRpm", "ggt0", "ptcgt0", "TCMRt", "tgt1", "igt1"), 0.05, 0.1),
            side = c("l"),
            fill = "grey95",
            col = "grey60",
            trim = FALSE,
            scale = "width"
        ) +
        geom_violin(
            inherit.aes = FALSE,
            data = data %>% dplyr::filter(Group == c("FU1")),
            mapping = aes(
                x = Group,
                y = value
            ),
            # bw = ifelse(data$variable[[1]] %in% c("ABMRpm", "ggt0", "ptcgt0", "TCMRt", "tgt1", "igt1"), 0.05, 0.1),
            fill = "#f2f2f2ec",
            col = "#cdcdcdb8",
            trim = FALSE,
            scale = "width"
        ) +
        geom_half_violin(
            inherit.aes = FALSE,
            data = data %>% dplyr::filter(Group %in% c("FU2")),
            mapping = aes(
                x = Group,
                y = value
            ),
            # bw = ifelse(data$variable[[1]] %in% c("ABMRpm", "ggt0", "ptcgt0", "TCMRt", "tgt1", "igt1"), 0.05, 0.1),
            side = c("r"),
            fill = "grey95",
            col = "grey60",
            trim = FALSE,
            scale = "width"
        ) +
        geom_point(
            inherit.aes = FALSE,
            data = data %>% dplyr::filter(Group %in% c("Index")),
            mapping = aes(
                x = Group,
                y = value,
                col = delta
            ),
            position = position_nudge(x = 0.1),
            size = 2, alpha = 0.75
        ) +
        geom_point(
            inherit.aes = FALSE,
            data = data %>% dplyr::filter(Group %in% c("FU1")),
            mapping = aes(
                x = Group,
                y = value,
                col = delta
            ),
            size = 2, alpha = 0.75
        ) +
        geom_point(
            inherit.aes = FALSE,
            data = data %>% dplyr::filter(Group %in% c("FU2")),
            mapping = aes(
                x = Group,
                y = value,
                col = delta2
            ),
            position = position_nudge(x = -0.1),
            size = 2, alpha = 0.75
        ) +
        geom_line(
            data = data %>% dplyr::filter(Group %in% c("Index", "FU1")),
            mapping = aes(
                x = Group,
                col = delta,
                group = Patient
            ),
            position = position_nudge(
                x = ifelse(data %>%
                    dplyr::filter(Group %in% c("Index", "FU1")) %>%
                    pull(Group) == "Index", 0.1, 0)
            ),
            linewidth = 0.5, alpha = 0.25,
            show.legend = FALSE
        ) +
        geom_line(
            data = data %>% dplyr::filter(Group %in% c("FU1", "FU2")),
            mapping = aes(
                x = Group,
                col = delta2,
                group = Patient
            ),
            position = position_nudge(
                x = ifelse(data %>%
                    dplyr::filter(Group %in% c("FU1", "FU2")) %>%
                    pull(Group) == "FU1", 0, -0.1)
            ),
            linewidth = 0.5, alpha = 0.25,
            show.legend = FALSE
        ) +
        stat_summary(
            fun = median, geom = "point", size = 4
        ) +
        stat_summary(
            aes(group = Felz),
            fun = median,
            geom = "line",
            linewidth = 1,
            linetype = "dashed"
        ) +
        labs(
            title = paste(
                paste(
                    "\u394\u394", delta_delta$Group_pairwise[[1]], "=",
                    delta_delta$median_delta_delta[[1]] %>% round(2),
                    "| FDR =",
                    delta_delta_p$FDR[[1]]
                ),
                paste(
                    "\u394\u394", delta_delta$Group_pairwise[[3]], "=",
                    delta_delta$median_delta_delta[[3]] %>% round(2),
                    "| FDR =",
                    delta_delta_p$FDR[[3]]
                ),
                paste(
                    "\u394\u394", delta_delta$Group_pairwise[[2]], "=",
                    delta_delta$median_delta_delta[[2]] %>% round(2),
                    "| FDR =",
                    delta_delta_p$FDR[[2]]
                ),
                sep = "\n"
            ),
            x = NULL,
            y = paste(annotation, score, sep = "\n"),
            col = "Individual patient response     ",
            parse = TRUE
        ) +
        scale_color_gradient2(
            low = "#00ff00bc",
            mid = "grey60",
            high = "red",
            midpoint = 0,
            breaks = c(min(data$delta), max(data$delta)),
            labels = c("improved", "worsened"),
            guide = guide_colorbar(
                title.position = "top",
                barwidth = 20,
                ticks = FALSE,
                label.hjust = c(-0.1, 1.1),
                label.vjust = 8,
                reverse = TRUE
            ),
        ) +
        # scale_x_discrete(labels = c("Index\n(n = 11)", "FU1\n(n = 11)", "FU2\n(n = 11)")) +
        scale_x_discrete(labels = c("Index", "FU1", "FU2")) +
        scale_y_continuous(expand = c(0, 0), trans = "identity") +
        coord_cartesian(
            xlim = c(NA, NA),
            ylim = c(NA,NA)
        ) +
        theme_bw() +
        theme(
            axis.title = element_text(size = 15),
            axis.text = element_text(size = 12, colour = "black"),
            panel.grid = element_blank(),
            strip.text = element_text(size = 12, colour = "black"),
            legend.position = "top",
            legend.title = element_text(size = 15, hjust = 0.66, vjust = 1, face = "bold"),
            legend.text = element_text(size = 15),
            plot.title = element_text(
                colour = "black",
                hjust = 0,
                size = 12,
                face = case_when(
                    delta_delta_p$adj.p.value[[1]] < 0.05 ~ "bold.italic",
                    delta_delta_p$adj.p.value[[2]] < 0.05 ~ "bold.italic",
                    TRUE ~ "italic"
                )
            )
        ) +
        facet_wrap(~Felz)
}
