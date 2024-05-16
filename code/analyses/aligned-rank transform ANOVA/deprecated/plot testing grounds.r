medians <- df_univariate_02 %>%
    dplyr::slice(1) %>%
    pull(medians) %>%
    pluck(1)

variable <- df_univariate_02 %>%
    dplyr::slice(1) %>%
    pull(variable)


df_univariate_02 %>%
    dplyr::slice(1) %>%
    pull(medians_delta) %>%
    pluck(1) %>%
    dplyr::filter(Group_pairwise != "Index - FU2") %>% 
    distinct(Group_pairwise, .keep_all = TRUE)


data <- df_univariate_02 %>%
    dplyr::slice(1) %>%
    pull(data) %>%
    pluck(1) %>%
    mutate(
        variable = variable,
        Felz = Felz %>% factor(labels = c("No Felzartamab", "Felzartamab"))
    ) %>%
    left_join(
        data %>%
            reframe(median = median(value), .by = c(Group, Felz, Patient)) %>%
            spread(Group, median) %>%
            mutate(delta = FU1 - Index, delta2 = FU2 - FU1) %>%
            pivot_longer(cols = c(Index, FU1, FU2), names_to = "Group", values_to = "value") %>%
            # mutate(delta = ifelse(delta < 0, "improved", "worsened")) %>%
            dplyr::select(Felz, Patient, Group, delta, delta2),
        by = c("Felz", "Patient", "Group")
    ) %>%
    mutate(
        Group = Group %>% factor(levels = c("Index", "FU1", "FU2")),
    )

data %>%
    ggplot(
        aes(
            x = Group,
            y = value,
            # group = Patient,
            # linetype = Felz
        )
    ) +
    geom_half_violin(
        inherit.aes = FALSE,
        data = data %>% dplyr::filter(Group %in% c("Index")),
        mapping = aes(
            x = Group,
            y = value
        ),
        # draw_quantiles = c(0.25, 0.75),
        side = c("l"),
        fill = "grey95",
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
        fill = "grey95",
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
        # draw_quantiles = c(0.25, 0.75),
        side = c("r"),
        fill = "grey95",
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
            col = delta,
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
    facet_wrap(~Felz)



range_colorbar <- df_univariate_01 %>%
    dplyr::select(medians_delta) %>%
    unnest(medians_delta) %>%
    pull(median_delta) %>%
    range()

panel_pairs <- plot_violin_pairs %>%
    dplyr::select(category, variable, gg_line) %>%
    mutate(
        gg_line = pmap(
            list(variable, gg_line),
            function(variable, gg_line) {
                if (variable == "igt1") {
                    gg_line +
                        scale_color_gradient2(
                            low = "#00ff00bc",
                            mid = "grey60",
                            high = "red",
                            midpoint = 0,
                            limits = range_colorbar,
                            breaks = range_colorbar,
                            labels = c("improved", "worsened"),
                            guide = guide_colorbar(
                                title.position = "top",
                                barwidth = 20,
                                ticks = FALSE,
                                label.hjust = c(1.1, -0.1),
                                label.vjust = 8,
                                reverse = TRUE
                            )
                        )
                } else if (variable != "igt1") {
                    gg_line +
                        scale_color_gradient2(guide = NULL)
                    # theme(legend.position = "none")
                }
            }
        )
    )






panel_pairs <- wrap_plots(
    panel_pairs %>%
        dplyr::filter(category %in% c("ABMR", "TCMR")) %>%
        pull(gg_line),
    nrow = 2,
    # guides = "collect"
) +
    plot_layout(guides = "collect") &
    theme(legend.position = "top")
