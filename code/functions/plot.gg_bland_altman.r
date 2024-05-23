gg_bland_altman <- function(data, variable, score, id = "Patient", group = "Felzartamab", split = "Followup", value = "value") {
    require(tidyverse)
    data <- data %>%
        dplyr::rename(
            id = id,
            group = group,
            split = split,
            value = value
        ) %>%
        dplyr::select(id, group, split, value) %>%
        tidyr::pivot_wider(names_from = split, values_from = value) %>%
        dplyr::mutate(
            average = (Baseline + Week24) / 2,
            diff = Week24 - Baseline,
            mean_diff = diff %>% mean(),
            lower = mean_diff - 1.96 * sd(diff),
            upper = mean_diff + 1.96 * sd(diff),
            .by = group
        ) %>%
        dplyr::select(
            id, group, Baseline, Week24, Week52,
            average, diff, mean_diff, lower, upper
        )
    data %>%
        ggplot2::ggplot(mapping = ggplot2::aes(x = Baseline, y = diff)) +
        ggplot2::geom_point(pch = 21, size = 7, fill = "grey95") +
        ggplot2::geom_hline(yintercept = 0, linetype = "dotted") +
        ggplot2::geom_hline(mapping = ggplot2::aes(yintercept = mean_diff)) +
        ggplot2::geom_hline(
            mapping = ggplot2::aes(yintercept = lower),
            color = "red", linetype = "dashed"
        ) +
        ggplot2::geom_hline(
            mapping = ggplot2::aes(yintercept = upper),
            color = "red", linetype = "dashed"
        ) +
        ggplot2::geom_text(
            mapping = ggplot2::aes(label = id),
            size = 4,
            col = "black",
            show.legend = FALSE
        ) +
        ggplot2::labs(
            x = paste("Baseline", score %>% stringr::str_replace("\\(", "\n("), sep = " "),
            y = paste("\u394", score %>% stringr::str_replace("\\(", "\n("), sep = " "),
            parse = TRUE
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            # panel.grid = element_blank(),
            axis.text = element_text(colour = "black")
            ) +
        ggplot2::facet_wrap(~group, scales = "free_x")
}
