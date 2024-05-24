# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(ggpubr) # install.packages("ggpubr")
library(patchwork) # install.packages("patchwork")
# Custom operators, functions, and datasets
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
# source plot function
# source("C:/slope/CD38-effect-of-treatment/code/functions/plot.gg_violin_interaction.r")
# load reference data
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/felzartamab_cfdna_cor_k1208")



# WRANGLE DATA FOR PLOTTING ####
# tmp <- felzartamab_cfdna_cor_k1208 %>%
#     dplyr::select(category:score, quantile_regression) %>%
#     unnest(
#         c(quantile_regression),
#         names_repair = tidyr_legacy
#     ) %>%
#     rename(
#         slope = rq_slope,
#         p = rq_slope_p
#     ) %>%
#     unnest(data) %>%
#     nest(
#         .by = c(
#             category, annotation, score, variable,
#             Followup_pairwise
#         )
#     )

data_plot <- felzartamab_cfdna_cor_k1208 %>%
    dplyr::select(category:score, quantile_regression, cor) %>%
    unnest(
        c(quantile_regression, cor),
        names_repair = tidyr_legacy
    ) %>%
    unnest(data) %>%
    nest(
        .by = c(
            category, annotation, score, variable,
            Followup_pairwise
        )
    )



# MAKE INDIVIDUAL PLOTS ####
-100 / -0.5

50 / 0.4

felzartamab_cfdna_qr_plots <- data_plot %>%
    mutate(
        plot_scatter = pmap(
            list(data, score),
            function(data, score) {
                qr <- data %>%
                    dplyr::select(Felzartamab, rq_slope, rq_slope_p) %>%
                    mutate(
                        slope = rq_slope %>% formatC(digits = 1, format = "f"),
                        p = dplyr::case_when(
                            rq_slope_p < 0.001 ~ rq_slope_p %>% formatC(digits = 1, format = "e"),
                            TRUE ~ rq_slope_p %>% formatC(digits = 2, format = "f")
                        ),
                        label = paste("slope = ", slope, " | p = ", p, sep = "")
                    ) %>%
                    distinct(Felzartamab, .keep_all = TRUE)
                cor <- data %>%
                    dplyr::select(Felzartamab, cor, p) %>%
                    mutate(
                        cor = cor %>% formatC(digits = 1, format = "f"),
                        p = dplyr::case_when(
                            p < 0.001 ~ p %>% formatC(digits = 1, format = "e"),
                            TRUE ~ p %>% formatC(digits = 2, format = "f")
                        ),
                        label = paste("scc = ", cor, " | p = ", p, sep = "")
                    ) %>%
                    distinct(Felzartamab, .keep_all = TRUE)
                # null_slope <- mean(c(
                #     min(data$delta_score) / min(data$delta_cfdna),
                #     max(data$delta_score) / max(data$delta_cfdna)
                # ))
                null_slope <- mean(c(
                    min(data$delta_cfdna) / min(data$delta_score),
                    max(data$delta_cfdna) / max(data$delta_score)
                ))
                ab_data <- tibble(y = data$delta_score * null_slope, x = data$delta_score)
                plot <- data %>%
                    ggplot2::ggplot(mapping = aes(x = delta_score, y = delta_cfdna)) +
                    ggplot2::geom_hline(yintercept = 0, linetype = "solid", col = "grey70") +
                    ggplot2::geom_vline(xintercept = 0, linetype = "solid", col = "grey70") +
                    ggplot2::geom_line(
                        inherit.aes = FALSE,
                        data = ab_data,
                        mapping = aes(x = x, y = y),
                        linetype = "dotted"
                    ) +
                    ggplot2::geom_quantile(
                        method = "rq",
                        quantiles = c(0.5),
                        formula = y ~ x,
                        colour = "#f52727",
                        linewidth = 1
                    ) +
                    # ggplot2::geom_point(shape = 21, size = 7, fill = "grey95", alpha = 0.5) +
                    ggplot2::geom_point(pch = 21, size = 3, fill = "grey90") +
                    # ggplot2::geom_text(
                    #     mapping = ggplot2::aes(label = Patient),
                    #     size = 4,
                    #     col = "black",
                    #     show.legend = FALSE
                    # ) +
                    ggrepel::geom_text_repel(
                        seed = 42,
                        mapping = ggplot2::aes(label = Patient),
                        size = 4,
                        col = "#000000",
                        segment.color = "black",
                        min.segment.length = 0.01,
                        segment.size = 0.25,
                        # direction = "y",
                        max.overlaps = Inf,
                        box.padding = 0.5,
                        point.padding = 0.35,
                        # nudge_y = 10,
                        show.legend = FALSE
                    ) +
                    ggplot2::geom_text(
                        inherit.aes = FALSE,
                        mapping = ggplot2::aes(
                            x = -Inf, y = Inf
                            
                        ),
                        label = paste("null slope = ", null_slope %>% formatC(digits = 1, format = "f"), sep = ""),
                        vjust = -2.75,
                        hjust = 0,
                        size = 4,
                        fontface = "italic"
                    ) +
                    ggplot2::geom_text(
                        data = qr,
                        mapping = ggplot2::aes(x = -Inf, y = Inf, label = label),
                        vjust = -4.25,
                        hjust = 0,
                        size = 4,
                        fontface = ifelse(as.numeric(qr$p) < 0.05, "bold.italic", "italic")
                    ) +
                    ggplot2::geom_text(
                        data = cor,
                        mapping = ggplot2::aes(x = -Inf, y = Inf, label = label),
                        vjust = -5.75,
                        hjust = 0,
                        size = 4,
                        fontface = ifelse(as.numeric(qr$p) < 0.05, "bold.italic", "italic")
                    ) +
                    ggplot2::labs(
                        x = paste("\u394", score %>% stringr::str_replace("\\(", "\n("), sep = " "),
                        y = "\u394 Donor-derived cell-free DNA\n(dd-cfDNA, cp/mL)",
                        parse = TRUE
                    ) +
                    ggplot2::coord_cartesian(clip = "off") +
                    ggplot2::theme_bw() +
                    ggplot2::theme(
                        plot.margin = unit(c(0.7, 0.2, 0.2, 0.2), "cm"),
                        axis.text = element_text(colour = "black"),
                        strip.text = element_text(size = 12, colour = "black")
                    ) +
                    facet_wrap(~Felzartamab, scales = "fixed")
            }
        )
    )
felzartamab_cfdna_qr_plots$plot_scatter[[1]]



# SAVE THE PLOT DATA ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/"
save(felzartamab_cfdna_qr_plots, file = paste(saveDir, "Felzartamab_cfdna_quantile_regression_plots_1208.RData", sep = ""))
