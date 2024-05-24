# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(ggpubr) # install.packages("ggpubr")
library(patchwork) # install.packages("patchwork")
# Custom operators, functions, and datasets
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
# source plot function
# source("C:/R/CD38-effect-of-treatment/code/functions/plot.gg_violin_interaction.r")
# load reference data
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/felzartamab_cfdna_cor_k1208")

felzartamab_cfdna_cor_k1208$data[[1]] %>%
    nest(.by = )


tmp <- felzartamab_cfdna_cor_k1208 %>%
    dplyr::select(category:score, cor) %>%
    unnest(cor) %>%
    unnest(data) %>%
    nest(
        .by = c(
            category, annotation, score, variable,
            Followup_pairwise
        )
    )



# WRANGLE DATA FOR PLOTTING ####
tmp <- felzartamab_cfdna_cor_k1208 %>%
    dplyr::select(category:score, cor, cor_foldchange_log2) %>%
    unnest(
        c(cor, cor_foldchange_log2),
        names_repair = tidyr_legacy
    ) %>%
    rename(cor_foldchange_log2 = cor1, p_folcdhange_log2 = p1) %>%
    dplyr::select(-data1) %>%
    unnest(data) %>%
    nest(
        .by = c(
            category, annotation, score, variable,
            Followup_pairwise
        )
    )


data_plot <- felzartamab_cfdna_cor_k1208 %>%
    dplyr::select(category:score, cor, cor_foldchange_log2) %>%
    unnest(
        c(cor, cor_foldchange_log2),
        names_repair = tidyr_legacy
    ) %>%
    rename(cor_foldchange_log2 = cor1, p_foldchange_log2 = p1) %>%
    dplyr::select(-data1) %>%
    unnest(data) %>%
    nest(
        .by = c(
            category, annotation, score, variable,
            Followup_pairwise
        )
    )

tmp <- data_plot$data[[1]] %>%
    dplyr::filter(Felzartamab == "Felzartamab")

cor(
    tmp$delta_foldchange_log2_cfdna,
    tmp$delta_score,
    use = "p", method = "spearman"
)

cor(
    tmp$delta_cfdna,
    tmp$delta_score,
    use = "p", method = "spearman"
)


# MAKE INDIVIDUAL PLOTS ####
felzartamab_cfdna_cor_plots <- data_plot %>%
    mutate(
        plot_scatter = pmap(
            list(data, score),
            function(data, score) {
                cor <- data %>%
                    dplyr::select(Felzartamab, cor, p) %>%
                    mutate(
                        cor = cor %>% round(2),
                        p = dplyr::case_when(
                            p < 0.001 ~ p %>% formatC(digits = 1, format = "e"),
                            TRUE ~ p %>% formatC(digits = 2, format = "f")
                        ),
                        label = paste("SCC = ", cor, " | p = ", p, sep = "")
                    ) %>%
                    distinct(Felzartamab, .keep_all = TRUE)
                data %>%
                    ggplot2::ggplot(mapping = aes(x = delta_score, y = delta_cfdna)) +
                    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", col = "grey30") +
                    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", col = "grey30") +
                    # ggplot2::geom_smooth(
                    #     method = "lm",
                    #     formula = y ~ x,
                    #     se = FALSE,
                    #     col = "#0c3571"
                    # ) +
                    ggplot2::geom_quantile(quantiles = c(0.5), formula = y ~ x, colour = "red") +
                    ggplot2::geom_point(shape = 21, size = 7, fill = "grey95") +
                    ggplot2::geom_text(
                        mapping = ggplot2::aes(label = Patient),
                        size = 4,
                        col = "black",
                        show.legend = FALSE
                    ) +
                    ggplot2::geom_text(
                        data = cor,
                        mapping = ggplot2::aes(x = -Inf, y = Inf, label = label),
                        vjust = -2.75,
                        hjust = 0,
                        size = 4,
                        fontface = ifelse(as.numeric(cor$p) < 0.05, "bold.italic", "italic")
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
        ),
        plot_scatter_foldchange = pmap(
            list(data, score),
            function(data, score) {
                cor <- data %>%
                    dplyr::select(Felzartamab, cor_foldchange_log2, p_foldchange_log2) %>%
                    mutate(
                        cor = cor_foldchange_log2 %>% round(2),
                        p = dplyr::case_when(
                            p_foldchange_log2 < 0.001 ~ p_foldchange_log2 %>% formatC(digits = 1, format = "e"),
                            TRUE ~ p_foldchange_log2 %>% formatC(digits = 2, format = "f")
                        ),
                        label = paste("SCC = ", cor, " | p = ", p, sep = "")
                    ) %>%
                    distinct(Felzartamab, .keep_all = TRUE)
                data %>%
                    ggplot2::ggplot(mapping = aes(x = delta_score, y = delta_foldchange_log2_cfdna)) +
                    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", col = "grey30") +
                    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", col = "grey30") +
                    # ggplot2::geom_smooth(
                    #     method = "lm",
                    #     formula = y ~ x,
                    #     se = FALSE,
                    #     col = "#0c3571"
                    # ) +
                    ggplot2::geom_quantile(quantiles = c(0.5), formula = y ~ x, colour = "red") +
                    ggplot2::geom_point(shape = 21, size = 7, fill = "grey95") +
                    ggplot2::geom_text(
                        mapping = ggplot2::aes(label = Patient),
                        size = 4,
                        col = "black",
                        show.legend = FALSE
                    ) +
                    ggplot2::geom_text(
                        data = cor,
                        mapping = ggplot2::aes(x = -Inf, y = Inf, label = label),
                        vjust = -2.75,
                        hjust = 0,
                        size = 4,
                        fontface = ifelse(as.numeric(cor$p) < 0.05, "bold.italic", "italic")
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
felzartamab_cfdna_cor_plots$plot_scatter[[1]]
felzartamab_cfdna_cor_plots$plot_scatter_foldchange[[2]]



# SAVE THE PLOT DATA ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/"
save(felzartamab_cfdna_cor_plots, file = paste(saveDir, "Felzartamab_cfdna_cor_plots_plots_1208.RData", sep = ""))
