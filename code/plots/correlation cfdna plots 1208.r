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





# MAKE INDIVIDUAL PLOTS ####
felzartamab_cfdna_cor_plots <- felzartamab_cfdna_cor_k1208 %>%
    dplyr::select(category:score, cor) %>%
    unnest(cor) %>%
    unnest(data) %>%
    nest(
        .by = c(
            category, annotation, score, variable,
            Followup_pairwise
        )
    ) %>%
    mutate(
        plot_scatter = pmap(
            list(data, score),
            function(data, score) {
                data %>%
                    ggplot2::ggplot(mapping = aes(x = delta_score, y = delta_cfdna)) +
                    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", col = "grey30") +
                    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", col = "grey30") +
                    ggplot2::geom_smooth(
                        method = "lm",
                        formula = y ~ x,
                        se = FALSE,
                        col = "#0c3571"
                    ) +
                    ggplot2::geom_point(shape = 21, size = 7, fill = "grey95") +
                    ggplot2::geom_text(
                        mapping = ggplot2::aes(label = Patient),
                        size = 4,
                        col = "black",
                        show.legend = FALSE
                    ) +
                    ggpubr::stat_cor(
                        geom = "label",
                        method = "spearman",
                        cor.coef.name = "SCC",
                        size = 4
                    ) +
                    ggplot2::labs(
                        x = paste("\u394", score %>% stringr::str_replace("\\(", "\n("), sep = " "),
                        y = "\u394 Donor-derived cell-free DNA\n(dd-cfDNA, cp/mL)",
                        parse = TRUE
                    ) +
                    ggplot2::theme_bw() +
                    ggplot2::theme(
                        # panel.grid = element_blank(),
                        axis.text = element_text(colour = "black")
                    ) +
                    facet_grid(~Felzartamab, scales = "free_x")
            }
        )
    )
felzartamab_cfdna_cor_plots$plot_scatter[[1]]


# SAVE THE PLOT DATA ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/"
save(felzartamab_cfdna_cor_plots, file = paste(saveDir, "Felzartamab_cfdna_cor_plots_plots_1208.RData", sep = ""))
