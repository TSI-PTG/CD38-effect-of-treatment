# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(patchwork) # install.packages("patchwork")
library(ggsci) # install.packages("ggsci")
# load DE results
load("natmed/data/tables_genesets_abmr.RData")


# WRANGLE DATA FOR PLOTTING ####
tables_genesets_abmr$table_deg[[1]]
plots_deg00 <- tables_genesets_abmr %>%
    mutate(
        design = design %>%
            factor(
                levels = c(
                    "Baseline_vs_Week24",
                    "Week24_vs_Week52",
                    "Baseline_vs_Week52"
                ),
                labels = c(
                    "Baseline - Week24",
                    "Week24 - Week52",
                    "Baseline - Week52"
                )
            )
    ) %>%
    dplyr::select(design, geneset, table_deg) %>%
    nest(.by = "geneset") %>%
    mutate(
        geneset = geneset %>%
            factor(
                levels = c(
                    "ABMR_activity",
                    "ABMR_IFNG",
                    "ABMR_NK",
                    "ABMR_endothelial"
                ),
                labels = c(
                    "ABMR activity genes",
                    "IFNG-inducible ABMR activity genes",
                    "NK cell-expressed ABMR activity genes",
                    "ABMR-associated endothelial genes"
                )
            ),
        data = map(data, unnest, everything())
    ) %>%
    arrange(geneset)


# UNIVERSAL PLOTTING PARAMETERS ####
# scales::show_col(pal_npg()(9))
col_dn <- pal_npg()(9)[2]
col_up <- pal_npg()(9)[1]


# MAKE DOT PLOTS ####
plots_deg <- plots_deg00 %>%
    mutate(
        plot_deg = pmap(
            list(geneset, data),
            function(geneset, data) {
                data <- data %>%
                    dplyr::mutate(
                        col = dplyr::case_when(
                            p < 0.05 & logFC < 0 ~ col_dn,
                            p < 0.05 & logFC > 0 ~ col_up,
                            TRUE ~ "grey30"
                        )
                    )
                order_symb <- data %>%
                    dplyr::filter(design == "Baseline - Week24") %>%
                    dplyr::arrange(logFC %>% desc()) %>%
                    dplyr::pull(Symb)
                data %>%
                    dplyr::mutate(Symb = Symb %>% factor(levels = order_symb)) %>%
                    ggplot2::ggplot(mapping = ggplot2::aes(y = Symb, x = logFC)) +
                    ggplot2::geom_vline(xintercept = 0, linetype = "solid", col = "grey10", linewidth = 0.1) +
                    ggplot2::geom_point(col = data$col) +
                    ggplot2::geom_errorbar(
                        mapping = ggplot2::aes(xmin = logFC - se, xmax = logFC + se),
                        col = data$col,
                        linewidth = 0.25,
                        width = 0.4
                    ) +
                    ggplot2::geom_segment(
                        mapping = ggplot2::aes(x = 0, xend = logFC, y = Symb),
                        col = data$col, linewidth = 0.25
                    ) +
                    ggplot2::labs(
                        x = "\u0394\u0394 logFC",
                        y = NULL,
                        title = geneset
                    ) +
                    ggplot2::coord_cartesian(xlim = c(-1, 1)) +
                    ggplot2::theme_bw() +
                    ggplot2::theme(
                        plot.title = element_text(size = 10),
                        strip.text = element_text(size = 7),
                        panel.grid.minor = element_blank(),
                        panel.grid.major = element_line(size = 0.25),
                        axis.text.x = element_text(colour = "black", size = 6),
                        axis.text.y = element_text(colour = "black", size = 7)
                    ) +
                    ggplot2::facet_wrap(~design, scales = "fixed", nrow = 1)
            }
        )
    )


# SAVE THE PLOT DATA ####
saveDir <- "natmed/results/"
save(plots_deg, file = paste0(saveDir, "plots_genesets_deg.RData"))
