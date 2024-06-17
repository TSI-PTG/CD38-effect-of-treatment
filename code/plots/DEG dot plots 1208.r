# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(patchwork) # install.packages("patchwork")
# load DEG results
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/gene_tables_limma_1208.RData")


# WRANGLE DATA FOR PLOTTING ####
DEG_plots00 <- gene_tables %>%
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
    dplyr::select(design, geneset, gene_tables) %>%
    nest(.by = "geneset") %>%
    mutate(
        geneset = geneset %>%
            factor(
                levels = c(
                    "ABMR_activity",
                    "IFNG",
                    "NK_ATAGC_U133", "NK_KTB18_RNAseq", "NK_LM22_U133", "NK_L765",
                    "Endothelial"
                ),
                labels = c(
                    "ABMR_activity genes",
                    "IFNG-inducible ABMR activity genes",
                    "NK_ATAGC_U133 genes", "NK_KTB18_RNAseq genes", "NK_LM22_U133 genes", "NK cell expressed ABMR activity genes",
                    "ABMR-associated endothelial genes"
                )
            ),
        data = map(data, unnest, everything())
    ) %>%
    arrange(geneset)

DEG_plots00$geneset[[1]] %>%
    droplevels() %>%
    levels()

col_dn <- "#005eff"
col_up <- "#ff0040"


# MAKE DOT PLOTS ####
DEG_plots <- DEG_plots00 %>%
    mutate(
        plot_long = pmap(
            list(geneset, data),
            function(geneset, data) {
                data <- data %>%
                    dplyr::rename(logFC = `<U+0394><U+0394> logFC`, p = `<U+0394><U+0394> p`) %>%
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
                    ggplot2::geom_vline(xintercept = 0, linetype = "solid", col = "grey10", size = 0.1) +
                    ggplot2::geom_point(col = data$col) +
                    ggplot2::geom_segment(
                        mapping = ggplot2::aes(x = 0, xend = logFC, y = Symb),
                        col = data$col, size = 0.25
                    ) +
                    ggplot2::labs(
                        x = "\u0394\u0394 logFC",
                        y = NULL,
                        # title = paste(geneset %>% str_replace("_", " "), "genes"),
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
        ),
        plot_wide = pmap(
            list(geneset, data),
            function(geneset, data) {
                data <- data %>%
                    dplyr::rename(logFC = `<U+0394><U+0394> logFC`, p = `<U+0394><U+0394> p`) %>%
                    dplyr::mutate(
                        col = dplyr::case_when(
                            p < 0.05 & logFC < 0 ~ col_dn,
                            p < 0.05 & logFC > 0 ~ col_up,
                            TRUE ~ "grey30"
                        )
                    )
                order_symb <- data %>%
                    dplyr::filter(design == "Baseline - Week24") %>%
                    dplyr::arrange(logFC) %>%
                    dplyr::pull(Symb)
                data %>%
                    dplyr::mutate(Symb = Symb %>% factor(levels = order_symb)) %>%
                    ggplot2::ggplot(mapping = ggplot2::aes(x = Symb, y = logFC)) +
                    ggplot2::geom_hline(yintercept = 0, linetype = "solid", col = "grey10", size = 0.1) +
                    ggplot2::geom_point(col = data$col) +
                    ggplot2::geom_segment(
                        mapping = ggplot2::aes(y = 0, yend = logFC, x = Symb),
                        col = data$col, size = 0.25
                    ) +
                    ggplot2::labs(
                        y = "\u0394\u0394 logFC",
                        x = NULL,
                        # title = paste(geneset %>% str_replace("_", " "), "genes"),
                        title = geneset
                    ) +
                    ggplot2::theme_bw() +
                    ggplot2::theme(
                        plot.title = element_text(size = 10),
                        strip.text = element_text(size = 7),
                        panel.grid.minor = element_blank(),
                        panel.grid.major = element_line(size = 0.25),
                        axis.text.y = element_text(colour = "black"),
                        axis.text.x = element_text(colour = "black", size = 7, angle = 90, hjust = 1, vjust = 0.5)
                    ) +
                    ggplot2::facet_wrap(~design, scales = "fixed", nrow = 1)
            }
        )
    )
DEG_plots$plot_long[[2]]


# SAVE THE PLOT DATA ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/"
save(DEG_plots, file = paste(saveDir, "gene_DEG_plots.RData", sep = ""))
