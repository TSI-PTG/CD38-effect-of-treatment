# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(ggrepel) # install.packages("ggrepel")
library(ggpubr) # install.packages("ggpubr")
library(patchwork) # install.packages("patchwork")
library(readxl) # install.packages("readxl")
# Custom operators, functions, and datasets
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
# load reference data
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/all probes limma 1208.RData")
# load SCC data
simplefile <- read_excel("Z:/MISC/Phil/AA All papers in progress/A GC papers/0000 simple XL files/Kidney 5086/MASTER COPY K5086 SimpleCorrAAInjRej 5AAInjNR 7AARej.xlsx")



# DEFINE PROBE CORRELATIONS WITH EABMR ####
probes_eabmr <- simplefile %>%
    dplyr::select(Affy, "pvalRej7AA4-EABMR") %>%
    rename(
        AffyID = Affy,
        EABMRcorrp = `pvalRej7AA4-EABMR`
    ) %>%
    mutate(ABMRrank = EABMRcorrp %>% rank()) %>%
    arrange(EABMRcorrp)

selected <- simplefile %>%
    dplyr::select(Affy, SYMB, "pvalRej7AA4-EABMR") %>%
    rename(
        AffyID = Affy,
        EABMRcorrp = `pvalRej7AA4-EABMR`
    ) %>%
    arrange(EABMRcorrp) %>%
    distinct(SYMB, .keep_all = TRUE) %>%
    slice(1:20) %>%
    pull(AffyID)


# IDENTIY GENES SUPPRESSED FROM BASELINE TO WEEK24 ####
probes_supressed <- limma_tables %>%
    dplyr::filter(design == "Baseline_vs_Week24") %>%
    pull(toptable) %>%
    pluck(1) %>%
    dplyr::filter(logFC < 0 & P.Value < 0.05) %>%
    # dplyr::filter(logFC < 0) %>%
    pull(AffyID)

probes_increased <- limma_tables %>%
    dplyr::filter(design == "Baseline_vs_Week24") %>%
    pull(toptable) %>%
    pluck(1) %>%
    dplyr::filter(logFC > 0 & P.Value < 0.05) %>%
    # dplyr::filter(logFC < 0) %>%
    pull(AffyID)


# DEFINE COLORS FOR TOP PROBES ####
col_up <- "#7accf9"
col_dn <- "#82fcbf"
col_null <- "grey20"


# WRANGLE THE DATA FOR PLOTTING ####
data_plot <- limma_tables %>%
    dplyr::select(design, table) %>%
    dplyr::filter(design != "Baseline_vs_Week52") %>%
    unnest(everything()) %>%
    left_join(probes_eabmr, by = "AffyID") %>%
    mutate(
        ABMRrank = EABMRcorrp %>% rank(),
        col = case_when(
            AffyID %in% probes_supressed ~ col_up,
            AffyID %in% probes_increased ~ col_dn,
            TRUE ~ col_null
        ),
        alpha = ifelse(col == col_null, 0.5, 0.8),
        order = col %>% factor(levels = c(col_null, col_up, col_dn))
    ) %>%
    arrange(order)


# MAKE VOLCANO PLOT ####
plot_volcano <- data_plot %>%
    ggplot(aes(x = -log10(P.Value), y = logFC)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed") +
    geom_point(shape = 21, fill = data_plot$col, alpha = data_plot$alpha, stroke = 0.125) +
    # geom_point(
    #     # aes(col = direction),
    #     data = data_plot %>% dplyr::filter(AffyID %in% selected),
    #     col = "blue",
    #     size = 2
    # ) +
    # geom_label_repel(
    #     mapping = aes(label = Symb),
    #     size = 2,
    #     max.overlaps = Inf,
    #     min.segment.length = 0.1,
    #     data = data_plot %>% dplyr::filter(AffyID %in% selected),
    # ) +
    labs(
        y = "Felzartamab effect (log2 fold change)",
        x = "-log10 pvalue",
        col = NULL
    ) +
    theme_bw() +
    theme(
        panel.grid = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 15),
        axis.text = element_text(colour = "black")
    ) +
    facet_wrap(~design)


# SAVE THE PLOT ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/output/"
# ggsave(
#     plot_volcano,
#     file = paste(saveDir, "volcano plot all probes 1208.png", sep = ""),
#     dpi = 300,
#     width = 20,
#     height = 12,
#     units = "cm",
#     bg = "white"
# )



# MAKE VOLCANO PLOT ###
min_p <- data_plot %>%
    slice_min(P.Value) %>%
    pull(P.Value) %>%
    log10() * -1.2

data_plot2 <- data_plot %>%
    mutate(
        p = ifelse(design == "Week24_vs_Week52", -log10(P.Value) + min_p, -log10(P.Value))
    )


data_plot2_trimmed <- data_plot2 %>%
    dplyr::filter(AffyID %in% c(probes_supressed, probes_increased)) %>%
    dplyr::select(AffyID, p, logFC, col, alpha) %>%
    arrange(AffyID)


plot_volcano2 <- data_plot2 %>%
    ggplot(aes(x = p, y = logFC)) +
    geom_vline(xintercept = 3.70, linetype = "solid", linewidth = 0.2) +
    geom_vline(xintercept = 3.85, linetype = "solid", linewidth = 0.2) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed") +
    geom_vline(xintercept = -log10(0.05) + min_p, linetype = "dashed") +
    geom_segment(
        inherit.aes = FALSE,
        mapping = aes(x = -Inf, xend = 3.70, y = -Inf),
        linewidth = 0.2
    ) +
    geom_segment(
        inherit.aes = FALSE,
        mapping = aes(x = 3.85, xend = Inf, y = -Inf),
        linewidth = 0.2
    ) +
    geom_segment(
        inherit.aes = FALSE,
        mapping = aes(x = -Inf, xend = 3.70, y = Inf),
        linewidth = 0.2
    ) +
    geom_segment(
        inherit.aes = FALSE,
        mapping = aes(x = 3.85, xend = Inf, y = Inf),
        linewidth = 0.2
    ) +
    geom_segment(
        inherit.aes = FALSE,
        mapping = aes(x = Inf, y = -Inf, yend = Inf),
        linewidth = 0.2
    ) +
    geom_point(
        data = data_plot2 %>% dplyr::filter(AffyID %nin% c(probes_supressed, probes_increased)),
        col = col_null,
        alpha = 0.25
    ) +
    geom_line(
        data = data_plot2_trimmed,
        mapping = aes(x = p, y = logFC, group = AffyID),
        col = data_plot2_trimmed$col,
        alpha = 0.125,
        linewidth = 0.5
    ) +
    geom_point(
        data = data_plot2_trimmed,
        shape = 21,
        fill = data_plot2_trimmed$col,
        alpha = data_plot2_trimmed$alpha,
        stroke = 0.125
    ) +
    geom_text(
        x = -Inf,
        y = Inf,
        vjust = -0.25,
        hjust = 0,
        label = "\u394\u394 Baseline - Week24",
        fontface = "bold.italic"
    ) +
    geom_text(
        x = 3.85,
        y = Inf,
        vjust = -0.25,
        hjust = 0,
        label = "\u394\u394 Week24 - Week52",
        fontface = "bold.italic"
    ) +
    scale_x_continuous(
        breaks = c(0, 1, 2, 3, min_p, 1 + min_p, 2 + min_p, 3 + min_p),
        labels = 10^-c(0, 1, 2, 3, 0, 1, 2, 3)
    ) +
    scale_y_continuous(expand = c(0, 0.05)) +
    coord_cartesian(clip = "off") +
    labs(
        y = "Treatment effect (log2 fold change)",
        x = "Treatment effect p-value",
        col = NULL
    ) +
    theme_bw() +
    theme(
        panel.border = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_line(linewidth = 0.2),
        panel.grid = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 12),
        axis.text = element_text(colour = "black"),
        plot.margin = unit(c(0.5, 0.1, 0.1, 0.1), "cm")
    )


# SAVE THE PLOT ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/output/"
ggsave(
    plot_volcano2,
    file = paste(saveDir, "volcano2 plot all probes 1208.png", sep = ""),
    dpi = 300,
    width = 18,
    height = 12,
    units = "cm",
    bg = "white"
) %>% suppressWarnings()
