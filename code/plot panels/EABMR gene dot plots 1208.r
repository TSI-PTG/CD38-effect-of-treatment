# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(ggpubr) # install.packages("ggpubr")
library(patchwork) # install.packages("patchwork")
# Custom operators, functions, and datasets
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
# load reference data
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/EABMRcorr probes limma 1208.RData")



# WRAGNLE THE DATA FOR PLOTTING ####
data_plot <- limma_tables %>%
    dplyr::select(design, toptable) %>%
    dplyr::filter(design != "Index_vs_Week52") %>%
    unnest(everything())


# MAKE VOLCANO PLOT ###
plot_volcanco <- data_plot %>%
    ggplot(aes(x = -log10(P.Value), y = logFC)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed") +
    geom_point() +
    # geom_point(
    #     aes(col = direction),
    #     data = df_volcano_topprobes,
    #     size = 3
    # ) +
    # geom_label_repel(
    #     aes(label = Symb),
    #     max.overlaps = 40,
    #     min.segment.length = 0.1,
    #     data = df_volcano_topprobes,
    # ) +
    # scale_color_manual(values = cols_topprobes) +
    labs(
        y = "log2 Fold Change",
        x = "-log10 pvalue",
        col = NULL
    ) +
    theme_bw() +
    theme(
        panel.grid = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 15)
    ) +
    facet_grid(~design)


# SAVE THE PLOT ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP2.4 BLAD TBB/Output/"
ggsave(
    plot_volcanco,
    file = paste(saveDir, "BLAD volcano plot (all time - time corrected).png", sep = ""),
    dpi = 300,
    width = 20,
    height = 20,
    units = "cm",
    bg = "white"
)
