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



# DEFINE PROBE CORRELATIONS WITH EABMR
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


# IDENTIY GENES SUPPRESSED FROM BASELINE TO WEEK24
probes_supressed <- limma_tables %>%
    dplyr::filter(design == "Index_vs_Week24") %>%
    pull(toptable) %>%
    pluck(1) %>%
    dplyr::filter(logFC < 0 & P.Value < 0.05) %>%
    # dplyr::filter(logFC < 0) %>%
    pull(AffyID)

probes_increased <- limma_tables %>%
    dplyr::filter(design == "Index_vs_Week24") %>%
    pull(toptable) %>%
    pluck(1) %>%
    dplyr::filter(logFC > 0 & P.Value < 0.05) %>%
    # dplyr::filter(logFC < 0) %>%
    pull(AffyID)



# WRANGLE THE DATA FOR PLOTTING ####
data_plot <- limma_tables %>%
    dplyr::select(design, table) %>%
    dplyr::filter(design != "Index_vs_Week52") %>%
    unnest(everything()) %>%
    left_join(probes_eabmr, by = "AffyID") %>%
    mutate(
        ABMRrank = EABMRcorrp %>% rank(),
        col = case_when(
            AffyID %in% probes_supressed ~ "blue",
            AffyID %in% probes_increased ~ "red",
            TRUE ~ "black"
        ),
        order = col %>% factor(levels = c("black", "red", "blue"))
    ) %>%
    arrange(order)


# DEFINE COLORS FOR TOP PROBES ####
# cols_topprobes <- c("up" = "red", "down" = "blue")


# MAKE VOLCANO PLOT ###
plot_volcano <- data_plot %>%
    ggplot(aes(x = -log10(P.Value), y = logFC)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed") +
    geom_point(col = data_plot$col) +
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
    scale_color_gradient(
        high = "blue",
        low = "black",
    ) +
    labs(
        y = "Felzartamab effect\n(log2 foldchange)",
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
ggsave(
    plot_volcano,
    file = paste(saveDir, "volcano plot all probes 1208.png", sep = ""),
    dpi = 300,
    width = 20,
    height = 12,
    units = "cm",
    bg = "white"
)
