# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(ggrepel) # install.packages("ggrepel")
library(ggpubr) # install.packages("ggpubr")
library(patchwork) # install.packages("patchwork")
# Custom operators and functions
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
# load volcano plot data
load("results/plots_volcano.RData")
load("results/plots_volcano_timeseries.RData")
load("results/plots_genesets_deg.RData")


# UNIUVERSAL PLOTTING PARAMETERS ####
size_title <- 15


# DEFINE THE GENESET PLOTS ####
df_plots_deg <- plots_deg  %>%
    dplyr::filter(
        geneset %in% c(
            "IFNG-inducible ABMR activity genes",
            "NK cell-expressed ABMR activity genes",
            "ABMR-associated endothelial genes"
        )
    )


# MAKE THE GENESET PLOT PANELS ####
panel_genesets <- wrap_plots(
    df_plots_deg$plot_deg,
    nrow = 1, ncol = 3
) +
    plot_annotation(tag_levels = list(c("D", "E", "F"))) +
    plot_layout(
        guides = "collect",
    ) &
    theme(
        axis.text.y = element_text(size = 6),
        legend.position = "top",
        plot.background = element_rect(fill = "white"),
        plot.tag = element_text(vjust = 1)
    )

panel_genesets <- panel_genesets %>%
    ggarrange(., NULL, ncol = 2, nrow = 1, widths = c(1, 0))


# MAKE THE VOLCANO TIMESERIES PANELS ####
panel_volcano_timeseries <- plot_volcano_timeseries +
    plot_annotation(tag_levels = list(c("C"))) &
    theme(
        legend.position = "none",
        plot.background = element_rect(fill = "white"),
        plot.tag = element_text(vjust = 1)
    )%>% suppressWarnings()

panel_volcano_timeseries2 <- panel_volcano_timeseries %>%
    ggarrange(NULL, ., NULL, ncol = 3, nrow = 1, widths = c(0, 1, 0)) %>% suppressWarnings()



# MAKE THE VOLCANO ENRICHMENT PANELS ####
panel_volcano <- wrap_plots(
    A = plot_volcano$plot_volcano_enrichment[[1]],
    B = plot_volcano$plot_volcano_enrichment[[2]],
    design = "AB"
) +
    plot_annotation(tag_levels = list(c("A", "B"))) +
    plot_layout(
    ) &
    theme(
        legend.position = "none",
        plot.background = element_rect(fill = "white", colour = "white"),
        plot.tag = element_text(vjust = 1, face = "plain")
    )

panel_volcano <- panel_volcano %>%
    ggarrange(., NULL, ncol = 2, nrow = 1, widths = c(1, 0))


# ARRANGE PANELS ####
figure3 <- ggarrange(
    ggarrange(
        panel_volcano,
        panel_volcano_timeseries2,
        nrow = 1, widths = c(1, 0.8)
    ),
    panel_genesets,
    nrow = 2, ncol = 1,
    heights = c(1, 0.8)
)


# SAVE THE PLOT ####
saveDir <- "png/"
ggsave(
    figure3,
    file = paste0(saveDir, "figure3.png"),
    dpi = 300,
    width = 30,
    height = 15,
    units = "cm",
    bg = "white"
) %>% suppressWarnings()
