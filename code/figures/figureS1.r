# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(ggpubr) # install.packages("ggpubr")
library(patchwork) # install.packages("patchwork")
library(ggprism) # install.packages("ggprism")
# Custom operators and functions
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
# load plot data
load("results/plots_artANOVA.RData")


# EXTRACT LEGEND FOR PLOTS ####
legend <- plots_artANOVA %>%
    dplyr::filter(variable == "AMAT1") %>%
    pull(plot_violin) %>%
    pluck(1) +
    guides(
        fill = guide_colorbar(
            title.position = "top",
            barwidth = 12,
            draw.ulim = FALSE,
            draw.llim = FALSE,
            label.hjust = c(-0.075, 1.5), # first value is improved, second value worsened (i.e., reverse = TRUE)
            label.vjust = 7.75,
            reverse = TRUE
        )
    ) +
    theme(legend.title = ggplot2::element_text(size = 15, hjust = 0.2, vjust = 1, face = "bold"))

panel_legend <- legend %>%
    ggpubr::get_legend() %>%
    ggpubr::as_ggplot()


# MAKE PANEL OF VIOLIN PLOTS ####
plot_cfdna <- plots_artANOVA %>%
    dplyr::filter(variable == c("cfDNA_cpml")) %>%
    pull(plot_violin) %>%
    pluck(1) +
    facet_wrap(~Felzartamab, nrow = 1, ncol = 2) +
    theme(
        legend.position = "none",
        axis.text = element_text(size = 12, colour = "black"),
        plot.tag = element_text(size = 20, face = "bold.italic", vjust = 1)
    )

panel_violin <- plot_cfdna +
    plot_layout(ncol = 1) +
    theme(
        # legend.position = "none",
        axis.text = element_text(size = 12, colour = "black"),
        plot.tag = element_text(size = 20, face = "bold.italic", vjust = 1),
        strip.text = element_text(size = 20)
    )

figureS1 <- ggarrange(
    panel_legend,
    panel_violin,
    nrow = 2, ncol = 1,
    heights = c(0.15, 1)
)


# SAVE THE PLOTS ####
saveDir <- "png/"
ggsave(
    filename = paste0(saveDir, "figureS1.png"),
    plot = figureS1,
    dpi = 300,
    width = 20,
    height = 16,
    units = "cm",
    bg = "white"
)

# END ####
