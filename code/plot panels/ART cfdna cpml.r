# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(ggpubr) # install.packages("ggpubr")
library(patchwork) # install.packages("patchwork")
library(ggprism) # install.packages("ggprism")
# Custom operators, functions, and datasets
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
# load reference data
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/Felzartamab_plots.RData")
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/Felzartamab_cfdna_percent_quantile_regression_plots_1208_wholepopulation.RData")


# EXTRACT LEGEND FOR PLOTS ####
legend <- felzartamab_plots %>%
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
plot_cfdna <- felzartamab_plots %>%
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
    # plot_annotation(
    #     tag_levels = list(c(LETTERS[1:15]))
    # ) &
    theme(
        # legend.position = "none",
        axis.text = element_text(size = 12, colour = "black"),
        plot.tag = element_text(size = 20, face = "bold.italic", vjust = 1),
        strip.text = element_text(size = 20)

    )

panels <- ggarrange(
    panel_legend,
    panel_violin,
    nrow = 2, ncol = 1,
    heights = c(0.15, 1)
)

# panels_legend <- ggarrange(
#     panel_legend,
#     panels,
#     nrow = 2, ncol = 1,
#     heights = c(0.075, 1)
# )


# SAVE THE PLOTS ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/output/"
ggsave(
    filename = paste(saveDir, "Felzartamab cfDNA cpml panel.png"),
    plot = panels,
    dpi = 600,
    width = 20,
    height = 16,
    units = "cm",
    bg = "white"
)


# END ####
