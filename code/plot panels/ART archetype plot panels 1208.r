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


# EXTRACT LEGEND FOR PLOTS ####
panel_legend <- felzartamab_plots %>%
    dplyr::filter(variable == "AMAT1") %>%
    pull(plot_violin) %>%
    ggpubr::get_legend() %>%
    ggpubr::as_ggplot() +
    theme(plot.margin = unit(c(0, 0, -1, 0), "cm"))


# MAKE PANEL OF VIOLIN PLOTS ####
panel_violin <- felzartamab_plots %>%
    dplyr::filter(category %in% c("archetypes")) %>%
    dplyr::filter(variable %in% c("RejAA_NR", "RejAA_FABMR", "RejAA_TCMR1"))  %>% 
    pull(plot_violin) %>%
    wrap_plots(nrow = 1, ncol = 3) +
    plot_annotation(tag_levels = list(c(LETTERS[1:15]))) &
    theme(
        legend.position = "none",
        axis.text = element_text(size = 10, colour = "black"), plot.tag = element_text(size = 20, face = "bold", vjust = 1)
    )

panels_violin_legend <- (panel_legend / panel_violin) +
    plot_layout(
        nrow = 2,
        heights = c(0.5, 1.5)
    )



# SAVE THE PLOTS ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/output/"
ggsave(
    filename = paste(saveDir, "Felzartamab violin plots archetype scores.png"),
    plot = panels_violin_legend,
    dpi = 300,
    width = 40,
    height = 11,
    units = "cm",
    bg = "white"
)
# ggsave(
#     filename = paste(saveDir, "Felzartamab patient plots ABMR and TCMR scores.png"),
#     plot = panel_patient,
#     dpi = 300,
#     width = 60,
#     height = 22,
#     units = "cm",
#     bg = "white"
# )


# END ####
