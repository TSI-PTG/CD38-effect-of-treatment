# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(patchwork) # install.packages("patchwork")
# Custom operators, functions, and datasets
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
# load reference data
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/injury_gene_DEG_plots_1208.RData")



# WRANGLE THE DOTPLOT DATA ####
data_plot <- DEG_plots %>%
    mutate(
        keep = map_lgl(
            data,
            function(data) {
                ifelse(nrow(data) > 0, TRUE, FALSE)
            }
        )
    ) %>%
    dplyr::filter(keep)



# MAKE PLOTS PANELS ####
plot_panels_long <- data_plot %>%
    dplyr::filter(cluster  %>% str_detect("New")) %>%
    pull(plot_long) %>%
    wrap_plots(nrow = 3) +
    # plot_layout(axes = "collect") +
    plot_annotation(tag_levels = "A")


# SAVE THE PLOTS ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/output/"
ggsave(
    filename = paste(saveDir, "Injury gene DEG panel long.png"),
    plot = plot_panels_long,
    dpi = 600,
    width = 70,
    height = 40,
    units = "cm",
    bg = "white"
)
