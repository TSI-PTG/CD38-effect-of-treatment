# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(ggrepel) # install.packages("ggrepel")
library(ggpubr) # install.packages("ggpubr")
library(patchwork) # install.packages("patchwork")
library(readxl) # install.packages("readxl")
library(ggsci) # install.packages("ggsci")
# Custom operators and functions
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
source("code/functions/plot.gg_volcano_timeseries.r")
# load DE results
load("results/deg.RData")


# WRANGLE THE DATA FOR PLOTTING TIMESERIES ####
data_plot_timeseries_00 <- deg %>%
    dplyr::select(design, table) %>%
    dplyr::filter(design != "Baseline_vs_Week52") %>%
    unnest(everything())

min_p <- data_plot_timeseries_00 %>%
    slice_min(p) %>%
    pull(p) %>%
    log10() * -1.2

data_plot_timeseries <- data_plot_timeseries_00 %>%
    mutate(
        p = ifelse(
            design == "Week24_vs_Week52",
            -log10(p) + min_p,
            -log10(p)
        )
    )


# UNIVERAL PARAMETERS FOR VOLCANO PLOTS ####
# scales::show_col(pal_npg()(9))
col_dn <- pal_npg()(9)[2]
col_up <- pal_npg()(9)[1]


# MAKE TIME SERIES VOLCANO PLOT ####
plot_volcano_timeseries <- data_plot_timeseries %>%
    gg_volcano_timeseries(
        xbreak1 = 3.15,
        xbreak2 = 3.30,
        point_size_null = 1.25,
        point_size = 2.5,
        labels_probes_size = 3,
        axis_label_position_x_1 = 1.5,
        col_dn = col_dn, col_up = col_up
    )


# SAVE THE PLOT DATA ####
saveDir <- "results/"
save(plot_volcano_timeseries, file = paste(saveDir, "plots_volcano_timeseries.RData", sep = ""))
