# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(ggrepel) # install.packages("ggrepel")
library(ggpubr) # install.packages("ggpubr")
library(patchwork) # install.packages("patchwork")
library(readxl) # install.packages("readxl")
# Custom operators, functions, and datasets
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
source("C:/R/CD38-effect-of-treatment/code/functions/plot.gg_volcano_timeseries.r")
# load reference data
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/IQR_filtered_probes_unique_genes_baseline_corrected_cortex_corrected_limma_1208.RData")


# WRANGLE THE DATA FOR PLOTTING TIMESERIES ####
data_plot_timeseries_00 <- limma_tables %>%
    dplyr::select(design, table) %>%
    dplyr::filter(design != "Baseline_vs_Week52") %>%
    unnest(everything())

min_p <- data_plot_timeseries_00 %>%
    slice_min(`<U+0394><U+0394> p`) %>%
    pull(`<U+0394><U+0394> p`) %>%
    log10() * -1.2

data_plot_timeseries <- data_plot_timeseries_00 %>%
    mutate(
        p = ifelse(
            design == "Week24_vs_Week52",
            -log10(`<U+0394><U+0394> p`) + min_p,
            -log10(`<U+0394><U+0394> p`)
        ),
        logFC = `<U+0394><U+0394> logFC`
    )
# data_plot_timeseries %>% dplyr::filter(AffyID %in% probes_abmr)


# MAKE TIME SERIES VOLCANO PLOT ####
plot_volcano_timeseries <- data_plot_timeseries %>%
    gg_volcano_timeseries(
        xbreak1 = 3.15,
        xbreak2 = 3.30,
        point_size_null = 1.25,
        point_size = 2.5,
        # labels_probes = probes_abmr,
        # labels_probes_n = 100,
        labels_probes_size = 3,
        axis_label_position_x_1 = 1.5
        # col_dn = "#005eff",
        # col_up = "#ff0040",
        # col_null = "grey30"
    )


# SAVE THE PLOT DATA ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/"
save(plot_volcano_timeseries, file = paste(saveDir, "Volcano timeseries plots IQR_filtered_probes_unique_genes_baseline_corrected_cortex_corrected 1208.RData", sep = ""))
