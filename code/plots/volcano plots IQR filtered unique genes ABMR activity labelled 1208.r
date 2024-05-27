# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(ggrepel) # install.packages("ggrepel")
library(ggpubr) # install.packages("ggpubr")
library(patchwork) # install.packages("patchwork")
library(readxl) # install.packages("readxl")
# Custom operators, functions, and datasets
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
source("C:/R/CD38-effect-of-treatment/code/functions/plot.gg_volcano.r")
source("C:/R/CD38-effect-of-treatment/code/functions/plot.gg_volcano_timeseries.r")
# load reference data
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/IQR_filtered_probes_unique_genes_limma_1208.RData")
# load SCC data
simplefile <- read_excel("Z:/MISC/Phil/AA All papers in progress/A GC papers/0000 simple XL files/Kidney 5086/MASTER COPY K5086 SimpleCorrAAInjRej 5AAInjNR 7AARej.xlsx")



# DEFINE PROBE CORRELATIONS WITH EABMR ####
genes_abmr <- simplefile %>%
    dplyr::select(Affy, SYMB, "corrRej7AA4-EABMR", "pvalRej7AA4-EABMR") %>%
    rename(
        AffyID = Affy,
        EABMRcorrp = `pvalRej7AA4-EABMR`
    ) %>%
    mutate(ABMRrank = EABMRcorrp %>% rank()) %>%
    dplyr::filter(`corrRej7AA4-EABMR` > 0) %>%
    arrange(EABMRcorrp) %>%
    distinct(SYMB, .keep_all = TRUE) %>%
    slice(1:20) %>%
    pull(SYMB)

probes_abmr <- simplefile %>%
    dplyr::filter(SYMB %in% genes_abmr) %>%
    pull(Affy)



# WRANGLE THE DATA FOR PLOTTING Baseline - Week52 ####
data_plot <- limma_tables %>%
    mutate(
        data_plot = map(
            table,
            function(table) {
                table %>%
                    mutate(
                        p = -log10(`<U+0394><U+0394> p`),
                        logFC = `<U+0394><U+0394> logFC`
                    ) %>%
                    suppressWarnings()
            }
        )
    )



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
data_plot_timeseries %>% dplyr::filter(AffyID %in% probes_abmr)


# MAKE SINGLE VOLCANO PLOTS ####
plot_volcano <- data_plot %>%
    mutate(
        plot_volcano = map2(
            data_plot, design,
            gg_volcano,
            x_break = 1.75,
            labels_probes = probes_abmr,
            point_size_null = 1.25,
            point_size = 2.5,
            labels_probes_size = 3,
        )
    )



plot_volcano$plot_volcano[[3]]



# MAKE TIME SERIES VOLCANO PLOT ####
plot_volcano_timeseries <- data_plot_timeseries %>%
    gg_volcano_timeseries(
        xbreak1 = 3.15,
        xbreak2 = 3.30,
        point_size_null = 1.25,
        point_size = 2.5,
        labels_probes = probes_abmr,
        # labels_probes_n = 100,
        labels_probes_size = 3,
        axis_label_position_x_1 = 1.5
        # col_dn = "#005eff",
        # col_up = "#ff0040",
        # col_null = "grey30"
    )



plot_volcano_timeseries_titled <- ggarrange(
    plot_volcano_timeseries +
        plot_annotation(tag_levels = list("A")) &
        theme(
            plot.tag = element_text(size = 15, face = "bold", vjust = 4),
            plot.background = ggplot2::element_rect(fill = "white", colour = " white")
        )
) %>%
    ggpubr::annotate_figure(
        top = text_grob(
            "Genome-wide Effect of Felzartamab and Relapse Post-Treatment",
            face = "bold.italic", size = 15, hjust = 0.68
        )
    ) %>%
    suppressWarnings()


plot_volcano_titled <- ggarrange(
    plot_volcano$plot_volcano[[3]] +
        plot_annotation(tag_levels = list("B")) &
        theme(plot.tag = element_text(size = 15, face = "bold", vjust = 4))
) %>%
    ggpubr::annotate_figure(
        top = text_grob(
            "Persistent Effect of Felzartamab",
            face = "bold.italic", size = 15, hjust = 0.71
        )
    ) %>%
    suppressWarnings()



plot_volcano_panel <- wrap_plots(
    plot_volcano_timeseries_titled,
    plot_volcano_titled,
    nrow = 1, ncol = 2, widths = c(2, 1.075)
)


# SAVE THE PLOT ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/output/"
ggsave(
    plot_volcano_panel,
    file = paste(saveDir, "volcanos IQR filtered probes unique genes ABMR activity labelled 1208.png", sep = ""),
    dpi = 300,
    width = 35,
    height = 14,
    units = "cm",
    bg = "white"
) %>% suppressWarnings()
