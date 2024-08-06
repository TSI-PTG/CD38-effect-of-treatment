# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(ggrepel) # install.packages("ggrepel")
library(ggpubr) # install.packages("ggpubr")
library(patchwork) # install.packages("patchwork")
library(ggforce) # install.packages("ggforce")
library(readxl) # install.packages("readxl")
library(clusterProfiler) # pak::pak("YuLab-SMU/clusterProfiler")
# Custom operators, functions, and datasets
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
source("C:/R/CD38-effect-of-treatment/code/functions/plot.gg_volcano.r")
# load reference data
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/IQR_filtered_probes_unique_genes_baseline_corrected_cortex_corrected_limma_1208.RData")
# load enrichment results
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/felzartamab_gsea_collated_k1208.RData")
# load simplified enrichment plots
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/simplified_enrichment_plots.RData")


# SUMMARIZE THE GSEA DATA ####
felzartamab_gsea_k1208 %>%
    dplyr::filter(design == "Baseline_vs_Week52") %>% 
    dplyr::select(design, gsea_table) %>%
    unnest(everything()) %>%
    mutate(Symb = core_enrichment %>% strsplit("/"), .by = ID) %>%
    unnest(Symb) %>%
    nest(.by = c("design", "interpretation")) %>%
    mutate(data = map(data, dplyr::distinct, Symb))



# FIND TOP GENES BY OVERLAP ####
overlap <- felzartamab_gsea_k1208 %>%
    dplyr::filter(design == "Baseline_vs_Week52") %>%
    dplyr::select(design, gsea_table) %>%
    unnest(everything()) %>%
    mutate(Symb = core_enrichment %>% strsplit("/"), .by = ID) %>%
    unnest(Symb) %>%
    dplyr::select(design, db, ID, interpretation, Symb) %>%
    nest(.by = c("design")) %>%
    mutate(
        data = map(
            data,
            function(data) {
                data %>%
                    nest(.by = c(Symb, interpretation)) %>%
                    mutate(count = map_dbl(data, nrow)) %>%
                    mutate(prop = count / max(count), .by = interpretation) %>%
                    dplyr::select(-data)
            }
        )
    )

# overlap %>%
#     unnest(data) %>%
#     slice_max(prop, n = 10, by = c("design", "interpretation"), with_ties = FALSE) %>%
#     arrange(interpretation, prop %>% desc()) %>%
#     print(n="all")



# WRANGLE THE GSEA RESULTS ####
gsea <- felzartamab_gsea_k1208 %>%
    dplyr::filter(design == "Baseline_vs_Week52") %>%
    mutate(
        gsea_table = map(
            gsea_table,
            function(gsea_table) {
                gsea_table %>%
                    mutate(Symb = core_enrichment %>% strsplit("/"), .by = ID) %>%
                    mutate(
                        Description = Description %>% as.character(),
                        col = "grey"
                    ) %>%
                    dplyr::select(db, Description, interpretation, Symb, NES) %>%
                    unnest(Symb) %>%
                    relocate(Symb) %>%
                    arrange(Symb)
            }
        )
    ) %>%
    dplyr::select(-genes_gsea)



# WRANGLE THE DATA FOR VOLCANO PLOTS ####
data_volcano <- limma_tables %>%
    dplyr::filter(design == "Baseline_vs_Week52") %>%
    dplyr::select(design, table) %>%
    mutate(
        data_plot = map(
            table,
            function(table) {
                table %>%
                    mutate(p = -log10(p)) %>%
                    suppressWarnings() %>%
                    dplyr::select(AffyID, Symb, logFC, p)
            }
        )
    )
data_volcano$data_plot


# WRANGLE THE DATA FOR ENRICHMENT DIAGRAMS ####
data_enrichment <- gsea %>%
    left_join(data_volcano, by = "design") %>%
    mutate(
        data_enrichment = pmap(
            list(gsea_table, data_plot),
            function(gsea_table, data_plot) {
                annotation <- gsea_table %>%
                    dplyr::slice_max(abs(NES), by = c("interpretation", "Symb"), with_ties = FALSE) %>%
                    dplyr::select(Symb, interpretation)
                data_plot %>%
                    dplyr::left_join(annotation, by = "Symb") %>%
                    tidyr::drop_na(interpretation)
            }
        )
    ) %>%
    dplyr::select(design, data_enrichment) %>%
    unnest(everything()) %>%
    nest(.by = c("design", "interpretation"))




# UNIVERAL PARAMETERS FOR VOLCANO PLOTS ####
ylim <- data_volcano %>%
    dplyr::select(data_plot) %>%
    unnest(data_plot) %>%
    pull(logFC) %>%
    range()
xlim <- c(0, data_volcano %>%
    dplyr::select(data_plot) %>%
    unnest(data_plot) %>%
    pull(p) %>%
    max() * 1.3)
x_end <- 3.25


# MAKE VOLCANO PLOT ####
plot_volcano


plot_volcano <- data_volcano %>%
    mutate(
        plot_volcano = pmap(
            list(data_plot, design),
            gg_volcano,
            point_size_null = 3,
            point_size = 5,
            ylim = ylim, 
            xlim = xlim
        )
    )
plot_volcano$plot_volcano[[1]]



# DEFINE CURVE GRADIENT LAYER ####
plot_volcano_curve_gradient <- ggplot2::geom_segment(
    inherit.aes = FALSE,
    data = dplyr::tibble(
        x0 = seq(
            -log10(0.05),
            x_end,
            length.out = 400
        ),
        xend = x0,
        y0 = -Inf,
        yend = Inf,
        col = x0
    ),
    mapping = ggplot2::aes(
        x = x0,
        xend = xend,
        y = y0,
        yend = yend,
        col = x0
    ),
    show.legend = FALSE
)


# DEFINE CURVE LAYERS ####
curves_baseline_vs_week24_immune_response <- ggplot2::geom_curve(
    data = data_enrichment %>%
        dplyr::filter(
            design == "Baseline_vs_Week24",
            interpretation == "immune response"
        ) %>%
        pull(data) %>%
        pluck(1),
    mapping = ggplot2::aes(
        x = p,
        y = logFC,
        xend = 3.1,
        yend = 0.15,
    ),
    col = "#5d00ff",
    curvature = 0.325,
    angle = 45,
    linewidth = 0.1,
    ncp = 1
)

curves_baseline_vs_week52_immune_response <- ggplot2::geom_curve(
    data = data_enrichment %>%
        dplyr::filter(
            design == "Baseline_vs_Week52",
            interpretation == "immune-related response to injury"
        ) %>%
        pull(data) %>%
        pluck(1),
    mapping = ggplot2::aes(
        x = p,
        y = logFC,
        xend = 2.95,
        yend = -0.75,
    ),
    col = "#d300ffff",
    curvature = 0.65,
    angle = 45,
    linewidth = 0.1,
    ncp = 1
)

curves_baseline_vs_week52_injury_response <- ggplot2::geom_curve(
    data = data_enrichment %>%
        dplyr::filter(
            design == "Baseline_vs_Week52",
            interpretation == "response to injury"
        ) %>%
        pull(data) %>%
        pluck(1),
    mapping = ggplot2::aes(
        x = p,
        y = logFC,
        xend = 3,
        yend = 0.45,
    ),
    col = "#ff9900",
    curvature = -0.15,
    angle = 90,
    linewidth = 0.1,
    ncp = 1
)



# DEFINE VOLCANO PLOTS ####
plot_volcano_baseline_vs_week52 <- plot_volcano %>%
    dplyr::filter(design == "Baseline_vs_Week52") %>%
    pull(plot_volcano) %>%
    pluck(1)


# RE-ARRANGE LAYERS OF GGPLOT ####
plot_volcano_baseline_vs_week52$layers <- c(
    curves_baseline_vs_week52_immune_response,
    curves_baseline_vs_week52_injury_response,
    plot_volcano_curve_gradient,
    plot_volcano_baseline_vs_week52 %>% pluck("layers", 1),
    plot_volcano_baseline_vs_week52 %>% pluck("layers", 2),
    plot_volcano_baseline_vs_week52 %>% pluck("layers", 3),
    plot_volcano_baseline_vs_week52 %>% pluck("layers", 4),
    plot_volcano_baseline_vs_week52 %>% pluck("layers", 5),
    plot_volcano_baseline_vs_week52 %>% pluck("layers", 6)
)


# FINALIZE PLOTS ####
plot_volcano_baseline_vs_week52 <- plot_volcano_baseline_vs_week52 +
    ggplot2::scale_colour_gradient(low = "#ffffff95", high = "#ffffff00")



# MAKE PLOT PANELS ####
plot_volcano_enrichment <- tibble(
    plot_volcano_enrichment = list(
        plot_volcano_baseline_vs_week52 %>%
            ggpubr::ggarrange(., NULL, ncol = 2, widths = c(1, 0))
    )
)



# SAVE THE PLOT DATA ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/"
save(plot_volcano_enrichment, file = paste(saveDir, "Volcano enrichment plots baseline to week 52 IQR_filtered_probes_unique_genes_baseline_corrected_cortex_corrected 1208.RData", sep = ""))



# SAVE THE PLOTS ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/output/"
# ggsave(
#     filename = paste(saveDir, "volcano enrichment draft.png"),
#     plot = plot_volcano_enrichment$plot_volcano_enrichment[[3]],
#     dpi = 300,
#     width = 20,
#     height = 20,
#     units = "cm",
#     bg = "white"
# )