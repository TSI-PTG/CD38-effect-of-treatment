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
source("C:/R/CD38-effect-of-treatment/code/functions/plot.gg_volcano_enrichment.r")
# load reference data
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/IQR_filtered_probes_unique_genes_baseline_corrected_cortex_corrected_limma_1208.RData")
# load enrichment results
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/felzartamab_gsea_baseline_corrected_cortex_corrected_k1208.RData")
# load simplified enrichment plots
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/simplified_enrichment_plots.RData")
# load simplified GO annotations
source("C:/R/CD38-effect-of-treatment/code/data management/functional enrichment annotation/simplified GO annotation.r")



# JOIN THE DE AND ENRICHMENT DATA ####
felzartamab_gsea_k1208$gsea_go_tables
data_joined_00 <- limma_tables %>%
    left_join(felzartamab_gsea_k1208 %>% dplyr::select(-toptable, -table, -gsea_go_flextables), by = "design") %>%
    dplyr::select(design, table, genes_gsea, gsea_go_tables)


# WRANGLE THE DATA FOR PLOTTING ####
data_joined_01 <- data_joined_00 %>%
    dplyr::select(design, table, genes_gsea, gsea_go_tables) %>%
    # expand_grid(direction = c("increased", "decreased")) %>%
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
        ),
        data_enrichment = map(
            gsea_go_tables,
            function(gsea_go_tables) {
                gsea_go_tables %>%
                    slice_min(pvalue, n = 20) %>%
                    mutate(
                        ID = ID %>% str_remove(":"),
                        GO = Description %>% as.character(),
                        col = "grey"
                    ) %>%
                    mutate(Symb = core_enrichment %>% strsplit("/"), .by = ID) %>%
                    dplyr::select(GO, Symb, NES) %>%
                    unnest(Symb) %>%
                    relocate(Symb) %>%
                    arrange(Symb) %>%
                    mutate(
                        group = case_when(
                            GO %>% str_detect(immune_response) ~ "immune response",
                            GO %>% str_detect(cell_cycle) ~ "cell cycling",
                            GO %>% str_detect(infection_response) ~ "response to infection",
                            GO %>% str_detect(inflammation) ~ "inflammation",
                            GO %>% str_detect(injury) ~ "injury response",
                            GO %>% str_detect(metabolic_response) ~ "metabolic response",
                            GO %>% str_detect(external_stimulus) ~ "response to external stimilus",
                            GO %>% str_detect(reg_cellular_processes) ~ "regulation of cellular processes",
                            GO %>% str_detect(cellular_development) ~ "cellular development",
                            GO %>% str_detect(cellular_communication) ~ "cellular communication",
                        ) %>% factor()
                    ) %>%
                    mutate(
                        count = n(),
                        .by = Symb,
                        .after = "Symb"
                    ) %>%
                    mutate(
                        prop = count / max(count),
                        .by = group,
                        .after = "count"
                    ) %>%
                    mutate(
                        n_group = group %>% unique() %>% length(),
                        groupID = group %>% as.numeric(),
                        group_mult = groupID / n_group
                    )
            }
        ),
        data_annotation = map2(
            data_plot, data_enrichment,
            function(data_plot, data_enrichment) {
                data_plot %>%
                    mutate(
                        p = -log10(`<U+0394><U+0394> p`),
                        logFC = `<U+0394><U+0394> logFC`
                    ) %>%
                    dplyr::select(AffyID, Symb, p, logFC) %>%
                    suppressWarnings() %>%
                    right_join(data_enrichment, by = "Symb")
            }
        )
    ) %>%
    dplyr::select(design, data_plot, data_enrichment, data_annotation)



# MAKE PLOTS ####
df_plot <- data_joined_01 %>%
    mutate(
        plot = pmap(
            list(data_plot, data_annotation, design),
            gg_volcano_enrichment
        )
    )





# ADD SIMPLIFIED PATHWAY ENRICHMENT PLOTS ####
plot_baseline_week24 <- df_plot %>%
    dplyr::filter(design == "Baseline_vs_Week24") %>%
    pull(plot) %>%
    pluck(1) +
    patchwork::inset_element(
        simplified_enrichment_plot %>%
            dplyr::filter(design == "Baseline_vs_Week24") %>%
            pull(plot_enrichment) %>%
            pluck(1) %>%
            pull(plot) %>%
            pluck(1),
        align_to = "plot",
        left = 0.78, right = 1,
        bottom = 0.275, top = 0.45
    ) +
    patchwork::inset_element(
        simplified_enrichment_plot %>%
            dplyr::filter(design == "Baseline_vs_Week24") %>%
            pull(plot_enrichment) %>%
            pluck(1) %>%
            pull(plot) %>%
            pluck(2),
        align_to = "plot",
        left = 0.78, right = 1,
        bottom = 0.075, top = 0.25
    )

plot_week24_week52 <- df_plot %>%
    dplyr::filter(design == "Week24_vs_Week52") %>%
    pull(plot) %>%
    pluck(1) +
    patchwork::inset_element(
        simplified_enrichment_plot %>%
            dplyr::filter(design == "Week24_vs_Week52") %>%
            pull(plot_enrichment) %>%
            pluck(1) %>%
            pull(plot) %>%
            pluck(1),
        align_to = "plot",
        left = 0.78, right = 1,
        bottom = 0.275, top = 0.45
    ) +
    patchwork::inset_element(
        simplified_enrichment_plot %>%
            dplyr::filter(design == "Week24_vs_Week52") %>%
            pull(plot_enrichment) %>%
            pluck(1) %>%
            pull(plot) %>%
            pluck(2),
        align_to = "plot",
        left = 0.78, right = 1,
        bottom = 0.075, top = 0.25
    )

plot_baseline_week52 <- df_plot %>%
    dplyr::filter(design == "Baseline_vs_Week52") %>%
    pull(plot) %>%
    pluck(1) +
    patchwork::inset_element(
        simplified_enrichment_plot %>%
            dplyr::filter(design == "Baseline_vs_Week52") %>%
            pull(plot_enrichment) %>%
            pluck(1) %>%
            pull(plot) %>%
            pluck(1),
        align_to = "plot",
        left = 0.68, right = 1,
        bottom = 0.425, top = 0.6
    ) +
    patchwork::inset_element(
        simplified_enrichment_plot %>%
            dplyr::filter(design == "Baseline_vs_Week52") %>%
            pull(plot_enrichment) %>%
            pluck(1) %>%
            pull(plot) %>%
            pluck(2),
        align_to = "plot",
        left = 0.78, right = 1,
        bottom = 0.325, top = 0.5
    )+
    patchwork::inset_element(
        simplified_enrichment_plot %>%
            dplyr::filter(design == "Baseline_vs_Week52") %>%
            pull(plot_enrichment) %>%
            pluck(1) %>%
            pull(plot) %>%
            pluck(3),
        align_to = "plot",
        left = 0.80, right = 1,
        bottom = 0.15, top = 0.325
    ) +
    patchwork::inset_element(
        simplified_enrichment_plot %>%
            dplyr::filter(design == "Baseline_vs_Week52") %>%
            pull(plot_enrichment) %>%
            pluck(1) %>%
            pull(plot) %>%
            pluck(4),
        align_to = "plot",
        left = 0.55, right = 1,
        bottom = 0.05, top = 0.225
    )




# MAKE PLOT PANELS ####
plot_volcano_enrichment <- df_plot %>%
    mutate(
        plot_volcano_enrichment = list(
            plot_baseline_week24,
            plot_week24_week52,
            plot_baseline_week52
        )
    )


# SAVE THE PLOT DATA ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/"
save(plot_volcano_enrichment, file = paste(saveDir, "Volcano enrichment plots IQR_filtered_probes_unique_genes_baseline_corrected_cortex_corrected 1208.RData", sep = ""))



# SAVE THE PLOTS ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/output/"
ggsave(
    filename = paste(saveDir, "volcano enrichment draft.png"),
    plot = plot_baseline_week52,
    dpi = 300,
    width = 20,
    height = 20,
    units = "cm",
    bg = "white"
)




# # DEFINE TEMPORARY DATA FOR DRAFTING PLOTS ####
# data_annotation <- data_joined_01$data_annotation[[2]]
# data <- data_joined_01$data_plot[[2]]
# df_enrichement <- simplified_enrichment_plot$plot_enrichment[[2]]


# # DEFINE LINES ####
# GO_lines <- data_annotation %>%
#     arrange(count %>% desc()) %>%
#     distinct(Symb, group, .keep_all = TRUE) %>%
#     mutate(
#         y_min = min(logFC),
#         y_max = max(logFC),
#         p2 = max(data$p) * 1.1,
#         logFC2 = ifelse(logFC < 0, y_min * group_mult, y_max * group_mult),
#         curvature = dplyr::case_when(
#             groupID == 1 ~ -0.25,
#             groupID == 2 ~ 0.25,
#             groupID == 3 ~ -0.25,
#             groupID == 4 ~ 0.25,
#             TRUE ~ 0.25
#         ),
#         conflict = dplyr::case_when(
#             logFC < 0 & NES > 0 ~ FALSE,
#             logFC > 0 & NES < 0 ~ FALSE,
#             TRUE ~ TRUE
#         )
#     ) %>%
#     dplyr::filter(conflict)





# GO_labels <- data_annotation %>%
#     distinct(Symb, group, .keep_all = TRUE) %>%
#     slice_max(prop, n = 5, by = c("group"), with_ties = FALSE) %>%
#     arrange(group, prop %>% desc()) %>%
#     left_join(GO_lines %>% distinct(Symb, group, .keep_all = TRUE))



# # PLOTTING GLOBALS ####
# ylim <- data$logFC %>% range() * 1.25
# xlim <- c(0, data$p %>% max() * 1.4)



# plot_temp <- data %>%
#     ggplot2::ggplot(mapping = ggplot2::aes(x = p, y = logFC)) +
#     Map(
#         function(i) {
#             geom_curve(
#                 data = GO_lines,
#                 mapping = ggplot2::aes(
#                     x = p,
#                     y = logFC,
#                     xend = p2,
#                     yend = logFC2,
#                     col = group,
#                     group = group,
#                 ),
#                 curvature = i,
#                 angle = 30,
#                 linewidth = 0.1
#             )
#         },
#         i = GO_lines$curvature
#     ) +
#     ggnewscale::new_scale_colour() +
#     ggplot2::geom_segment(
#         inherit.aes = FALSE,
#         data = tibble(
#             x0 = seq(
#                 -log10(0.05),
#                 GO_lines$p2 %>% max(),
#                 length.out = 750
#             ),
#             xend = x0,
#             y0 = -Inf,
#             yend = Inf,
#             col = x0
#         ),
#         mapping = ggplot2::aes(
#             x = x0,
#             xend = xend,
#             y = y0,
#             yend = yend,
#             col = x0
#         ),
#         show.legend = FALSE
#     ) +
#     ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
#     ggplot2::geom_vline(xintercept = -log10(0.05), linetype = "dashed") +
#     ggplot2::geom_point(col = ifelse(data$p < -log10(0.05), "grey70", "black")) +
#     ggplot2::geom_point(
#         show.legend = FALSE,
#         data = GO_labels %>% distinct(group, .keep_all = TRUE),
#         aes(x = p2, y = logFC2, fill = group), shape = 21, size = 4, col = "black"
#     ) +
#     scale_colour_gradient(low = "#ffffff95", high = "#ffffff00") +
#     coord_cartesian(
#         xlim = xlim,
#         ylim = ylim
#     ) +
#     theme_bw() +
#     theme(
#         legend.position = "none",
#         panel.grid = element_blank()
#     )
