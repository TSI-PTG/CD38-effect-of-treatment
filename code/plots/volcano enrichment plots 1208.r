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
source("C:/R/CD38-effect-of-treatment/code/data management/functional enrichment annotation/very simplified GO annotation.r")



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
                            GO %>% str_detect(infection_response) ~ "response to infection",
                            GO %>% str_detect(response_to_stimulus) ~ "response to exogenous/endogenous stimulus",
                            GO %>% str_detect(inflammation) ~ "inflammation",
                            GO %>% str_detect(injury) ~ "injury response",
                            GO %>% str_detect(cell_signalling_and_RNA_transcription) ~ "cell signalling and RNA transcription",
                            GO %>% str_detect(cellular_development_and_metabolism) ~ "cell development,\nmobilization,\nand metabolism",
                            TRUE ~ GO
                        ) %>% factor(levels = GO_annotation_levels_truncated),
                        col_group = case_when(
                            group == "immune response" ~ "#5d00ff",
                            group == "response to infection" ~ "#ff0000",
                            group == "response to exogenous/endogenous stimulus" ~ "#ff9900",
                            group == "inflammation" ~ "#ff9900",
                            group == "injury response" ~ "#5d00ff",
                            group == "cell signalling and RNA transcription" ~ "#ff00ee",
                            group == "cell development,\nmobilization,\nand metabolism" ~ "#4dff00"
                        )
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
                        groupID = group %>% factor() %>% as.numeric(),
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


# MAKE ONE-OFF EDITS TO GO ANNOTATIONS ####
data_joined_01 <- data_joined_01 %>%
    mutate(
        data_enrichment = pmap(
            list(design, data_enrichment),
            function(design, data_enrichment) {
                data_enrichment %>%
                    mutate(
                        group = case_when(
                            design == "Baseline_vs_Week24" ~ "immune response",
                            TRUE ~ group
                        ) %>% factor(levels = GO_annotation_levels_truncated),
                        # n_group = group %>% unique() %>% length(),
                        # groupID = group %>% as.numeric(),
                        # group_mult = groupID / n_group
                    )
            }
        ),
        data_annotation = pmap(
            list(design, data_annotation),
            function(design, data_annotation) {
                data_annotation %>%
                    mutate(
                        group = case_when(
                            design == "Baseline_vs_Week24" ~ "immune response",
                            TRUE ~ group
                        ) %>% factor(levels = GO_annotation_levels_truncated),
                        # n_group = group %>% unique() %>% length(),
                        # groupID = group %>% as.numeric(),
                        # group_mult = groupID / n_group
                    )
            }
        )
    )

data_joined_01$data_annotation[[3]] %>%
    dplyr::select(GO, group, n_group, groupID, group_mult)






# DEFINE THE PATHWAY ENRICHMENT LINES ####
data_joined_02 <- data_joined_01 %>%
    mutate(
        GO_lines = pmap(
            list(data_plot, data_annotation),
            function(data_plot, data_annotation) {
                xlim <- c(0, data_plot$p %>% max() * 1.4)
                y_min <- data_plot %>%
                    dplyr::pull(logFC) %>%
                    min()
                y_max <- data_plot %>%
                    dplyr::pull(logFC) %>%
                    max()
                y_diff <- c(y_min, y_max) %>% diff()
                x_end <- data_plot %>%
                    dplyr::pull(p) %>%
                    max() * 1.1
                n_group <- data_annotation %>%
                    pull(n_group) %>%
                    unique()
                interval <- y_diff / (n_group + 1)
                interval_prop <- 1 / (n_group + 1)
                # interval <- y_diff / (n_group)
                # interval_prop <- 1 / (n_group)
                data_annotation %>%
                    dplyr::arrange(count %>% dplyr::desc()) %>%
                    dplyr::distinct(Symb, group, .keep_all = TRUE) %>%
                    dplyr::mutate(
                        y_min = y_min,
                        y_max = y_max,
                        y_diff = y_diff,
                        interval = interval,
                        interval_prop = interval_prop,
                        x_end = x_end,
                        y_end = y_min + (groupID * interval),
                        y_end_prop = groupID * interval_prop,
                        x_end_prop = x_end / xlim[2],
                        curvature = dplyr::case_when(
                            groupID == 1 ~ -0.25,
                            groupID == 2 ~ 0.25,
                            groupID == 3 ~ -0.25,
                            groupID == 4 ~ 0.25,
                            TRUE ~ 0.25
                        ),
                        conflict = dplyr::case_when(
                            logFC < 0 & NES > 0 ~ FALSE,
                            logFC > 0 & NES < 0 ~ FALSE,
                            TRUE ~ TRUE
                        )
                    ) %>%
                    dplyr::filter(conflict) %>%
                    dplyr::select(-AffyID)
            }
        )
    )
data_joined_02$GO_lines[[3]]
# data_joined_02$GO_lines[[2]] %>%
#     dplyr::select(GO, group)  %>% print(n="all")


# MAKE PLOTS ####
df_plot <- data_joined_02 %>%
    mutate(
        plot = pmap(
            list(data_plot, data_annotation, GO_lines, design),
            gg_volcano_enrichment
        )
    )
# df_plot$plot[[1]]



# ADD SIMPLIFIED ENRICHMENT ANNOTATION TO Baseline_vs_Week24 PLOTS ####
plot_baseline_week24 <- df_plot %>%
    dplyr::filter(design == "Baseline_vs_Week24") %>%
    pull(plot) %>%
    pluck(1)
    #  +
    # patchwork::inset_element(
    #     simplified_enrichment_plot %>%
    #         dplyr::filter(design == "Baseline_vs_Week24") %>%
    #         pull(plot_enrichment) %>%
    #         pluck(1) %>%
    #         pull(plot) %>%
    #         pluck(1),
    #     align_to = "plot",
    #     left = df_plot %>%
    #         dplyr::filter(design == "Baseline_vs_Week24") %>%
    #         pull(GO_lines) %>%
    #         pluck(1) %>%
    #         dplyr::filter(groupID == 1) %>%
    #         pull(x_end_prop) %>%
    #         unique(),
    #     right = 0.995,
    #     bottom = df_plot %>%
    #         dplyr::filter(design == "Baseline_vs_Week24") %>%
    #         pull(GO_lines) %>%
    #         pluck(1) %>%
    #         dplyr::filter(groupID == 1) %>%
    #         pull(y_end_prop) %>%
    #         unique() - 0.1,
    #     top = df_plot %>%
    #         dplyr::filter(design == "Baseline_vs_Week24") %>%
    #         pull(GO_lines) %>%
    #         pluck(1) %>%
    #         dplyr::filter(groupID == 1) %>%
    #         pull(y_end_prop) %>%
    #         unique() + 0.1
    # )



# ADD SIMPLIFIED ENRICHMENT ANNOTATION TO Week24_vs_Week52 PLOTS ####
plot_week24_week52 <- df_plot %>%
    dplyr::filter(design == "Week24_vs_Week52") %>%
    pull(plot) %>%
    pluck(1) 
    # +
    # patchwork::inset_element(
    #     simplified_enrichment_plot %>%
    #         dplyr::filter(design == "Week24_vs_Week52") %>%
    #         pull(plot_enrichment) %>%
    #         pluck(1) %>%
    #         pull(plot) %>%
    #         pluck(1),
    #     align_to = "plot",
    #     left = df_plot %>%
    #         dplyr::filter(design == "Week24_vs_Week52") %>%
    #         pull(GO_lines) %>%
    #         pluck(1) %>%
    #         dplyr::filter(groupID == 1) %>%
    #         pull(x_end_prop) %>%
    #         unique(),
    #     right = 0.995,
    #     bottom = df_plot %>%
    #         dplyr::filter(design == "Week24_vs_Week52") %>%
    #         pull(GO_lines) %>%
    #         pluck(1) %>%
    #         dplyr::filter(groupID == 1) %>%
    #         pull(y_end_prop) %>%
    #         unique() - 0.1,
    #     top = df_plot %>%
    #         dplyr::filter(design == "Week24_vs_Week52") %>%
    #         pull(GO_lines) %>%
    #         pluck(1) %>%
    #         dplyr::filter(groupID == 1) %>%
    #         pull(y_end_prop) %>%
    #         unique() + 0.1
    # )



# ADD SIMPLIFIED ENRICHMENT ANNOTATION TO Baseline_vs_Week52 PLOTS ####
plot_baseline_week52 <- df_plot %>%
    dplyr::filter(design == "Baseline_vs_Week52") %>%
    pull(plot) %>%
    pluck(1) 
    # +
    # patchwork::inset_element(
    #     simplified_enrichment_plot %>%
    #         dplyr::filter(design == "Baseline_vs_Week52") %>%
    #         pull(plot_enrichment) %>%
    #         pluck(1) %>%
    #         pull(plot) %>%
    #         pluck(1),
    #     align_to = "plot",
    #     left = df_plot %>%
    #         dplyr::filter(design == "Baseline_vs_Week52") %>%
    #         pull(GO_lines) %>%
    #         pluck(1) %>%
    #         dplyr::filter(groupID == 1) %>%
    #         pull(x_end_prop) %>%
    #         unique(),
    #     right = 0.995,
    #     bottom = df_plot %>%
    #         dplyr::filter(design == "Baseline_vs_Week52") %>%
    #         pull(GO_lines) %>%
    #         pluck(1) %>%
    #         dplyr::filter(groupID == 1) %>%
    #         pull(y_end_prop) %>%
    #         unique() - 0.1,
    #     top = df_plot %>%
    #         dplyr::filter(design == "Baseline_vs_Week52") %>%
    #         pull(GO_lines) %>%
    #         pluck(1) %>%
    #         dplyr::filter(groupID == 1) %>%
    #         pull(y_end_prop) %>%
    #         unique() + 0.1
    # ) +
    # patchwork::inset_element(
    #     simplified_enrichment_plot %>%
    #         dplyr::filter(design == "Baseline_vs_Week52") %>%
    #         pull(plot_enrichment) %>%
    #         pluck(1) %>%
    #         pull(plot) %>%
    #         pluck(2),
    #     align_to = "plot",
    #     left = df_plot %>%
    #         dplyr::filter(design == "Baseline_vs_Week52") %>%
    #         pull(GO_lines) %>%
    #         pluck(1) %>%
    #         dplyr::filter(groupID == 2) %>%
    #         pull(x_end_prop) %>%
    #         unique(),
    #     right = 0.995,
    #     bottom = df_plot %>%
    #         dplyr::filter(design == "Baseline_vs_Week52") %>%
    #         pull(GO_lines) %>%
    #         pluck(1) %>%
    #         dplyr::filter(groupID == 2) %>%
    #         pull(y_end_prop) %>%
    #         unique() - 0.1,
    #     top = df_plot %>%
    #         dplyr::filter(design == "Baseline_vs_Week52") %>%
    #         pull(GO_lines) %>%
    #         pluck(1) %>%
    #         dplyr::filter(groupID == 2) %>%
    #         pull(y_end_prop) %>%
    #         unique() + 0.1
    # ) +
    # patchwork::inset_element(
    #     simplified_enrichment_plot %>%
    #         dplyr::filter(design == "Baseline_vs_Week52") %>%
    #         pull(plot_enrichment) %>%
    #         pluck(1) %>%
    #         pull(plot) %>%
    #         pluck(3),
    #     align_to = "plot",
    #     left = df_plot %>%
    #         dplyr::filter(design == "Baseline_vs_Week52") %>%
    #         pull(GO_lines) %>%
    #         pluck(1) %>%
    #         dplyr::filter(groupID == 3) %>%
    #         pull(x_end_prop) %>%
    #         unique(),
    #     right = 0.995,
    #     bottom = df_plot %>%
    #         dplyr::filter(design == "Baseline_vs_Week52") %>%
    #         pull(GO_lines) %>%
    #         pluck(1) %>%
    #         dplyr::filter(groupID == 3) %>%
    #         pull(y_end_prop) %>%
    #         unique() - 0.1,
    #     top = df_plot %>%
    #         dplyr::filter(design == "Baseline_vs_Week52") %>%
    #         pull(GO_lines) %>%
    #         pluck(1) %>%
    #         dplyr::filter(groupID == 3) %>%
    #         pull(y_end_prop) %>%
    #         unique() + 0.1
    # ) +
    # patchwork::inset_element(
    #     simplified_enrichment_plot %>%
    #         dplyr::filter(design == "Baseline_vs_Week52") %>%
    #         pull(plot_enrichment) %>%
    #         pluck(1) %>%
    #         pull(plot) %>%
    #         pluck(4),
    #     align_to = "plot",
    #     left = df_plot %>%
    #         dplyr::filter(design == "Baseline_vs_Week52") %>%
    #         pull(GO_lines) %>%
    #         pluck(1) %>%
    #         dplyr::filter(groupID == 4) %>%
    #         pull(x_end_prop) %>%
    #         unique(),
    #     right = 0.995,
    #     bottom = df_plot %>%
    #         dplyr::filter(design == "Baseline_vs_Week52") %>%
    #         pull(GO_lines) %>%
    #         pluck(1) %>%
    #         dplyr::filter(groupID == 4) %>%
    #         pull(y_end_prop) %>%
    #         unique() - 0.1,
    #     top = df_plot %>%
    #         dplyr::filter(design == "Baseline_vs_Week52") %>%
    #         pull(GO_lines) %>%
    #         pluck(1) %>%
    #         dplyr::filter(groupID == 4) %>%
    #         pull(y_end_prop) %>%
    #         unique() + 0.1
    # )


# MAKE PLOT PANELS ####
plot_volcano_enrichment <- df_plot %>%
    mutate(
        plot_volcano_enrichment = list(
            plot_baseline_week24 %>%
                ggpubr::ggarrange(., NULL, ncol = 2, widths = c(1, 0)),
            plot_week24_week52 %>%
                ggpubr::ggarrange(., NULL, ncol = 2, widths = c(1, 0)),
            plot_baseline_week52 %>%
                ggpubr::ggarrange(., NULL, ncol = 2, widths = c(1, 0))
        )
    )


# SAVE THE PLOT DATA ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/"
save(plot_volcano_enrichment, file = paste(saveDir, "Volcano enrichment plots IQR_filtered_probes_unique_genes_baseline_corrected_cortex_corrected 1208.RData", sep = ""))



# SAVE THE PLOTS ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/output/"
ggsave(
    filename = paste(saveDir, "volcano enrichment draft.png"),
    plot = plot_volcano_enrichment$plot_volcano_enrichment[[3]],
    dpi = 300,
    width = 20,
    height = 20,
    units = "cm",
    bg = "white"
)
