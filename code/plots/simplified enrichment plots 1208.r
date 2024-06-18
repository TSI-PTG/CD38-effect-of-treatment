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
source("C:/R/CD38-effect-of-treatment/code/functions/plot.gg_volcano_timeseries.r")
# load reference data
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/IQR_filtered_probes_unique_genes_baseline_corrected_cortex_corrected_limma_1208.RData")
# load enrichment results
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/felzartamab_gsea_baseline_corrected_cortex_corrected_k1208.RData")
# load simplified GO annotations
source("C:/R/CD38-effect-of-treatment/code/data management/functional enrichment annotation/simplified GO annotation.r")


# JOIN THE DE AND ENRICHMENT DATA ####
felzartamab_gsea_k1208$gsea_go_tables
data_joined_00 <- limma_tables %>%
    left_join(felzartamab_gsea_k1208 %>% dplyr::select(-toptable, -table, -gsea_go_flextables), by = "design") %>%
    dplyr::select(design, table, genes_gsea, gsea_go_tables)


# FURTHER ANNOTATE THE ENRICHMENT RESULTS ####
data_enrichment <- data_joined_00 %>%
    dplyr::select(design, table, genes_gsea, gsea_go_tables) %>%
    mutate(
        data_enrichment = map(
            gsea_go_tables,
            function(gsea_go_tables) {
                gsea_go_tables %>%
                    # slice_min(pvalue, n = 20) %>%
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
                            GO %>% str_detect(external_stimulus) ~ "response to external stimilus",
                            GO %>% str_detect(inflammation) ~ "inflammation",
                            GO %>% str_detect(injury) ~ "injury response",
                            GO %>% str_detect(cell_cycle) ~ "cell cycling",
                            GO %>% str_detect(cell_signalling) ~ "cell signalling",
                            GO %>% str_detect(cell_mobilization) ~ "cell mobilization",
                            GO %>% str_detect(cellular_development) ~ "cell development",
                            GO %>% str_detect(cellular_regulation) ~ "cellular regulation",
                            # GO %>% str_detect(protein_metabolism) ~ "protein metabolism",
                            # GO %>% str_detect(nitrogen_metabolism) ~ "nitrogen metabolism",
                            # GO %>% str_detect(xenobiotic_metabolism) ~ "xenobiotic metabolism",
                            GO %>% str_detect(general_metabolic_response) ~ "metabolic response"
                            # ) %>% factor(levels = GO_annotation_levels),
                        ) %>% factor(levels = GO_annotation_levels_truncated),
                        col_group = case_when(
                            group == "immune response" ~ "#ffb700",
                            group == "response to infection" ~ "#ff0000",
                            group == "response to external stimilus" ~ "#00ff91",
                            group == "inflammation" ~ "#ff9900",
                            group == "injury response" ~ "#5d00ff",
                            group == "cell cycling" ~ "#00ff33",
                            group == "cell signalling" ~ "#00ff33",
                            group == "cell mobilization" ~ "#00ff33",
                            group == "cell development" ~ "#4dff00",
                            group == "cellular regulation" ~ "#ff00ea",
                            # group == "protein metabolism" ~ "#7b00ff",
                            # group == "nitrogen metabolism" ~ "#7b00ff",
                            # group == "xenobiotic metabolism" ~ "#7b00ff",
                            group == "metabolic response" ~ "#7b00ff"
                        )
                    ) %>%
                    mutate(
                        count = n(),
                        .by = c(group, Symb),
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
                    ) %>%
                    mutate(n_in_group = GO %>% unique() %>% length(), .by = group, .after = n_group)
            }
        )
    ) %>%
    dplyr::select(design, data_enrichment)


# WRANGLE DATA FOR PLOTTING ####
data_enrichment_plot <- data_enrichment %>%
    mutate(
        data_plot = map(
            data_enrichment,
            function(data_enrichment) {
                data_enrichment %>%
                    dplyr::distinct(Symb, group, .keep_all = TRUE) %>%
                    dplyr::slice_max(prop, n = 5, by = c("group"), with_ties = FALSE) %>%
                    dplyr::arrange(group, prop %>% dplyr::desc())
            }
        )
    )
# data_enrichment_plot$data_plot[[3]] %>% print(n = "all")


# MAKE SIMPLIFIED ENRICHEMENT PLOTS ####
simplified_enrichment_plot <- data_enrichment_plot %>%
    mutate(
        plot_enrichment = map(
            data_plot,
            function(data_plot) {
                data_plot %>%
                    nest(.by = group) %>%
                    mutate(
                        plot = pmap(
                            list(group, data),
                            function(group, data) {
                                data %>%
                                    ggplot2::ggplot(mapping = ggplot2::aes(x = 1, y = Symb, label = Symb)) +
                                    ggplot2::geom_label(
                                        hjust = "left",
                                        size = 2
                                    ) +
                                    labs(
                                        title = paste(" ", group, "\n (", data$n_in_group, " GO term(s) enriched)", sep = "")
                                    ) +
                                    ggplot2::theme_void() +
                                    ggplot2::coord_cartesian(xlim = c(1, 1.01), expand = 0.1) +
                                    ggplot2::theme(
                                        plot.title = element_text(hjust = 0, size = 6.5),
                                        plot.margin = unit(c(0, 0, 0, 0), "cm"),
                                        plot.background = element_rect(
                                            colour = data$col_group[[1]]
                                        )
                                    )
                            }
                        )
                    )
            }
        )
    )
# simplified_enrichment_plot$plot_enrichment[[1]]$plot


# SAVE THE PLOT DATA ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/"
save(simplified_enrichment_plot, file = paste(saveDir, "simplified_enrichment_plots.RData", sep = ""))
