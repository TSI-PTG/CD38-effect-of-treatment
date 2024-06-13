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
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/IQR_filtered_probes_unique_genes_limma_1208.RData")
# load enrichment results
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/felzartamab_gsea_k1208.RData")



# PATHWAY KEYS ####
immune_response <- paste(
    c(
        "immune", "immunity", "cytokine", "leukocyte", "cell activation",
        "response to", "interaction", "virus", "symbiont", "defense response"
    ),
    collapse = "|"
)
cell_cycle <- paste(c("cycle"), collapse = "|")
inflammation <- paste(c("inflam"), collapse = "|")
injury <- paste(c("injury"), collapse = "|")
external_stimulus <- paste(c("response to", "interaction"), collapse = "|")
reg_cellular_processes <- paste(c("regulation of"), collapse = "|")
cellular_development <- paste(c(
    "chromosome", "organelle fission", "organization", "segregation", "division",
    "development", "neurogenesis", "generation", "morphogenesis", "differentiation", "component"
), collapse = "|")
cellular_communication <- paste(c("communication", "signal", "signalling"), collapse = "|")
infection_response <- paste(c("virus", "symbiont", "defense response"), collapse = "|")
metabolic_response <- paste(c("metabolism", "metabolic", "catabolic"), collapse = "|")


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



# DEFINE TEMPORARY DATA FOR DRAFTING PLOTS ####
df_plot_enrichment <- data_enrichment_plot$data_plot[[1]]




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
                                # data  %>% print
                                data %>%
                                    ggplot2::ggplot(mapping = ggplot2::aes(x = 1, y = Symb, label = Symb)) +
                                    ggplot2::geom_label(
                                        hjust = "left",
                                        size = 3
                                    ) +
                                    labs(
                                        title = paste(group, "\n(", data$n_in_group, " GO term(s) enriched)", sep = "")
                                    ) +
                                    ggplot2::theme_void() +
                                    ggplot2::coord_cartesian(xlim = c(1, 1.01), expand = 0.1) +
                                    ggplot2::theme(
                                        plot.title = element_text(hjust = 0, size = 10),
                                        plot.margin = unit(c(0, 0, 0, 0), "cm")
                                    )
                            }
                        )
                    )
            }
        )
    )



# SAVE THE PLOT DATA ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/"
save(simplified_enrichment_plot, file = paste(saveDir, "simplified_enrichment_plots.RData", sep = ""))



