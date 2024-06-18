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


# PATHWAY KEYS ####
immune_response <- paste(c("immune", "immunity", "cytokine", "leukocyte", "cell activation", "response to", "interaction", "virus", "symbiont", "defense response"), collapse = "|")
cell_cycle <- paste(c("cycle"), collapse = "|")
inflammation <- paste(c("inflam"), collapse = "|")
injury <- paste(c("injury"), collapse = "|")
external_stimulus <- paste(c("response to", "interaction"), collapse = "|")
cellular_regulation <- paste(c("regulation of"), collapse = "|")
cellular_development <- paste(c(
    "chromosome", "organelle fission", "organization", "segregation", "division",
    "development", "neurogenesis", "generation", "morphogenesis", "differentiation", "component"
), collapse = "|")
cell_signalling <- paste(c("communication", "signal", "signalling"), collapse = "|")
infection_response <- paste(c("virus", "symbiont", "defense response"), collapse = "|")
metabolic_response <- paste(c("metabolism", "metabolic", "catabolic"), collapse = "|")


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
        data_DE = map(
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
                            GO %>% str_detect(infection_response) ~ "response to infection",
                            GO %>% str_detect(external_stimulus) ~ "response to external stimilus",
                            GO %>% str_detect(inflammation) ~ "inflammation",
                            GO %>% str_detect(injury) ~ "injury response",
                            GO %>% str_detect(cell_cycle) ~ "cell cycling",
                            GO %>% str_detect(cell_signalling) ~ "cell signalling",
                            GO %>% str_detect(cell_mobilization) ~ "cellular mobilization",
                            GO %>% str_detect(cellular_development) ~ "cellular development",
                            GO %>% str_detect(cellular_regulation) ~ "cellular regulation",
                            GO %>% str_detect(protein_metabolism) ~ "protein metabolism",
                            GO %>% str_detect(nitrogen_metabolism) ~ "nitrogen metabolism",
                            GO %>% str_detect(xenobiotic_metabolism) ~ "xenobiotic metabolism",
                            GO %>% str_detect(metabolic_response) ~ "metabolic response",
                        ) %>% factor(levels = GO_annotation_levels),
                        col_group = case_when(
                            group == "immune response" ~ "#ffb700",
                            group == "response to infection" ~ "#ff0000",
                            group == "response to external stimilus" ~ "#00ff91",
                            group == "inflammation" ~ "#ff9900",
                            group == "injury response" ~ "#5d00ff",
                            group == "cell cycling" ~ "#00ff33",
                            group == "cell signalling" ~ "#00ff33",
                            group == "cellular mobilization" ~ "#00ff33",
                            group == "cellular development" ~ "#4dff00",
                            group == "cellular regulation" ~ "#ff00ea",
                            group == "protein metabolism" ~ "#7b00ff",
                            group == "nitrogen metabolism" ~ "#7b00ff",
                            group == "xenobiotic metabolism" ~ "#7b00ff",
                            group == "metabolic response" ~ "#7b00ff",
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
                        groupID = group %>% as.numeric(),
                        group_mult = groupID / n_group
                    )
            }
        ),
        data_joined = map2(
            data_DE, data_enrichment,
            function(data_DE, data_enrichment) {
                data_DE %>%
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
    dplyr::select(design, data_DE, data_enrichment, data_joined)

# data_joined_01$data_DE[[1]]
# data_joined_01$data_enrichment[[1]]
# data_joined_01$data_joined[[1]] %>%
#     arrange(count %>% desc()) %>%
#     print(n = "all")



# DEFINE LINES ####
df_lines <- data_joined_01$data_joined[[3]] %>%
    arrange(count %>% desc()) %>%
    distinct(Symb, group, .keep_all = TRUE) %>%
    mutate(,
        y_min = min(logFC),
        y_max = max(logFC),
        p2 = max(p) * 1.1,
        logFC2 = ifelse(logFC < 0, y_min * group_mult, y_max * group_mult),
        curvature = dplyr::case_when(
            groupID == 1 ~ -0.25,
            groupID == 2 ~ 0.25,
            groupID == 3 ~ -0.25,
            groupID == 4 ~ 0.25,
            TRUE ~ 0.25
        )
    )

df_labels <- data_joined_01$data_joined[[3]] %>%
    distinct(Symb, group, .keep_all = TRUE) %>%
    slice_max(prop, n = 5, by = c("group"), with_ties = FALSE) %>%
    arrange(group, prop %>% desc()) %>%
    left_join(df_lines %>% distinct(Symb, group, .keep_all = TRUE))



# PLOTTING GLOBALS ####
ylim <- data_joined_01$data_DE[[3]]$logFC  %>% range * 1.25
xlim <- c(0, data_joined_01$data_DE[[3]]$p  %>% max *1.25)


# MAKE PLOTS ####
data_joined_01$data_DE[[3]] %>%
    ggplot2::ggplot(mapping = ggplot2::aes(x = p, y = logFC)) +
    Map(
        function(i) {
            geom_curve(
                data = df_lines,
                mapping = ggplot2::aes(
                    x = p,
                    y = logFC,
                    xend = p2,
                    yend = logFC2,
                    col = group,
                    group = group,
                ),
                curvature = i,
                angle = 30,
                line_size = 0.125
            )
        },
        i = df_lines$curvature
    ) +
    ggnewscale::new_scale_colour() +
    ggplot2::geom_segment(
        inherit.aes = FALSE,
        data = tibble(
            x0 = seq(
                -log10(0.05),
                df_lines$p2 %>% max(),
                length.out = 1000
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
    ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::geom_vline(xintercept = -log10(0.05), linetype = "dashed") +
    ggplot2::geom_point() +
    ggplot2::geom_point(
        show.legend = FALSE,
        data = df_labels %>% distinct(group, .keep_all = TRUE),
        aes(x = p2, y = logFC2, fill = group), shape = 21, size = 4, col = "black"
    ) +
    ggrepel::geom_label_repel(
        data = df_labels,
        aes(x = p2, y = logFC2, fill = group, label = Symb),
        size = 0.5,
        nudge_x = 0.05,
        hjust = "left",
        direction = "y",
        min.segment.length = 10,
        seed = 42,
        label.padding = 0.1,
        max.overlaps = Inf,
        show.legend = FALSE
    ) +
    scale_colour_gradient(low = "#ffffff95", high = "#ffffff00") +
    coord_cartesian(
        xlim = xlim,
        ylim = ylim
    ) +
    theme_bw() +
    theme(
        legend.position = "none",
        panel.grid = element_blank()
    )









# # MAKE SINGLE VOLCANO PLOTS ####
# plot_volcano <- data_plot %>%
#     mutate(
#         plot_volcano_enrichment = map2(
#             data_plot, design,
#             gg_volcano,
#             x_break = 1.75,
#             labels_probes = probes_abmr,
#             point_size_null = 1.25,
#             point_size = 2.5,
#             labels_probes_size = 3,
#         )
#     )



# plot_volcano$plot_volcano[[3]]


