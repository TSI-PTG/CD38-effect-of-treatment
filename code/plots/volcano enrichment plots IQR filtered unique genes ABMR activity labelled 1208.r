# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(ggrepel) # install.packages("ggrepel")
library(ggpubr) # install.packages("ggpubr")
library(patchwork) # install.packages("patchwork")
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
cell_cycle <- paste(c("cycle", "division"), collapse = "|")
inflammation <- paste(c("inflam"), collapse = "|")
injury <- paste(c("injury"), collapse = "|")
external_stimulus <- paste(c("response to", "interaction"), collapse = "|")
reg_cellular_processes <- paste(c("regulation of"), collapse = "|")
cellular_development <- paste(c(
    "chromosome", "organelle fission", "organization", "segregation",
    "development", "neurogenesis", "generation", "morphogenesis", "differentiation", "component"
), collapse = "|")
cellular_communication <- paste(c("communication", "signal", "signalling"), collapse = "|")
infection_response <- paste(c("virus", "symbiont", "defense response"), collapse = "|")
metabolic_response <- paste(c("metabolism", "metabolic", "catabolic"), collapse = "|")


# JOIN THE DE AND ENRICHMENT DATA ####
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
                    dplyr::select(GO, Symb) %>%
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
                        )
                    ) %>%
                    mutate(count = n(), .by = Symb, .after = "Symb")
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

data_joined_01$data_DE[[1]]
data_joined_01$data_enrichment[[1]]
data_joined_01$data_joined[[1]] %>%
    arrange(count %>% desc()) %>%
    print(n = "all")

labels_tmp <- data_joined_01$data_joined[[1]] %>%
    arrange(count %>% desc()) %>%
    distinct(Symb, group, .keep_all = TRUE) %>%
    mutate(
        p2 = 4, 
        logFC2 = ifelse(group  %>% str_detect("immune"), -0.75, -0.25)
    )


data_joined_01$data_DE[[1]] %>%
    ggplot2::ggplot(mapping = ggplot2::aes(x = p, y = logFC)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::geom_vline(xintercept = -log10(0.05), linetype = "dashed") +
    ggplot2::geom_point() +
    geom_curve(
        data = labels_tmp,
        mapping = ggplot2::aes(
            x = p,
            y = logFC,
            xend = p2,
            yend = logFC2,
            col = p, 
            group = group
        ),
        curvature = 0.25
    )




# MAKE SINGLE VOLCANO PLOTS ####
plot_volcano <- data_plot %>%
    mutate(
        plot_volcano_enrichment = map2(
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
