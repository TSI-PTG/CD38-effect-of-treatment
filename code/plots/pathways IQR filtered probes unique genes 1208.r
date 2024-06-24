# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(clusterProfiler) # pak::pak("YuLab-SMU/clusterProfiler")
library(ggpubr) # install.packages("ggpubr")
library(patchwork) # install.packages("patchwork")
# Custom operators, functions, and datasets
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
# load enrichment results
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/felzartamab_gsea_k1208.RData")


# PATHWAY KEYS ####
immune_response <- paste(c("immune", "immunity", "cytokine", "leukocyte", "cell activation", "response to", "interaction", "virus", "symbiont", "defense response"), collapse = "|")
cell_cycle <- paste(c("cycle", "division"), collapse = "|")
inflammation <- paste(c("inflam"), collapse = "|")
injury <- paste(c("injury"), collapse = "|")
exogenous_stimulus <- paste(c("response to", "interaction"), collapse = "|")
cellular_regulation <- paste(c("regulation of"), collapse = "|")
cellular_development <- paste(c(
    "chromosome", "organelle fission", "organization", "segregation",
    "development", "neurogenesis", "generation", "morphogenesis", "differentiation", "component"
), collapse = "|")
cell_signalling <- paste(c("communication", "signal", "signalling"), collapse = "|")
infection_response <- paste(c("virus", "symbiont", "defense response"), collapse = "|")
metabolic_response <- paste(c("metabolism", "metabolic", "catabolic"), collapse = "|")



# WRANGLE DATA FOR PLOTTING ####
df_plot_00 <- felzartamab_gsea_k1208 %>%
    dplyr::select(design, genes_gsea, gsea_go_tables) %>%
    expand_grid(direction = c("increased", "decreased")) %>%
    mutate(
        data_plot = map2(
            gsea_go_tables,direction,
            function(gsea_go_tables,direction) {
                if (direction == "increased") {
                    dat <- gsea_go_tables %>%
                        dplyr::filter(NES > 0)
                } else if (direction == "decreased") {
                    dat <- gsea_go_tables %>%
                        dplyr::filter(NES < 0)
                }
                dat %>%
                    # slice_min(pvalue, n = 20) %>%
                    mutate(
                        ID = ID %>% str_remove(":"),
                        Description = Description %>% as.character(),
                        col = "grey"
                    ) %>%
                    mutate(Genes = core_enrichment %>% strsplit("/"), .by = ID) %>%
                    dplyr::select(Description, Genes) %>%
                    unnest(Genes) %>%
                    relocate(Genes) %>%
                    arrange(Genes) %>%
                    mutate(
                        group = case_when(
                            Description %>% str_detect(immune_response) ~ "immune response",
                            Description %>% str_detect(cell_cycle) ~ "cell cycling",
                            Description %>% str_detect(infection_response) ~ "response to infection",
                            Description %>% str_detect(inflammation) ~ "inflammation",
                            Description %>% str_detect(injury) ~ "injury response",
                            Description %>% str_detect(metabolic_response) ~ "metabolic response",
                            Description %>% str_detect(exogenous_stimulus) ~ "response to exogenous stimulus",
                            Description %>% str_detect(cellular_regulation) ~ "regulation of cellular processes",
                            Description %>% str_detect(cellular_development) ~ "cellular development",
                            Description %>% str_detect(cell_signalling) ~ "cell signalling",
                        )
                    )
            }
        )
    )
df_plot_00$data_plot[[3]]



# DEFINE PATHWAYS ####
df_plot_01 <- df_plot_00 %>%
    mutate(
        pathways = map(
            data_plot,
            function(data_plot) {
                data_plot %>%
                    distinct(Description) %>%
                    pull()
            }
        )
    ) %>%
    dplyr::filter(length(pathways) > 0)


# DEFINE GROUPS FOR PATHWAYS ####
df_plot_02 <- df_plot_01 %>%
    mutate(
        link_group = map(
            data_plot, function(data_plot) {
                data_plot %>%
                    pivot_longer(cols = c("Genes", "Description"), values_to = "Description") %>%
                    mutate(
                        group = case_when(
                            GO %>% str_detect(immune_response) ~ "immune response",
                            GO %>% str_detect(infection_response) ~ "response to infection",
                            GO %>% str_detect(exogenous_stimulus) ~ "response to exogenous stimulus",
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
                            group == "response to exogenous stimulus" ~ "#00ff91",
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
                    distinct(Description, .keep_all = TRUE) %>%
                    dplyr::pull(group, Description)
            }
        )
    )

df_plot_02$link_group



# MAKE CHORD PLOTS ####
df_plot_02 <- df_plot_01 %>%
    mutate(
        link_group = map(
            data_plot, function(data_plot) {
                data_plot %>%
                    pivot_longer(cols = c("Genes", "Description"), values_to = "Description") %>%
                    mutate(
                        group = case_when(
                            name == "Genes" ~ "Genes",
                            Description %>% str_detect(immune_response) ~ "immune response",
                            Description %>% str_detect(cell_cycle) ~ "cell cycling",
                            Description %>% str_detect(infection_response) ~ "response to infection",
                            Description %>% str_detect(inflammation) ~ "inflammation",
                            Description %>% str_detect(injury) ~ "injury response",
                            Description %>% str_detect(metabolic_response) ~ "metabolic response",
                            Description %>% str_detect(exogenous_stimulus) ~ "response to exogenous stimulus",
                            Description %>% str_detect(cellular_regulation) ~ "regulation of cellular processes",
                            Description %>% str_detect(cellular_development) ~ "cellular development",
                            Description %>% str_detect(cell_signalling) ~ "cell signalling"
                        )
                    ) %>%
                    distinct(Description, .keep_all = TRUE) %>%
                    dplyr::pull(group, Description)
            }
        ),
        start_degree = case_when(
            design == "Baseline_vs_Week24" ~ -100,
            design == "Week24_vs_Week52" ~ 50,
            design == "Baseline_vs_Week52" ~ -50,
            TRUE ~ 0
        )
    )


# SAVE CHORD PLOTS ####
with(
    df_plot_02,
    pmap(
        list(design, direction, data_plot, genes_gsea, link_group, start_degree),
        enrichment_chord,
        saveDir = "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/output/"
    )
)
