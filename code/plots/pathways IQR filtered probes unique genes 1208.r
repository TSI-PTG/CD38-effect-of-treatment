# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(clusterProfiler) # pak::pak("YuLab-SMU/clusterProfiler")
library(ggpubr) # install.packages("ggpubr")
library(patchwork) # install.packages("patchwork")
# Custom operators, functions, and datasets
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
# load reference data
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/felzartamab_gsea_k1208.RData")


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



# WRANGLE DATA FOR PLOTTING ####
df_plot_00 <- felzartamab_gsea_k1208 %>%
    dplyr::select(design, genes_gsea, gsea_go_tables) %>%
    mutate(
        data_plot = map(
            gsea_go_tables,
            function(gsea_go_tables) {
                gsea_go_tables %>%
                    slice_min(pvalue, n = 20) %>%
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
                            Description %>% str_detect(external_stimulus) ~ "response to external stimilus",
                            Description %>% str_detect(reg_cellular_processes) ~ "regulation of cellular processes",
                            Description %>% str_detect(cellular_development) ~ "cellular development",
                            Description %>% str_detect(cellular_communication) ~ "cellular communication",
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
    )


# DEFINE GROUPS FOR PATHWAYS ####
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
                            Description %>% str_detect(external_stimulus) ~ "response to external stimilus",
                            Description %>% str_detect(reg_cellular_processes) ~ "regulation of cellular processes",
                            Description %>% str_detect(cellular_development) ~ "cellular development",
                            Description %>% str_detect(cellular_communication) ~ "cellular communication"
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
                            Description %>% str_detect(external_stimulus) ~ "response to external stimilus",
                            Description %>% str_detect(reg_cellular_processes) ~ "regulation of cellular processes",
                            Description %>% str_detect(cellular_development) ~ "cellular development",
                            Description %>% str_detect(cellular_communication) ~ "cellular communication"
                        )
                    ) %>%
                    distinct(Description, .keep_all = TRUE) %>%
                    dplyr::pull(group, Description)
            }
        )
    )
df_plot_02$link_group[[1]][df_plot_02$link_group[[1]] %>% is.na()]
df_plot_02$data_plot[[1]]




# SAVE CHORD PLOTS ####
with(
    df_plot_02,
    pmap(
        list(design, data_plot, genes_gsea, link_group),
        enrichment_chord,
        saveDir = "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/output/"
    )
)



enrichment_chord(
    data = df_plot_02$data_plot[[1]],
    genes = df_plot_02$genes_gsea[[1]],
    link_group = df_plot_02$link_group[[1]],
    # fileName = "",
    saveDir = "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/output/"
)
