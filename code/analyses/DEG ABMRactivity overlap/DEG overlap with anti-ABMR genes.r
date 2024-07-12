# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(flextable) # install.packages("flextable")
library(officer) # install.packages("officer")
library(openxlsx) # install.packages("openxlsx")
# Custom operators, functions, and datasets
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
# load affymap
load("Z:/DATA/Datalocks/Other data/affymap219_21Oct2019_1306_JR.RData")
# load limma results
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/IQR_filtered_probes_unique_genes_baseline_corrected_cortex_corrected_limma_1208.RData")
# load gene lists
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/ABMRactivity_geneset_probe_scc_K5086.RData")





# REDUCE TABLES TO ONLY INCREASED GENES ####
limma_tables_DEG <- limma_tables %>%
    dplyr::select(design, table) %>%
    mutate(table = map(table, filter, p < 0.05, logFC > 0))


# WRANGLE THE TOP GENES ####
geneset_map <- data_scc %>%
    dplyr::select(score, top100_negative) %>%
    unnest(everything()) %>%
    dplyr::select(score, AffyID) %>%
    nest(.by = AffyID) %>%
    mutate(
        annotation = map_chr(
            data,
            function(data) {
                data %>%
                    pull() %>%
                    str_remove("_SCC") %>%
                    paste(collapse = ",")
            }
        )
    ) %>%
    dplyr::select(-data)



# ANNOTATE THE GENES BY MEMBERSHIP IN TOP100 NEGATIVE GENES BY ABMR ACTIVITY GENESETS ####
DEG_annotation <- limma_tables_DEG %>%
    mutate(
        table = map(
            table,
            function(table) {
                table %>%
                    dplyr::select(-t) %>%
                    left_join(geneset_map, by = "AffyID") %>%
                    relocate(annotation, .after = "PBT")
            }
        )
    )


# QUICK TEST TO SEE ALL OVERLAP ####
deg <- limma_tables_DEG %>%
    dplyr::select(table) %>%
    unnest(table) %>%
    pull(AffyID)


deg %in% geneset_map$AffyID



DEG_annotation$table[[1]] %>%
    print(n="all")


DEG_annotation$table[[3]] %>%
    print(n = "all")
