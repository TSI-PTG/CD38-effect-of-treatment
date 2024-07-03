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
limma_tables_te <- limma_tables
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/ONLY_FELZ_IQR_filtered_probes_unique_genes_baseline_corrected_cortex_corrected_limma_1208.RData")
limma_tables_fo <- limma_tables
# load gene lists
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/Hinze_injury_markers.RData")



# WRANGLE THE DATA ####
limma_tables_te$table[[1]]
limma_tables_fo$table[[1]]





data <- bind_rows(
    limma_tables_te %>%
        dplyr::select(design, table) %>%
        mutate(contrast = "te", .before = 1) %>%
        mutate(table = map(table, rename, p = `<U+0394><U+0394> p`)),
    limma_tables_fo %>% dplyr::select(design, table) %>%
        mutate(contrast = "fo", .before = 1)
) %>%
    mutate(
        table = map(
            table,
            function(table) {
                table %>% dplyr::filter(p < 0.05)
            }
        )
    )



