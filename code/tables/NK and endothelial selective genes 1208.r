# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(flextable) # install.packages("flextable") #for table outputs
library(officer) # install.packages("officer")
library(openxlsx) # install.packages("openxlsx")
# Bioconductor libraries
library(Biobase) # BiocManager::install("Biobase")
# Custom operators, functions, and datasets
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
# load reference set
# load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/data_expressionset_k1208.RData")
# load affymap
load("Z:/DATA/Datalocks/Other data/affymap219_21Oct2019_1306_JR.RData")
# load limma results
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/all_probes_limma_1208.RData")
# load gene lists
# load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/NK cell selective genes.RData")
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/NK_genes_L765.RData")
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/ABMR_endothelial_genes.RData")


# FILTER THE GENE TABLES ####
gene_tables <- limma_tables %>%
    dplyr::select(design, table) %>%
    mutate(
        tables_nk = map(
            table,
            function(table) {
                table %>%
                    dplyr::filter(AffyID %in% genes_NK$AffyID) %>%
                    dplyr::slice_min(`<U+0394><U+0394> p`, by = "Symb")
            }
        ),
        tables_endothelial = map(
            table,
            function(table) {
                table %>%
                    dplyr::filter(AffyID %in% genes_ABMR_endothelial$AffyID) %>%
                    dplyr::slice_min(`<U+0394><U+0394> p`, by = "Symb")
            }
        )
    )
names(gene_tables$tables_nk) <- gene_tables$design
names(gene_tables$tables_endothelial) <- gene_tables$design


gene_tables$tables_nk[[1]]
gene_tables$tables_endothelial[[1]]


# EXPORT THE DATA AS .RData FILE ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/"
save(gene_tables, file = paste(saveDir, "gene_tables_limma_1208.RData", sep = ""))


# EXPORT THE DATA AS AN EXCEL SHEET ####
saveDir1 <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/output/"
openxlsx::write.xlsx(gene_tables$tables_nk,
    asTable = TRUE,
    file = paste(saveDir1, "NK_genes_limma_1208_11Jun24",
        # Sys.Date(),
        # format(Sys.time(), "_%I%M%p"),
        ".xlsx",
        sep = ""
    )
)
openxlsx::write.xlsx(gene_tables$tables_endothelial,
    asTable = TRUE,
    file = paste(saveDir1, "Endothelial_genes_limma_1208_11Jun24",
        # Sys.Date(),
        # format(Sys.time(), "_%I%M%p"),
        ".xlsx",
        sep = ""
    )
)