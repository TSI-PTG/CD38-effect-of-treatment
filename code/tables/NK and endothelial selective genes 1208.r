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
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/ABMR_activity_genes.RData")


# FILTER THE GENE TABLES ####
gene_tables <- limma_tables %>%
    dplyr::select(design, table) %>%
    expand_grid(geneset = c("ABMR_activity", "NK", "Endothelial")) %>%
    mutate(
        genes = case_when(
            geneset == "ABMR_activity" ~ genes_ABMR_activity$AffyID %>% list(),
            geneset == "NK" ~ genes_NK$AffyID %>% list(),
            geneset == "Endothelial" ~ genes_ABMR_endothelial$AffyID %>% list()
        )
    ) %>%
    relocate(geneset, genes, .after = design) %>%
    mutate(
        gene_tables = pmap(
            list(table, genes),
            function(table, genes) {
                table %>%
                    dplyr::filter(AffyID %in% genes) %>%
                    dplyr::slice_min(`<U+0394><U+0394> p`, by = "Symb")
            }
        )
    )

names(gene_tables$gene_tables) <-gene_tables$design
gene_tables$gene_tables


# EXPORT THE DATA AS .RData FILE ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/"
save(gene_tables, file = paste(saveDir, "gene_tables_limma_1208.RData", sep = ""))


# SPLIT GENE TABLES BY GENESET FOR EXPORTING TABLES ####
gene_tables_ABMR_activity <- gene_tables %>% dplyr::filter(geneset == "ABMR_activity")
gene_tables_NK <- gene_tables %>% dplyr::filter(geneset == "NK")
gene_tables_Endothelial <- gene_tables %>% dplyr::filter(geneset == "Endothelial")



# EXPORT THE DATA AS AN EXCEL SHEET ####
saveDir1 <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/output/"
openxlsx::write.xlsx(gene_tables_ABMR_activity$gene_tables,
    asTable = TRUE,
    file = paste(saveDir1, "ABMR_activity_genes_limma_1208_12Jun24",
        # Sys.Date(),
        # format(Sys.time(), "_%I%M%p"),
        ".xlsx",
        sep = ""
    )
)
openxlsx::write.xlsx(gene_tables_NK$gene_tables,
    asTable = TRUE,
    file = paste(saveDir1, "NK_genes_limma_1208_12Jun24",
        # Sys.Date(),
        # format(Sys.time(), "_%I%M%p"),
        ".xlsx",
        sep = ""
    )
)
openxlsx::write.xlsx(gene_tables_Endothelial$gene_tables,
    asTable = TRUE,
    file = paste(saveDir1, "Endothelial_genes_limma_1208_12Jun24",
        # Sys.Date(),
        # format(Sys.time(), "_%I%M%p"),
        ".xlsx",
        sep = ""
    )
)




