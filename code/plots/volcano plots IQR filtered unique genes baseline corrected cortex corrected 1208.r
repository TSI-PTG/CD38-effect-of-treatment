# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(ggrepel) # install.packages("ggrepel")
library(ggpubr) # install.packages("ggpubr")
library(patchwork) # install.packages("patchwork")
library(readxl) # install.packages("readxl")
# Custom operators, functions, and datasets
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
source("C:/R/CD38-effect-of-treatment/code/functions/plot.gg_volcano.r")
source("C:/R/CD38-effect-of-treatment/code/functions/plot.gg_volcano_timeseries.r")
# load reference data
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/IQR_filtered_probes_unique_genes_baseline_corrected_cortex_corrected_limma_1208.RData")
# load gene lists
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/genes_NK_GEP.RData")
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/NK_genes_TBB896.RData")
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/ABMR_endothelial_genes.RData")
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/ABMR_activity_genes.RData")
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/IFNG_genes.RData")



# DEFINE PROBE CORRELATIONS WITH EABMR ####
genes_NK_ATAGC_U133 <- genes_NK_GEP %>%
    dplyr::filter(panel == "ATAGC_U133") %>%
    pull(data) %>%
    pluck(1) %>%
    dplyr::select(AffyID, Symb)

genes_NK_KTB18_RNAseq <- genes_NK_GEP %>%
    dplyr::filter(panel == "KTB18_RNAseq") %>%
    pull(data) %>%
    pluck(1) %>%
    dplyr::select(AffyID, Symb)

genes_NK_LM22_U133 <- genes_NK_GEP %>%
    dplyr::filter(panel == "LM22_U133") %>%
    pull(data) %>%
    pluck(1) %>%
    dplyr::select(AffyID, Symb)


labels_probes <- list(
    genes_IFNG %>%
        dplyr::select(AffyID) %>%
        mutate(label = "IFNG-inducible ABMR activity genes"),
    genes_NK_TBB896 %>% dplyr::select(AffyID) %>%
        mutate(label = "NK cell expressed ABMR activity genes"),
    genes_ABMR_endothelial %>%
        dplyr::select(AffyID) %>%
        mutate(label = "ABMR-associated endothelial genes")
)


# WRANGLE THE DATA FOR PLOTTING Baseline - Week52 ####
data_plot <- limma_tables %>%
    mutate(
        data_plot = map(
            table,
            function(table) {
                table %>%
                    mutate(
                        p = -log10(`<U+0394><U+0394> p`),
                        logFC = `<U+0394><U+0394> logFC`
                    ) %>%
                    suppressWarnings()
            }
        )
    )


# MAKE SINGLE VOLCANO PLOTS ####
plot_volcano <- data_plot %>%
    mutate(
        plot_volcano = map2(
            data_plot, design,
            gg_volcano,
            xlim = c(0, 3),
            ylim = c(-1, 1.15),
            x_break = 1.75,
            # show_labels = TRUE,
            labels_probes = labels_probes,
            # labels_probes_names = c("NK genes", "IFNG genes"),
            col_labels = c("#409f40", "#a1a128", "#167d86"),
            point_size_null = 1.25,
            point_size = 2.5,
            labels_probes_size = 3,
        )
    ) %>%
    dplyr::select(design, plot_volcano)
# plot_volcano$plot_volcano[[1]]




# SAVE THE PLOT DATA ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/"
save(plot_volcano, file = paste(saveDir, "Volcano plots IQR_filtered_probes_unique_genes_baseline_corrected_cortex_corrected 1208.RData", sep = ""))



