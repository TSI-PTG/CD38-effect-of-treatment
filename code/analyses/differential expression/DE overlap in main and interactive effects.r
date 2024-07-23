# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(ggvenn) # pak::pak("yanlinlin82/ggvenn")
library(ggrepel) # install.packages("ggrepel")
library(patchwork) # pak::pak("patchwork")
# Custom operators, functions, and datasets
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
# load affymap
load("Z:/DATA/Datalocks/Other data/affymap219_21Oct2019_1306_JR.RData")
# load limma results
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/IQR_filtered_probes_unique_genes_baseline_corrected_cortex_corrected_limma_1208.RData")
limma_tables_te <- limma_tables
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/ONLY_FELZ_IQR_filtered_probes_unique_genes_baseline_corrected_cortex_corrected_limma_1208.RData")
limma_tables_fo <- limma_tables
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/ONLY_PLACEBO_IQR_filtered_probes_unique_genes_baseline_corrected_cortex_corrected_limma_1208.RData")
limma_tables_po <- limma_tables


# load gene lists
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/Hinze_injury_markers.RData")


# WRANGLE THE DATA ####
data <- bind_rows(
    limma_tables_te %>%
        dplyr::select(design, table) %>%
        mutate(contrast = "te", .before = 1),
    limma_tables_fo %>%
        dplyr::select(design, table) %>%
        mutate(contrast = "fo", .before = 1),
    limma_tables_po %>%
        dplyr::select(design, table) %>%
        mutate(contrast = "po", .before = 1)
) %>%
    mutate(
        table = map(
            table,
            function(table) {
                table %>%
                    dplyr::filter(p < 0.05) %>%
                    dplyr::select(AffyID, Symb, logFC, p)
            }
        ),
        probes = map(table, pull, AffyID),
        genes = map(table, pull, Symb)
    )


# SUMMARIZE OVERLAP IN DIFFERENTIALLY EXPRESSED GENES ####
data_overlap_genes <- data %>%
    dplyr::select(contrast, design, genes) %>%
    pivot_wider(names_from = contrast, values_from = genes) %>%
    mutate(
        genes = pmap(
            list(te, fo, po),
            function(te, fo, po) {
                c(te, fo, po) %>% unique()
            }
        ),
        te = pmap(
            list(te, genes),
            function(te, genes) {
                genes %in% te
            }
        ),
        fo = pmap(
            list(fo, genes),
            function(fo, genes) {
                genes %in% fo
            }
        ),
        po = pmap(
            list(po, genes),
            function(po, genes) {
                genes %in% po
            }
        )
    ) %>%
    unnest(everything()) %>%
    nest(.by = "design", venn_data_genes = -design)
# data_overlap_genes %>%
#     unnest(everything()) %>%
#     dplyr::filter(te == TRUE & fo == TRUE) %>%
#     print(n = "all")


# MAKE VENN DIAGRAMS ####
plot_venn_genes <- data_overlap_genes %>%
    mutate(
        venn = pmap(
            list(design, venn_data_genes),
            function(design, venn_data_genes) {
                gene_overlap <- venn_data_genes %>%
                    dplyr::filter(te == TRUE & fo == TRUE & po == TRUE) %>%
                    dplyr::select(genes) %>%
                    dplyr::mutate(x = 0, y = -1)
                venn_data_genes %>%
                    ggplot(aes(A = po, B = fo, C = te)) +
                    # ggplot2::geom_segment(
                    #     mapping = ggplot2::aes(x = 0, xend = 0, y = 0, yend = -1),
                    #     col = "grey70", linewidth = 0.1
                    # ) +
                    # ggplot2::geom_point(
                    #     inherit.aes = FALSE,
                    #     mapping = ggplot2::aes(x = 0, y = -1),
                    #     col = "grey30", size = 0.25
                    # ) +
                    # ggplot2::geom_point(
                    #     inherit.aes = FALSE,
                    #     mapping = ggplot2::aes(x = 0, y = 0),
                    #     col = "grey30", size = 0.25
                    # ) +
                    ggvenn::geom_venn(
                        show_percentage = FALSE,
                        set_names = c(A = "Temporal effect\nin placebo patients", B = "Temporal effect\nin felzartamab patients", C = "Interactive effect of time and treatment"),
                        set_name_size = 3,
                        text_size = 5
                    ) +
                    # ggrepel::geom_label_repel(
                    #     seed = 42,
                    #     inherit.aes = FALSE,
                    #     data = gene_overlap,
                    #     mapping = ggplot2::aes(x = x, y = y, label = genes),
                    #     max.overlaps = Inf,
                    #     nudge_y = -1,
                    #     size = 2.5,
                    #     label.padding = 0.1,
                    #     box.padding = 0.1,
                    #     point.padding = 0,
                    #     label.size = 0.1,
                    #     segment.size = 0.1,
                    #     segment.color = "grey70",
                    #     force = 10
                    # ) +
                    labs(title = design %>% str_replace_all("_", " ")) +
                    theme_void() +
                    theme(
                        plot.title = element_text(hjust = 0.5, face = "bold.italic"),
                        panel.grid = element_blank(),
                        axis.title = element_blank(),
                        axis.text = element_blank(),
                        axis.ticks = element_blank()
                    )
            }
        )
    )


# MAKE PLOT PANELS ####
plot_panels_genes <- wrap_plots(
    A = plot_spacer(),
    B = plot_venn_genes$venn[[1]],
    C = plot_spacer(),
    D = plot_venn_genes$venn[[2]],
    E = plot_venn_genes$venn[[3]],
    design = c(
        "ABBC
        DDEE"
    )
) 
# +
#     plot_annotation(tag_levels = "A")


# SAVE THE PLOTS ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/output/"
ggsave(
    filename = paste(saveDir, "Venn diagram DE genes main vs interaction genes2.png"),
    plot = plot_panels_genes,
    dpi = 300,
    width = 20,
    height = 20,
    units = "cm",
    bg = "white"
)
