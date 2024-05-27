enrichment_chord <- function(
    design, data, genes, link_group, seed = 42,
    track_width = 3, track_nudge = -1.2, sector_title_text_size = 0.75,
    title = "differentially expressed genes",
    filename = paste(design, "tmp.png", sep = ""), saveDir = NULL) {
    require(clusterProfiler) # pak::pak("YuLab-SMU/clusterProfiler")
    require(circlize) # install.packages("circlize") pak::pak("jokergoo/circlize")

    # DEFINE SECTORS ####
    sector_genes <- data %>%
        dplyr::distinct(Genes) %>%
        dplyr::pull(Genes)

    sector_cellcycle <- data %>%
        dplyr::filter(group == "cell cycling") %>%
        dplyr::distinct(Description, .keep_all = TRUE) %>%
        dplyr::pull(Description)

    sector_immune <- data %>%
        dplyr::filter(group == "immune response") %>%
        dplyr::distinct(Description, .keep_all = TRUE) %>%
        dplyr::pull(Description)

    sector_metabolism <- data %>%
        dplyr::filter(group == "metabolic response") %>%
        dplyr::distinct(Description, .keep_all = TRUE) %>%
        dplyr::pull(Description)

    sector_infection <- data %>%
        dplyr::filter(group == "response to infection") %>%
        dplyr::distinct(Description, .keep_all = TRUE) %>%
        dplyr::pull(Description)

    sector_inflammation <- data %>%
        dplyr::filter(group == "inflammation") %>%
        dplyr::distinct(Description, .keep_all = TRUE) %>%
        dplyr::pull(Description)

    sector_injury <- data %>%
        dplyr::filter(group == "injury response") %>%
        dplyr::distinct(Description, .keep_all = TRUE) %>%
        dplyr::pull(Description)

    sector_stimulus <- data %>%
        dplyr::filter(group == "response to external stimilus") %>%
        dplyr::distinct(Description, .keep_all = TRUE) %>%
        dplyr::pull(Description)

    sector_cellregulation <- data %>%
        dplyr::filter(group == "regulation of cellular processes") %>%
        dplyr::distinct(Description, .keep_all = TRUE) %>%
        dplyr::pull(Description)

    sector_celldevelopment <- data %>%
        dplyr::filter(group == "cellular development") %>%
        dplyr::distinct(Description, .keep_all = TRUE) %>%
        dplyr::pull(Description)

    sector_cellcommunication <- data %>%
        dplyr::filter(group == "cellular communication") %>%
        dplyr::distinct(Description, .keep_all = TRUE) %>%
        dplyr::pull(Description)


    # DEFINE SECTOR COLOURS ####
    cols_group <- c(
        "genes" = "grey20",
        "immune" = "#f04141",
        "infection" = "#ff9040",
        "injury" = "#8a41f0",
        "inflammation" = "#ffe083",
        "cellcycle" = "#4a4aff",
        "stimulus" = "#7cfe7c",
        "metabolism" = "#7cfe7c",
        "cellregulation" = "#fe7ce4",
        "development" = "#62ff4a",
        "cellcommunication" = "#fe7ca9"
    )


    # DEFINE GO PATHWAY COLOURS ####
    df_cols <- data %>%
        tidyr::pivot_longer(cols = c("Genes", "Description")) %>%
        dplyr::distinct(group, value, .keep_all = TRUE) %>%
        tidyr::nest(.by = c("name", "group")) %>%
        dplyr::mutate(
            hue = dplyr::case_when(
                name == "Genes" ~ "grey20",
                group == "cell cycling" ~ "blue",
                group == "immune response" ~ "red",
                group == "response to infection" ~ "orange",
                group == "inflammation" ~ "yellow",
                group == "injury response" ~ "purple",
                group == "response to external stimilus" ~ "green",
                group == "metabolic response" ~ "green",
                group == "regulation of cellular processes" ~ "pink",
                group == "cellular development" ~ "green",
                group == "cellular communication" ~ "pink",
            ),
            data = purrr::pmap(
                list(name, group, data, hue),
                function(name, group, data, hue) {
                    n <- data %>% nrow()
                    if (name == "Genes") {
                        cols <- "grey20"
                    } else if (name != "Genes") {
                        set.seed(seed)
                        cols <- circlize::rand_color(
                            n = n,
                            hue = hue
                            # transparency = 0.9
                        )
                    }
                    data$col <- cols
                    data
                }
            )
        ) %>%
        dplyr::pull(data) %>%
        dplyr::bind_rows()

    grid_col <- data %>%
        tidyr::pivot_longer(cols = c("Genes", "Description")) %>%
        dplyr::left_join(df_cols, by = "value") %>%
        dplyr::pull(col, value)

    link_col <- data %>%
        dplyr::left_join(df_cols %>% dplyr::rename(Description = value), by = "Description") %>%
        dplyr::pull(col)

    # PLOTTING GLOBALS ####
    n_genes <- genes %>%
        length()
    n_genes_enriched <- data %>%
        dplyr::distinct(Genes) %>%
        nrow()


    # SAVE PLOT ####
    png(
        paste(saveDir, filename),
        units = "cm", width = 15, height = 15, res = 600,
        bg = "transparent"
    )
    circlize::circos.clear()
    circlize::circos.par(
        circle.margin = c(0.01, 0.225),
        start.degree = -11,
        track.height = 1
    )
    # MAKE CHORD PLOT ####
    data %>%
        circlize::chordDiagram(
            group = link_group,
            col = link_col,
            grid.col = grid_col,
            transparency = 0.25,
            annotationTrack = "grid",
            # preAllocateTracks = list(
            # transparency = 0.8,
            # link.lwd = 1,
            # link.lty = 1,
            # link.border = 1,
            # track.height = 1
            # ),
            big.gap = 10,
            small.gap = 0.1,
            scale = FALSE
        )
    # title(title, cex = 0.6, family = "TT Arial")


    # HIGHLIGHT SECTORS ####
    # genes
    circlize::highlight.sector(
        sector.index = sector_genes,
        track.index = 1,
        col = cols_group["genes"],
        text = paste("enriched genes (n = ", n_genes_enriched, " out of ", n_genes, ")", sep = ""),
        text.col = "white",
        facing = "bending.inside",
        # niceFacing = TRUE,
        cex = 1.5,
        padding = c(track_nudge, 0, track_width, 0)
    )
    circlize::highlight.sector(
        sector.index = sector_genes,
        track.index = 1,
        col = NA,
        text = title,
        text.col = "black",
        facing = "bending.outside",
        niceFacing = TRUE,
        cex = 2,
        padding = c(track_nudge, 0, track_width, 0),
        text.vjust = "10mm"
    )
    # immune
    if (length(sector_immune) > 0) {
        circlize::highlight.sector(
            sector.index = sector_immune,
            track.index = 1,
            col = cols_group["immune"],
            padding = c(track_nudge, 0, track_width, 0),
            text = "immune response",
            cex = sector_title_text_size,
            text.col = "black",
            facing = "bending.outside",
            niceFacing = TRUE,
            text.vjust = "7mm"
        )
    }
    # infection
    if (length(sector_infection) > 0) {
        circlize::highlight.sector(
            sector.index = sector_infection,
            track.index = 1,
            col = cols_group["infection"],
            padding = c(track_nudge, 0, track_width, 0),
            text = "resp. to infection",
            cex = sector_title_text_size,
            text.col = "black",
            facing = "bending.outside",
            niceFacing = TRUE,
            text.vjust = "7mm"
        )
    }
    # injury
    if (length(sector_injury) > 0) {
        circlize::highlight.sector(
            sector.index = sector_injury,
            track.index = 1,
            col = cols_group["injury"],
            padding = c(track_nudge, 0, track_width, 0),
            text = "resp. to injury",
            cex = sector_title_text_size,
            text.col = "black",
            facing = "bending.outside",
            niceFacing = TRUE,
            text.vjust = "7mm"
        )
    }
    # metabolism
    if (length(sector_metabolism) > 0) {
        circlize::highlight.sector(
            sector.index = sector_metabolism,
            track.index = 1,
            col = cols_group["metabolism"],
            padding = c(track_nudge, 0, track_width, 0),
            text = "metabolic response",
            cex = sector_title_text_size,
            text.col = "black",
            facing = "bending.outside",
            niceFacing = TRUE,
            text.vjust = "7mm"
        )
    }
    # cellcycle
    if (length(sector_cellcycle) > 0) {
        circlize::highlight.sector(
            sector.index = sector_cellcycle,
            track.index = 1,
            col = cols_group["cellcycle"],
            padding = c(track_nudge, 0, track_width, 0),
            text = "cell cycle",
            cex = sector_title_text_size,
            text.col = "black",
            facing = "bending.outside",
            niceFacing = TRUE,
            text.vjust = "7mm"
        )
    }
    # inflammation
    if (length(sector_inflammation) > 0) {
        circlize::highlight.sector(
            sector.index = sector_inflammation,
            track.index = 1,
            col = cols_group["inflammation"],
            padding = c(track_nudge, 0, track_width, 0),
            text = "inflammation",
            cex = sector_title_text_size,
            text.col = "black",
            facing = "bending.outside",
            niceFacing = TRUE,
            text.vjust = "7mm"
        )
    }
    # stimulus
    if (length(sector_stimulus) > 0) {
        circlize::highlight.sector(
            sector.index = sector_stimulus,
            track.index = 1,
            col = cols_group["stimulus"],
            padding = c(track_nudge, 0, track_width, 0),
            text = "resp. to external stimilus",
            cex = sector_title_text_size,
            text.col = "black",
            facing = "bending.outside",
            niceFacing = TRUE,
            text.vjust = "7mm"
        )
    }
    # cellregulation
    if (length(sector_celldevelopment) > 0) {
        circlize::highlight.sector(
            sector.index = sector_celldevelopment,
            track.index = 1,
            col = cols_group["development"],
            padding = c(track_nudge, 0, track_width, 0),
            text = "cellular development",
            cex = sector_title_text_size,
            text.col = "black",
            facing = "bending.outside",
            niceFacing = TRUE,
            text.vjust = "7mm"
        )
    }
    # cellregulation
    if (length(sector_cellregulation) > 0) {
        circlize::highlight.sector(
            sector.index = sector_cellregulation,
            track.index = 1,
            col = cols_group["cellregulation"],
            padding = c(track_nudge, 0, track_width, 0),
            text = "reg. of cell. proc.",
            cex = sector_title_text_size,
            text.col = "black",
            facing = "bending.outside",
            niceFacing = TRUE,
            text.vjust = "7mm"
        )
    }
    # cellcommunication
    if (length(sector_cellcommunication) > 0) {
        circlize::highlight.sector(
            sector.index = sector_cellcommunication,
            track.index = 1,
            col = cols_group["cellcommunication"],
            padding = c(track_nudge, 0, track_width, 0),
            text = "cellular communication",
            cex = sector_title_text_size,
            text.col = "black",
            facing = "bending.outside",
            niceFacing = TRUE,
            text.vjust = "7mm"
        )
    }


    dev.off()
}
