# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(gghalves) # install.packages("gghalves")
library(ggrepel)
# Bioconductor libraries
library(Biobase) # BiocManager::install("Biobase")
library(limma) # BiocManager::install("limma")
# Custom operators, functions, and datasets
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
# load reference set
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0Georg Felz CD38 Vienna/G_Rstuff/data/vienna_5086_6Mar24.RData")
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/data_expressionset_k1208.RData")
# load affymap
load("Z:/DATA/Datalocks/Other data/affymap219_21Oct2019_1306_JR.RData")
# load cortex estimator
load("Z:/Genome-Archive/RefData/Kidney Report 2023 for Kashi/data/glm.Cort.RData")


# DEFINE THE SET ####
set <- data_expressionset_k1208[, data_expressionset_k1208$Patient %nin% c(15, 18)]


# WRANGLE THE PHENOTYPE DATA ####
data_phenotype <- set %>%
    pData()


# WRANGLE THE EXPRESSION DATA ####
data_exprs <- set %>%
    exprs() %>%
    as_tibble(rownames = "AffyID") %>%
    right_join(affymap219 %>% dplyr::select(AffyID, Symb), ., by = "AffyID") %>%
    tibble()


# ISOLATE PODOCIN EXPRESSION DATA ####
data_exprs_podocin <- data_exprs %>%
    dplyr::filter(Symb == "NPHS2") %>%
    pivot_longer(
        cols = -1:-2,
        names_to = "CEL",
        values_to = "podocin"
    )


# JOIN THE EXPRESSION DATA WITH PHENOTYPE DATA ####
data_podocin <- data_exprs_podocin %>%
    left_join(data_phenotype, by = "CEL") %>%
    dplyr::filter(Group != "FU1b", Patient %nin% c(15, 18)) %>%
    mutate(
        cortex = predict.glm(glm.Cort, data.frame(NPHS2 = data_podocin$podocin), type = "response"),
        .after = podocin
    )

data_geomean <- data_podocin %>%
    summarize(mean = 2^(podocin %>% mean()), .by = c("Felzartamab_Followup", "Followup", "Felzartamab"))



# PLOT THE DISTRIBUTION OF PODOCIN ####
plot_podocin <- data_podocin %>%
    ggplot(mapping = aes(x = Followup, y = 2^podocin)) +
    geom_point(
        inherit.aes = FALSE,
        mapping = aes(x = Followup, y = mean),
        data = data_geomean, col = "red", size = 3
    ) +
    geom_point(
        # position = position_jitter(0.3, seed = 42),
        # shape = 21, fill = "grey95", size = 8
    ) +
    # geom_text(
    #     mapping = aes(label = Patient),
    #     position = position_jitter(0.3, seed = 42),
    #     # min.segment.length = 0.01,
    #     # nudge_x = 0.1,
    #     # direction = "x"
    # ) +
    geom_label_repel(
        mapping = aes(label = Patient),
        # position = position_jitter(0.3, seed = 42),
        min.segment.length = 0.01,
        # nudge_x = 0.1,
        # direction = "x"
    ) +
    labs(y = "Podocin expression (NPHS2)") +
    coord_cartesian(ylim = c(0, NA)) +
    theme_bw() +
    facet_wrap(~Felzartamab)


# PLOT THE DISTRIBUTION OF cortex ####
plot_cortex <- data_podocin %>%
    ggplot(mapping = aes(x = Followup, y = cortex)) +
    geom_point(
        # position = position_jitter(0.3, seed = 42),
        # shape = 21, fill = "grey95", size = 8
    ) +
    # geom_text(
    #     mapping = aes(label = Patient),
    #     position = position_jitter(0.3, seed = 42),
    #     # min.segment.length = 0.01,
    #     # nudge_x = 0.1,
    #     # direction = "x"
    # ) +
    geom_label_repel(
        mapping = aes(label = Patient),
        # position = position_jitter(0.3, seed = 42),
        min.segment.length = 0.01,
        max.overlaps = Inf
        # nudge_x = 0.1,
        # direction = "x"
    ) +
    labs(y = "Estimate cortex (%)") +
    coord_cartesian(ylim = c(0, NA)) +
    theme_bw() +
    facet_wrap(~Felzartamab)



# SAVE THE PLOTS ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/output/"
ggsave(
    filename = paste(saveDir, "podocin.png"),
    plot = plot_podocin,
    dpi = 300,
    width = 20,
    height = 10,
    units = "cm",
    bg = "white"
)
ggsave(
    filename = paste(saveDir, "cortex.png"),
    plot = plot_cortex,
    dpi = 300,
    width = 20,
    height = 10,
    units = "cm",
    bg = "white"
)
