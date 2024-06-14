# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(ggrepel) # install.packages("ggrepel")
library(ggpubr) # install.packages("ggpubr")
library(patchwork) # install.packages("patchwork")
library(readxl) # install.packages("readxl")
# Custom operators, functions, and datasets
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
source("C:/R/CD38-effect-of-treatment/code/functions/plot.gg_volcano_timeseries.r")
# load reference data
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/all probes limma 1208.RData")
# load SCC data
simplefile <- read_excel("Z:/MISC/Phil/AA All papers in progress/A GC papers/0000 simple XL files/Kidney 5086/MASTER COPY K5086 SimpleCorrAAInjRej 5AAInjNR 7AARej.xlsx")


# DEFINE PROBE CORRELATIONS WITH EABMR ####
probes_eabmr <- simplefile %>%
    dplyr::select(Affy, "pvalRej7AA4-EABMR") %>%
    rename(
        AffyID = Affy,
        EABMRcorrp = `pvalRej7AA4-EABMR`
    ) %>%
    mutate(ABMRrank = EABMRcorrp %>% rank()) %>%
    arrange(EABMRcorrp)

selected <- simplefile %>%
    dplyr::select(Affy, SYMB, "pvalRej7AA4-EABMR") %>%
    rename(
        AffyID = Affy,
        EABMRcorrp = `pvalRej7AA4-EABMR`
    ) %>%
    arrange(EABMRcorrp) %>%
    distinct(SYMB, .keep_all = TRUE) %>%
    slice(1:20) %>%
    pull(AffyID)


# IDENTIY GENES SUPPRESSED FROM BASELINE TO WEEK24 ####
probes_supressed <- limma_tables %>%
    dplyr::filter(design == "Baseline_vs_Week24") %>%
    pull(toptable) %>%
    pluck(1) %>%
    dplyr::filter(logFC < 0 & P.Value < 0.05) %>%
    # dplyr::filter(logFC < 0) %>%
    pull(AffyID)

probes_increased <- limma_tables %>%
    dplyr::filter(design == "Baseline_vs_Week24") %>%
    pull(toptable) %>%
    pluck(1) %>%
    dplyr::filter(logFC > 0 & P.Value < 0.05) %>%
    # dplyr::filter(logFC < 0) %>%
    pull(AffyID)


# DEFINE COLORS FOR TOP PROBES ####
col_up <- "#7accf9"
col_dn <- "#82fcbf"
col_null <- "grey20"


# WRANGLE THE DATA FOR PLOTTING ####
data_plot <- limma_tables %>%
    dplyr::select(design, table) %>%
    dplyr::filter(design != "Baseline_vs_Week52") %>%
    unnest(everything()) %>%
    left_join(probes_eabmr, by = "AffyID") %>%
    mutate(
        ABMRrank = EABMRcorrp %>% rank(),
        col = case_when(
            AffyID %in% probes_supressed ~ col_up,
            AffyID %in% probes_increased ~ col_dn,
            TRUE ~ col_null
        ),
        alpha = ifelse(col == col_null, 0.5, 0.8),
        order = col %>% factor(levels = c(col_null, col_up, col_dn))
    ) %>%
    arrange(order)

min_p <- data_plot %>%
    slice_min(`<U+0394><U+0394> p`) %>%
    pull(`<U+0394><U+0394> p`) %>%
    log10() * -1.2

data_plot2 <- data_plot %>%
    mutate(
        p = ifelse(
            design == "Week24_vs_Week52",
            -log10(`<U+0394><U+0394> p`) + min_p,
            -log10(`<U+0394><U+0394> p`)
        ),
        logFC = `<U+0394><U+0394> logFC`
    )

data_plot2_trimmed <- data_plot2 %>%
    dplyr::filter(AffyID %in% c(probes_supressed, probes_increased)) %>%
    dplyr::select(AffyID, p, logFC, col, alpha) %>%
    arrange(AffyID)


# MAKE VOLCANO PLOT ####
plot_volcano2 <- data_plot2 %>%
    gg_volcano_timeseries(
        data_plot2_trimmed,
        xbreak1 = 3.7,
        xbreak2 = 3.85
    )

plot_volcano <- ggarrange(plot_volcano2) %>%
    ggpubr::annotate_figure(
        top = text_grob(
            "Genome-wide Effect of Felzartamab Treatment and Refractory Response Post-Treatment",
            face = "bold.italic", size = 15, hjust = 0.51
        )
    ) %>% suppressWarnings()


# SAVE THE PLOT ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/output/"
ggsave(
    plot_volcano,
    file = paste(saveDir, "volcano timeseries all probes 1208.png", sep = ""),
    dpi = 300,
    width = 23,
    height = 14,
    units = "cm",
    bg = "white"
) %>% suppressWarnings()
