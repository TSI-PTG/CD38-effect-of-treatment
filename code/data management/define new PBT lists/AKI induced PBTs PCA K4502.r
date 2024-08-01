# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(FactoMineR) # install.packages('FactoMineR')
library(flextable) # install.packages('flextable')
library(patchwork) # install.packages('patchwork')
library(gghalves) # pak::pak("erocoar/gghalves")
# Custom operators, functions, and datasets
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
pbt.generate <- function(data_CEL, reference_set, pbtlist) {
    for (i in 1:length(names(pbtlist))) {
        Biobase::pData(data_CEL)[names(pbtlist)[i]] <- NA
        conavg <- apply(Biobase::exprs(reference_set)[pbtlist[[i]], ], 1, mean)
        smallset <- data_CEL[pbtlist[[i]], ]
        dog <- sweep(Biobase::exprs(smallset), 1, conavg)
        dog <- apply(dog, 2, mean)
        sel <- which(names(Biobase::pData(data_CEL)) == names(pbtlist)[i])
        Biobase::pData(data_CEL)[, sel] <- dog
    }
    return(data_CEL)
}
# load datalock
if (!exists("K5086")) {
    load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/K5086_29jul2024.RData")
}
# load nephrectomies
load("Z:/Genome-Archive/RefData/KidneyReports2021REDCap/data/Conset.RData") # 4 Controls used for pbt calculations
# load new pbt list
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/PBTlist219_AKIinduced_26Jul2024_PTG.RData")


# DEFINE CORE SET ####
set <- K5086[, which(K5086$Cort > 0.1)]
pData(set) <- set %>%
    pData() %>%
    tibble() %>%
    mutate(
        injAA = case_when(
            InjAA5Clust == 1 ~ "MildCKD",
            InjAA5Clust == 2 ~ "CKDAKI",
            InjAA5Clust == 3 ~ "AKI1",
            InjAA5Clust == 4 ~ "AKI2",
            InjAA5Clust == 5 ~ "Normal"
        )
    ) %>%
    dplyr::select(CEL, injAA)


# ADD pbt SCORES TO REFERENCE SET ####
set1 <- pbt.generate(set, Conset, PBTlist219_AKIinduced)


# DEFINE SUBSETS FOR EACH pbt ####
df_set <- tibble(
    pbt = PBTlist219_AKIinduced %>% names(),
    datExpr = map(
        pbt,
        function(pbt) {
            pbt %>% print()
            set[featureNames(set) %in% c(PBTlist219_AKIinduced %>% purrr::pluck(pbt)), ] %>%
                Biobase::exprs() %>%
                t() %>%
                as.data.frame()
        }
    )
)
# df_set$datExpr[[1]]


# FILTER PBTS USING PC METHOD ####
df_pca <- df_set %>%
    mutate(
        pc_scores = map2(
            pbt, datExpr,
            function(pbt, datExpr) {
                datExpr %>%
                    PCA(graph = FALSE) %>%
                    pluck("ind", "coord") %>%
                    as_tibble(rownames = "CEL") %>%
                    rename_at(vars(contains("Dim")), ~ str_replace(., "Dim.", paste(pbt, "PC", sep = "_"))) %>%
                    dplyr::select(CEL, contains("_PC1"), contains("_PC2"))
            }
        )
    )
df_pca$pc_scores[[1]]


# WRANGLE THE NEW PBTS ####
df_pbt_pc <- reduce(
    df_pca %>% pull(pc_scores),
    left_join,
    by = "CEL"
)


# CORRELATE THE PC SCORES WITH pbt SCORES IN ATTEMPT TO QUANTIFY THE DIRECTION OF CHANGE ####
df_pbt_score <- set1 %>%
    pData() %>%
    dplyr::select(-injAA) %>%
    pivot_longer(-CEL, names_to = "pbt", values_to = "pbt_score") %>%
    nest(.by = pbt, pbt_scores = -pbt)

df_pbt_pc_scores <- df_pbt_pc %>%
    pivot_longer(-CEL, names_to = "pbt", values_to = "pc_score") %>%
    dplyr::filter(pbt %>% str_detect("_PC1")) %>%
    nest(.by = pbt, pbt_pc_scores = -pbt) %>%
    mutate(pbt = pbt %>% str_remove("_PC1|_PC2"))

df_scc <- df_pbt_score %>%
    left_join(df_pbt_pc_scores, by = "pbt") %>%
    mutate(data = map2(pbt_scores, pbt_pc_scores, left_join, by = "CEL")) %>%
    mutate(
        SCC = map_dbl(
            data,
            function(data) {
                cor(data$pbt_score, data$pc_score, method = "spearman")
            }
        )
    ) %>%
    mutate_if(is.numeric, ~ formatC(., digits = 3, format = "f"))

df_scc %>% print(n = "all")

# TABULATE NUMBER OF GENES IN EACH pbt ####
n_genes <- df_set %>%
    mutate("n genes in pbt" = map_dbl(datExpr, ncol)) %>%
    dplyr::select(-datExpr)


# WRANGLE THE DATA TO MAKE SUMMARY TABLES ####
set_pbts <- set %>%
    pData() %>%
    left_join(df_pbt_pc, by = "CEL")


set_pbts_means <- set_pbts %>%
    dplyr::select(-CEL) %>%
    nest(.by = injAA) %>%
    mutate(means = map(data, summarise_all, mean)) %>%
    dplyr::select(-data) %>%
    unnest(means) %>%
    pivot_longer(cols = -injAA, names_to = "pbt") %>%
    pivot_wider(names_from = injAA, values_from = value) %>%
    mutate_if(is.numeric, ~ formatC(., digits = 3, format = "f")) %>%
    mutate(
        score = pbt %>% str_extract("PC1|PC2"),
        pbt = pbt %>% str_remove("_PC1|_PC2")
    ) %>%
    left_join(n_genes, by = "pbt") %>%
    left_join(df_scc %>% dplyr::select(pbt, SCC), by = "pbt") %>%
    dplyr::select(pbt, "n genes in pbt", score, SCC, MildCKD, CKDAKI, AKI1, AKI2, Normal)



# DEFINE COLOUR GRADIENT FOR INTERPRETATION ####
# colormat <-



# TABULATE pbt SCORES ####
set_pbts_means %>%
    dplyr::filter(score == "PC1") %>%
    flextable::flextable() %>%
    flextable::add_header_row(
        top = TRUE,
        values = rep("Mean PC score of AKI markers in K4502", flextable::ncol_keys(.))
    ) %>%
    flextable::merge_h(part = "header") %>%
    flextable::merge_v(j = 1, part = "body") %>%
    flextable::border_remove() %>%
    flextable::border(part = "all", border = officer::fp_border()) %>%
    flextable::border(part = "footer", border.left = officer::fp_border(), border.right = officer::fp_border()) %>%
    flextable::border(i = 1, part = "footer", border.bottom = officer::fp_border()) %>%
    flextable::align(align = "center", part = "all") %>%
    flextable::font(fontname = "Arial", part = "all") %>%
    flextable::fontsize(size = 12, part = "all") %>%
    flextable::fontsize(i = 1, size = 15, part = "header") %>%
    flextable::bold(part = "header") %>%
    flextable::bg(bg = "white", part = "all") %>%
    flextable::padding(padding = 0, part = "all") %>%
    flextable::width(., width = dim(.)$widths * 30 / (flextable::flextable_dim(.)$widths), unit = "cm") %>%
    print(preview = "pptx")


# MAKE PLOTS FOR ALL PBTS ####
df_plot <- df_scc %>%
    dplyr::select(pbt, data) %>%
    unnest(data) %>%
    pivot_longer(c(pbt_score, pc_score), names_to = "method", values_to = "score") %>%
    left_join(set1 %>% pData() %>% dplyr::select(CEL, injAA), by = "CEL") %>%
    nest(.by = c(pbt, method))


df_joint_plot <- df_plot %>%
    mutate(
        plot = pmap(
            list(pbt, method, data),
            function(pbt, method, data) {
                dotsize <- case_when(
                    method %>% str_detect("pbt") ~ 1,
                    method %>% str_detect("pc") ~ 10,
                    TRUE ~ 6
                )
                stackratio <- case_when(
                    method %>% str_detect("pbt") ~ 1,
                    method %>% str_detect("pc") ~ 0.5,
                    TRUE ~ 1
                )


                mean_normal <- data %>%
                    dplyr::filter(injAA == "Normal") %>%
                    pull(score) %>%
                    mean()
                method_name <- case_when(
                    method %>% str_detect("pbt") ~ "Method 1\nPBT",
                    method %>% str_detect("pc") ~ "Method 2\nPCA transformation"
                )
                data %>%
                    mutate(
                        method = case_when(
                            method %>% str_detect("pbt") ~ "Method 1\nPBT",
                            method %>% str_detect("pc") ~ "Method 2\nPCA transformation"
                        ),
                        pbt = pbt,
                        group = injAA %>% factor(levels = c("MildCKD", "CKDAKI", "AKI1", "AKI2", "Normal") %>% rev())
                    ) %>%
                    arrange(method) %>%
                    ggplot(aes(x = group, y = score, fill = group)) +
                    geom_hline(yintercept = mean_normal) +
                    geom_half_violin(
                        trim = FALSE,
                        alpha = 0.85,
                        side = "r",
                        scale = "width",
                        linewidth = 0.25
                    ) +
                    # geom_half_dotplot(
                    #     alpha = 0.85,
                    #     dotsize = dotsize,
                    #     stroke = 0.25,
                    #     stackratio = stackratio,
                    #     binwidth = 0.01,
                    #     stackdir = "down",
                    #     show.legend = FALSE,
                    #     col = "grey60"
                    # ) +
                    stat_summary(
                        fun = mean,
                        geom = "point",
                        size = 3,
                        col = "black",
                        fill = "white",
                        pch = 21,
                        stroke = 0.1,
                        show.legend = FALSE
                    ) +
                    stat_summary(
                        aes(col = group),
                        fun = mean,
                        geom = "point",
                        size = 1.5,
                        show.legend = FALSE
                    ) +
                    scale_fill_manual(values = c("#8b008b", "#50C878", "#E67451", "#7b68ee", "grey70") %>% rev()) +
                    scale_colour_manual(values = c("#8b008b", "#50C878", "#E67451", "#7b68ee", "grey70") %>% rev()) +
                    coord_flip() +
                    labs(
                        fill = NULL,
                        y = "score",
                        x = NULL
                    ) +
                    theme_bw() +
                    facet_wrap(~ method + pbt) +
                    theme(
                        panel.grid = element_blank(),
                        # axis.text.y = element_text(colour = "black", size = 8),
                        axis.text.y = element_blank(),
                        axis.ticks.y = element_blank(),
                        legend.position = "top",
                        strip.text = element_text(face = "bold", size = 12)
                    )
            }
        )
    )


df_joint_plot$plot[[2]]


# MAKE PLOT PANELS FOR PT pbt ####
plot_PT_New1 <- df_joint_plot %>%
    dplyr::filter(pbt %>% str_detect("PT_New1")) %>%
    arrange(pbt) %>%
    pull(plot) %>%
    wrap_plots(nrow = 1) +
    plot_layout(guides = "collect") &
    theme(legend.position = "top") &
    guides(fill = guide_legend(nrow = 1, reverse = TRUE))

plot_PT_New2 <- df_joint_plot %>%
    dplyr::filter(pbt %>% str_detect("PT_New2")) %>%
    arrange(pbt) %>%
    pull(plot) %>%
    wrap_plots(nrow = 1) +
    plot_layout(guides = "collect") &
    theme(legend.position = "top") &
    guides(fill = guide_legend(nrow = 1, reverse = TRUE))

plot_PT_New3 <- df_joint_plot %>%
    dplyr::filter(pbt %>% str_detect("PT_New3")) %>%
    arrange(pbt) %>%
    pull(plot) %>%
    wrap_plots(nrow = 1) +
    plot_layout(guides = "collect") &
    theme(legend.position = "top") &
    guides(fill = guide_legend(nrow = 1, reverse = TRUE))

plot_PT_New4 <- df_joint_plot %>%
    dplyr::filter(pbt %>% str_detect("PT_New4")) %>%
    arrange(pbt) %>%
    pull(plot) %>%
    wrap_plots(nrow = 1) +
    plot_layout(guides = "collect") &
    theme(legend.position = "top") &
    guides(fill = guide_legend(nrow = 1, reverse = TRUE))


# MAKE PLOT PANELS FOR TL pbt ####
plot_TL_New1 <- df_joint_plot %>%
    dplyr::filter(pbt %>% str_detect("TL_New1")) %>%
    arrange(pbt) %>%
    pull(plot) %>%
    wrap_plots(nrow = 1) +
    plot_layout(guides = "collect") &
    theme(legend.position = "top") &
    guides(fill = guide_legend(nrow = 1, reverse = TRUE))


# MAKE PLOT PANELS FOR TAL pbt ####
plot_TAL_New2 <- df_joint_plot %>%
    dplyr::filter(pbt %>% str_detect("TAL_New2")) %>%
    arrange(pbt) %>%
    pull(plot) %>%
    wrap_plots(nrow = 1) +
    plot_layout(guides = "collect") &
    theme(legend.position = "top") &
    guides(fill = guide_legend(nrow = 1, reverse = TRUE))

plot_TAL_New3 <- df_joint_plot %>%
    dplyr::filter(pbt %>% str_detect("TAL_New3")) %>%
    arrange(pbt) %>%
    pull(plot) %>%
    wrap_plots(nrow = 1) +
    plot_layout(guides = "collect") &
    theme(legend.position = "top") &
    guides(fill = guide_legend(nrow = 1, reverse = TRUE))

plot_TAL_New4 <- df_joint_plot %>%
    dplyr::filter(pbt %>% str_detect("TAL_New4")) %>%
    arrange(pbt) %>%
    pull(plot) %>%
    wrap_plots(nrow = 1) +
    plot_layout(guides = "collect") &
    theme(legend.position = "top") &
    guides(fill = guide_legend(nrow = 1, reverse = TRUE))


# MAKE PLOT PANELS FOR DCT pbt ####
plot_DCT_New1 <- df_joint_plot %>%
    dplyr::filter(pbt %>% str_detect("DCT_New1")) %>%
    arrange(pbt) %>%
    pull(plot) %>%
    wrap_plots(nrow = 1) +
    plot_layout(guides = "collect") &
    theme(legend.position = "top") &
    guides(fill = guide_legend(nrow = 1, reverse = TRUE))

plot_DCT_New2 <- df_joint_plot %>%
    dplyr::filter(pbt %>% str_detect("DCT_New2")) %>%
    arrange(pbt) %>%
    pull(plot) %>%
    wrap_plots(nrow = 1) +
    plot_layout(guides = "collect") &
    theme(legend.position = "top") &
    guides(fill = guide_legend(nrow = 1, reverse = TRUE))

plot_DCT_New3 <- df_joint_plot %>%
    dplyr::filter(pbt %>% str_detect("DCT_New3")) %>%
    arrange(pbt) %>%
    pull(plot) %>%
    wrap_plots(nrow = 1) +
    plot_layout(guides = "collect") &
    theme(legend.position = "top") &
    guides(fill = guide_legend(nrow = 1, reverse = TRUE))

plot_DCT_New4 <- df_joint_plot %>%
    dplyr::filter(pbt %>% str_detect("DCT_New4")) %>%
    arrange(pbt) %>%
    pull(plot) %>%
    wrap_plots(nrow = 1) +
    plot_layout(guides = "collect") &
    theme(legend.position = "top") &
    guides(fill = guide_legend(nrow = 1, reverse = TRUE))


# MAKE PLOT PANELS FOR CNT pbt ####
plot_CNT_New1 <- df_joint_plot %>%
    dplyr::filter(pbt %>% str_detect("CNT_New1")) %>%
    arrange(pbt) %>%
    pull(plot) %>%
    wrap_plots(nrow = 1) +
    plot_layout(guides = "collect") &
    theme(legend.position = "top") &
    guides(fill = guide_legend(nrow = 1, reverse = TRUE))

plot_CNT_New2 <- df_joint_plot %>%
    dplyr::filter(pbt %>% str_detect("CNT_New2")) %>%
    arrange(pbt) %>%
    pull(plot) %>%
    wrap_plots(nrow = 1) +
    plot_layout(guides = "collect") &
    theme(legend.position = "top") &
    guides(fill = guide_legend(nrow = 1, reverse = TRUE))

plot_CNT_New3 <- df_joint_plot %>%
    dplyr::filter(pbt %>% str_detect("CNT_New3")) %>%
    arrange(pbt) %>%
    pull(plot) %>%
    wrap_plots(nrow = 1) +
    plot_layout(guides = "collect") &
    theme(legend.position = "top") &
    guides(fill = guide_legend(nrow = 1, reverse = TRUE))


# MAKE PLOT PANELS FOR CD_IC pbt ####
plot_CD_IC_New1 <- df_joint_plot %>%
    dplyr::filter(pbt %>% str_detect("CD_IC_New1")) %>%
    arrange(pbt) %>%
    pull(plot) %>%
    wrap_plots(nrow = 1) +
    plot_layout(guides = "collect") &
    theme(legend.position = "top") &
    guides(fill = guide_legend(nrow = 1, reverse = TRUE))

plot_CD_IC_New2 <- df_joint_plot %>%
    dplyr::filter(pbt %>% str_detect("CD_IC_New2")) %>%
    arrange(pbt) %>%
    pull(plot) %>%
    wrap_plots(nrow = 1) +
    plot_layout(guides = "collect") &
    theme(legend.position = "top") &
    guides(fill = guide_legend(nrow = 1, reverse = TRUE))


# MAKE PLOT PANELS FOR CD_PC pbt ####
plot_CD_PC_New1 <- df_joint_plot %>%
    dplyr::filter(pbt %>% str_detect("CD_PC_New1")) %>%
    arrange(pbt) %>%
    pull(plot) %>%
    wrap_plots(nrow = 1) +
    plot_layout(guides = "collect") &
    theme(legend.position = "top") &
    guides(fill = guide_legend(nrow = 1, reverse = TRUE))

plot_CD_PC_New2 <- df_joint_plot %>%
    dplyr::filter(pbt %>% str_detect("CD_PC_New2")) %>%
    arrange(pbt) %>%
    pull(plot) %>%
    wrap_plots(nrow = 1) +
    plot_layout(guides = "collect") &
    theme(legend.position = "top") &
    guides(fill = guide_legend(nrow = 1, reverse = TRUE))


# MAKE PLOT PANELS FOR EC pbt ####
plot_EC_New1 <- df_joint_plot %>%
    dplyr::filter(pbt %>% str_detect("EC_New1")) %>%
    arrange(pbt) %>%
    pull(plot) %>%
    wrap_plots(nrow = 1) +
    plot_layout(guides = "collect") &
    theme(legend.position = "top") &
    guides(fill = guide_legend(nrow = 1, reverse = TRUE))



# SAVE THE PLOTS ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/output/"


# PT PLOTS ####
ggsave(
    plot_PT_New1,
    file = paste(saveDir, "PT_New1 PBT vs PC.png"),
    dpi = 300,
    height = 15,
    width = 30,
    units = "cm"
)

ggsave(
    plot_PT_New2,
    file = paste(saveDir, "PT_New2 PBT vs PC.png"),
    dpi = 300,
    height = 15,
    width = 30,
    units = "cm"
)

ggsave(
    plot_PT_New3,
    file = paste(saveDir, "PT_New3 PBT vs PC.png"),
    dpi = 300,
    height = 15,
    width = 30,
    units = "cm"
)

ggsave(
    plot_PT_New4,
    file = paste(saveDir, "PT_New4 PBT vs PC.png"),
    dpi = 300,
    height = 15,
    width = 30,
    units = "cm"
)


# TL PLOTS ####
ggsave(
    plot_TL_New1,
    file = paste(saveDir, "TL_New1 PBT vs PC.png"),
    dpi = 300,
    height = 15,
    width = 30,
    units = "cm"
)


# TAL PLOTS ####
ggsave(
    plot_TAL_New2,
    file = paste(saveDir, "TAL_New2 PBT vs PC.png"),
    dpi = 300,
    height = 15,
    width = 30,
    units = "cm"
)

ggsave(
    plot_TAL_New3,
    file = paste(saveDir, "TAL_New3 PBT vs PC.png"),
    dpi = 300,
    height = 15,
    width = 30,
    units = "cm"
)

ggsave(
    plot_TAL_New4,
    file = paste(saveDir, "TAL_New4 PBT vs PC.png"),
    dpi = 300,
    height = 15,
    width = 30,
    units = "cm"
)


# DCT PLOTS ####
ggsave(
    plot_DCT_New1,
    file = paste(saveDir, "DCT_New1 PBT vs PC.png"),
    dpi = 300,
    height = 15,
    width = 30,
    units = "cm"
)

ggsave(
    plot_DCT_New2,
    file = paste(saveDir, "DCT_New2 PBT vs PC.png"),
    dpi = 300,
    height = 15,
    width = 30,
    units = "cm"
)

ggsave(
    plot_DCT_New3,
    file = paste(saveDir, "DCT_New3 PBT vs PC.png"),
    dpi = 300,
    height = 15,
    width = 30,
    units = "cm"
)

ggsave(
    plot_DCT_New4,
    file = paste(saveDir, "DCT_New4 PBT vs PC.png"),
    dpi = 300,
    height = 15,
    width = 30,
    units = "cm"
)


# CNT PLOTS ####
ggsave(
    plot_CNT_New1,
    file = paste(saveDir, "CNT_New1 PBT vs PC.png"),
    dpi = 300,
    height = 15,
    width = 30,
    units = "cm"
)

ggsave(
    plot_CNT_New2,
    file = paste(saveDir, "CNT_New2 PBT vs PC.png"),
    dpi = 300,
    height = 15,
    width = 30,
    units = "cm"
)

ggsave(
    plot_CNT_New3,
    file = paste(saveDir, "CNT_New3 PBT vs PC.png"),
    dpi = 300,
    height = 15,
    width = 30,
    units = "cm"
)


# CD_IC PLOTS ####
ggsave(
    plot_CD_IC_New1,
    file = paste(saveDir, "CD_IC_New1 PBT vs PC.png"),
    dpi = 300,
    height = 15,
    width = 30,
    units = "cm"
)

ggsave(
    plot_CD_IC_New2,
    file = paste(saveDir, "CD_IC_New2 PBT vs PC.png"),
    dpi = 300,
    height = 15,
    width = 30,
    units = "cm"
)


# CD_PC PLOTS ####
ggsave(
    plot_CD_PC_New1,
    file = paste(saveDir, "CD_PC_New1 PBT vs PC.png"),
    dpi = 300,
    height = 15,
    width = 30,
    units = "cm"
)

ggsave(
    plot_CD_PC_New2,
    file = paste(saveDir, "CD_PC_New2 PBT vs PC.png"),
    dpi = 300,
    height = 15,
    width = 30,
    units = "cm"
)


# EC PLOTS ####
ggsave(
    plot_EC_New1,
    file = paste(saveDir, "EC_New1 PBT vs PC.png"),
    dpi = 300,
    height = 15,
    width = 30,
    units = "cm"
)
