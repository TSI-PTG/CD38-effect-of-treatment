# HOUSEKEEPING ####
# cran packages
library(tidyverse) # install.packages('tidyverse')
library(patchwork) # install.packages('patchwork')
library(gghalves) # pak::pak("erocoar/gghalves")
# bioconductor packages
library(Biobase) # BiocManager::install('biobroom')
# LOAD FUNCTIONS AND CUSTOM OPERATORS
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
# GenPBTscores IS A FUNCTION TO CALCULATE MEAN LOG2 FOLD CHANGE IN EXPRESSION OF ALL GENES WITHIN PBTS FOR EACH BIOPSY VERSUS THE CONTROL SET
GenPBTscores <- function(set, ctrl, pbtlist) {
    for (i in 1:length(names(pbtlist))) {
        pData(set)[names(pbtlist)[i]] <- NA
        conavg <- apply(Biobase::exprs(ctrl)[pbtlist[[i]], ], 1, mean)
        smallset <- set[pbtlist[[i]], ]
        dog <- sweep(Biobase::exprs(smallset), 1, conavg)
        dog <- apply(dog, 2, mean)
        sel <- which(names(pData(set)) == names(pbtlist)[i])
        pData(set)[, sel] <- dog
    }
    return(set)
}
# load data
if (!exists("K5086")) {
    load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/K5086_29jul2024.RData")
}
# load nephrectomies
load("Z:/Genome-Archive/RefData/KidneyReports2021REDCap/data/Conset.RData") # 4 Controls used for PBT calculations
# load new PBT list
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/PBTlist219_AKIinduced_refsplit.RData")
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/PBTlist219_AKIinduced_subnetwork.RData")
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/PBTlist219_AKIinduced_PCsplit.RData")
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/PBTlist219_AKIinduced_26Jul2024_PTG.RData")


# WRANGLE REFERENCE SET ####
K4502 <- K5086[, which(K5086$Cort > 0.1)]
pData(K4502) <- K4502 %>%
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
    dplyr::select(injAA)



# RENAME PBTS FOR PLOTTING ####
names(PBTlist219_AKIinduced_refsplit) <- PBTlist219_AKIinduced_refsplit %>%
    names() %>%
    str_replace_all("up", "control_up") %>%
    str_replace_all("dn", "control_dn")

names(PBTlist219_AKIinduced_PCsplit) <- PBTlist219_AKIinduced_PCsplit %>%
    names() %>%
    str_replace_all("up", "PC_up") %>%
    str_replace_all("dn", "PC_dn")

names(PBTlist219_AKIinduced_subnetwork) <- PBTlist219_AKIinduced_subnetwork %>%
    names() %>%
    str_replace_all("blue", "net_blue") %>%
    str_replace_all("turquoise", "net_turquoise") %>%
    str_replace_all("brown", "net_brown") %>%
    str_replace_all("yellow", "net_yellow") %>%
    str_replace_all("green", "net_green")


# ADD PBT SCORES TO REFERENCE SET ####
K4502 <- GenPBTscores(
    K4502,
    Conset,
    c(
        PBTlist219_AKIinduced,
        PBTlist219_AKIinduced_refsplit,
        PBTlist219_AKIinduced_PCsplit,
        PBTlist219_AKIinduced_subnetwork
    )
)

PBT_names <- K4502 %>%
    pData() %>%
    dplyr::select(-injAA) %>%
    colnames()


# WRANGLE THE DATA FOR PLOTTING ####
df_injAA <- K4502 %>%
    pData() %>%
    dplyr::select(injAA, all_of(PBT_names)) %>%
    pivot_longer(-c(injAA), names_to = "PBT", values_to = "score") %>%
    dplyr::rename(group = "injAA") %>%
    mutate(
        PBT = PBT %>% factor(
            levels = PBT_names,
            labels = c(PBT_names)
        ),
        group = group %>% factor(levels = c("MildCKD", "CKDAKI", "AKI1", "AKI2", "Normal") %>% rev())
    ) %>%
    nest(.by = PBT)



# PLOT PBTS ####
df_joint_plot <- df_injAA %>%
    mutate(
        plot = map2(
            PBT, data,
            function(PBT, data) {
                dotsize <- case_when(
                    PBT %>% str_detect("_net_") ~ 1,
                    PBT %>% str_detect("_control_") ~ 1,
                    PBT %>% str_detect("_PC_up|_PC_dn") ~ 1,
                    TRUE ~ 1
                )
                method <- case_when(
                    PBT %>% str_detect("_control_") ~ "Method 3\ncontrol method",
                    PBT %>% str_detect("_PC_up|_PC_dn") ~ "Method 4\nPC method",
                    PBT %>% str_detect("_net_") ~ "Method 5\nSigned sub-networks",
                    TRUE ~ "Method 1\ndo nothing"
                )
                mean_normal <- data %>%
                    dplyr::filter(group == "Normal") %>%
                    pull(score) %>%
                    mean()
                data %>%
                    mutate(
                        method = method %>% factor(levels = c(
                            "Method 1\ndo nothing",
                            "Method 3\ncontrol method",
                            "Method 4\nPC method",
                            "Method 5\nSigned sub-networks"
                        )),
                        PBT = PBT
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
                    geom_half_dotplot(
                        alpha = 0.85,
                        dotsize = dotsize,
                        stroke = 0.25,
                        stackratio = 1,
                        binwidth = 0.01,
                        stackdir = "down",
                        show.legend = FALSE,
                        col = "grey60"
                    ) +
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
                        y = "PBT score",
                        x = NULL
                    ) +
                    theme_bw() +
                    facet_wrap(~ method + PBT) +
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

df_joint_plot$plot[[5]]


# MAKE PLOT PANELS FOR PT PBT ####
plot_PT_New1 <- df_joint_plot %>%
    dplyr::filter(PBT %>% str_detect("PT_New1")) %>%
    arrange(PBT) %>%
    pull(plot) %>%
    wrap_plots(nrow = 1) +
    plot_layout(guides = "collect") &
    theme(legend.position = "top") &
    guides(fill = guide_legend(nrow = 1, reverse = TRUE))

plot_PT_New2 <- df_joint_plot %>%
    dplyr::filter(PBT %>% str_detect("PT_New2")) %>%
    arrange(PBT) %>%
    pull(plot) %>%
    wrap_plots(nrow = 1) +
    plot_layout(guides = "collect") &
    theme(legend.position = "top") &
    guides(fill = guide_legend(nrow = 1, reverse = TRUE))

plot_PT_New3 <- df_joint_plot %>%
    dplyr::filter(PBT %>% str_detect("PT_New3")) %>%
    arrange(PBT) %>%
    pull(plot) %>%
    wrap_plots(nrow = 1) +
    plot_layout(guides = "collect") &
    theme(legend.position = "top") &
    guides(fill = guide_legend(nrow = 1, reverse = TRUE))

plot_PT_New4 <- df_joint_plot %>%
    dplyr::filter(PBT %>% str_detect("PT_New4")) %>%
    arrange(PBT) %>%
    pull(plot) %>%
    wrap_plots(nrow = 1) +
    plot_layout(guides = "collect") &
    theme(legend.position = "top") &
    guides(fill = guide_legend(nrow = 1, reverse = TRUE))


# MAKE PLOT PANELS FOR TL PBT ####
plot_TL_New1 <- df_joint_plot %>%
    dplyr::filter(PBT %>% str_detect("TL_New1")) %>%
    arrange(PBT) %>%
    pull(plot) %>%
    wrap_plots(nrow = 1) +
    plot_layout(guides = "collect") &
    theme(legend.position = "top") &
    guides(fill = guide_legend(nrow = 1, reverse = TRUE))


# MAKE PLOT PANELS FOR TAL PBT ####
plot_TAL_New2 <- df_joint_plot %>%
    dplyr::filter(PBT %>% str_detect("TAL_New2")) %>%
    arrange(PBT) %>%
    pull(plot) %>%
    wrap_plots(nrow = 1) +
    plot_layout(guides = "collect") &
    theme(legend.position = "top") &
    guides(fill = guide_legend(nrow = 1, reverse = TRUE))

plot_TAL_New3 <- df_joint_plot %>%
    dplyr::filter(PBT %>% str_detect("TAL_New3")) %>%
    arrange(PBT) %>%
    pull(plot) %>%
    wrap_plots(nrow = 1) +
    plot_layout(guides = "collect") &
    theme(legend.position = "top") &
    guides(fill = guide_legend(nrow = 1, reverse = TRUE))

plot_TAL_New4 <- df_joint_plot %>%
    dplyr::filter(PBT %>% str_detect("TAL_New4")) %>%
    arrange(PBT) %>%
    pull(plot) %>%
    wrap_plots(nrow = 1) +
    plot_layout(guides = "collect") &
    theme(legend.position = "top") &
    guides(fill = guide_legend(nrow = 1, reverse = TRUE))


# MAKE PLOT PANELS FOR DCT PBT ####
plot_DCT_New1 <- df_joint_plot %>%
    dplyr::filter(PBT %>% str_detect("DCT_New1")) %>%
    arrange(PBT) %>%
    pull(plot) %>%
    wrap_plots(nrow = 1) +
    plot_layout(guides = "collect") &
    theme(legend.position = "top") &
    guides(fill = guide_legend(nrow = 1, reverse = TRUE))

plot_DCT_New2 <- df_joint_plot %>%
    dplyr::filter(PBT %>% str_detect("DCT_New2")) %>%
    arrange(PBT) %>%
    pull(plot) %>%
    wrap_plots(nrow = 1) +
    plot_layout(guides = "collect") &
    theme(legend.position = "top") &
    guides(fill = guide_legend(nrow = 1, reverse = TRUE))

plot_DCT_New3 <- df_joint_plot %>%
    dplyr::filter(PBT %>% str_detect("DCT_New3")) %>%
    arrange(PBT) %>%
    pull(plot) %>%
    wrap_plots(nrow = 1) +
    plot_layout(guides = "collect") &
    theme(legend.position = "top") &
    guides(fill = guide_legend(nrow = 1, reverse = TRUE))

plot_DCT_New4 <- df_joint_plot %>%
    dplyr::filter(PBT %>% str_detect("DCT_New4")) %>%
    arrange(PBT) %>%
    pull(plot) %>%
    wrap_plots(nrow = 1) +
    plot_layout(guides = "collect") &
    theme(legend.position = "top") &
    guides(fill = guide_legend(nrow = 1, reverse = TRUE))


# MAKE PLOT PANELS FOR CNT PBT ####
plot_CNT_New1 <- df_joint_plot %>%
    dplyr::filter(PBT %>% str_detect("CNT_New1")) %>%
    arrange(PBT) %>%
    pull(plot) %>%
    wrap_plots(nrow = 1) +
    plot_layout(guides = "collect") &
    theme(legend.position = "top") &
    guides(fill = guide_legend(nrow = 1, reverse = TRUE))

plot_CNT_New2 <- df_joint_plot %>%
    dplyr::filter(PBT %>% str_detect("CNT_New2")) %>%
    arrange(PBT) %>%
    pull(plot) %>%
    wrap_plots(nrow = 1) +
    plot_layout(guides = "collect") &
    theme(legend.position = "top") &
    guides(fill = guide_legend(nrow = 1, reverse = TRUE))

plot_CNT_New3 <- df_joint_plot %>%
    dplyr::filter(PBT %>% str_detect("CNT_New3")) %>%
    arrange(PBT) %>%
    pull(plot) %>%
    wrap_plots(nrow = 1) +
    plot_layout(guides = "collect") &
    theme(legend.position = "top") &
    guides(fill = guide_legend(nrow = 1, reverse = TRUE))


# MAKE PLOT PANELS FOR CD_IC PBT ####
plot_CD_IC_New1 <- df_joint_plot %>%
    dplyr::filter(PBT %>% str_detect("CD_IC_New1")) %>%
    arrange(PBT) %>%
    pull(plot) %>%
    wrap_plots(nrow = 1) +
    plot_layout(guides = "collect") &
    theme(legend.position = "top") &
    guides(fill = guide_legend(nrow = 1, reverse = TRUE))

plot_CD_IC_New2 <- df_joint_plot %>%
    dplyr::filter(PBT %>% str_detect("CD_IC_New2")) %>%
    arrange(PBT) %>%
    pull(plot) %>%
    wrap_plots(nrow = 1) +
    plot_layout(guides = "collect") &
    theme(legend.position = "top") &
    guides(fill = guide_legend(nrow = 1, reverse = TRUE))


# MAKE PLOT PANELS FOR CD_PC PBT ####
plot_CD_PC_New1 <- df_joint_plot %>%
    dplyr::filter(PBT %>% str_detect("CD_PC_New1")) %>%
    arrange(PBT) %>%
    pull(plot) %>%
    wrap_plots(nrow = 1) +
    plot_layout(guides = "collect") &
    theme(legend.position = "top") &
    guides(fill = guide_legend(nrow = 1, reverse = TRUE))

plot_CD_PC_New2 <- df_joint_plot %>%
    dplyr::filter(PBT %>% str_detect("CD_PC_New2")) %>%
    arrange(PBT) %>%
    pull(plot) %>%
    wrap_plots(nrow = 1) +
    plot_layout(guides = "collect") &
    theme(legend.position = "top") &
    guides(fill = guide_legend(nrow = 1, reverse = TRUE))


# MAKE PLOT PANELS FOR EC PBT ####
plot_EC_New1 <- df_joint_plot %>%
    dplyr::filter(PBT %>% str_detect("EC_New1")) %>%
    arrange(PBT) %>%
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
    file = paste(saveDir, "PT_New1 purification.png"),
    dpi = 300,
    height = 20,
    width = 40,
    units = "cm"
)

ggsave(
    plot_PT_New2,
    file = paste(saveDir, "PT_New2 purification.png"),
    dpi = 300,
    height = 20,
    width = 40,
    units = "cm"
)

ggsave(
    plot_PT_New3,
    file = paste(saveDir, "PT_New3 purification.png"),
    dpi = 300,
    height = 20,
    width = 40,
    units = "cm"
)

ggsave(
    plot_PT_New4,
    file = paste(saveDir, "PT_New4 purification.png"),
    dpi = 300,
    height = 20,
    width = 40,
    units = "cm"
)


# TL PLOTS ####
ggsave(
    plot_TL_New1,
    file = paste(saveDir, "TL_New1 purification.png"),
    dpi = 300,
    height = 20,
    width = 40,
    units = "cm"
)


# TAL PLOTS ####
ggsave(
    plot_TAL_New2,
    file = paste(saveDir, "TAL_New2 purification.png"),
    dpi = 300,
    height = 20,
    width = 40,
    units = "cm"
)

ggsave(
    plot_TAL_New3,
    file = paste(saveDir, "TAL_New3 purification.png"),
    dpi = 300,
    height = 20,
    width = 40,
    units = "cm"
)

ggsave(
    plot_TAL_New4,
    file = paste(saveDir, "TAL_New4 purification.png"),
    dpi = 300,
    height = 20,
    width = 40,
    units = "cm"
)


# DCT PLOTS ####
ggsave(
    plot_DCT_New1,
    file = paste(saveDir, "DCT_New1 purification.png"),
    dpi = 300,
    height = 20,
    width = 40,
    units = "cm"
)

ggsave(
    plot_DCT_New2,
    file = paste(saveDir, "DCT_New2 purification.png"),
    dpi = 300,
    height = 20,
    width = 40,
    units = "cm"
)

ggsave(
    plot_DCT_New3,
    file = paste(saveDir, "DCT_New3 purification.png"),
    dpi = 300,
    height = 20,
    width = 60,
    units = "cm"
)

ggsave(
    plot_DCT_New4,
    file = paste(saveDir, "DCT_New4 purification.png"),
    dpi = 300,
    height = 20,
    width = 40,
    units = "cm"
)


# CNT PLOTS ####
ggsave(
    plot_CNT_New1,
    file = paste(saveDir, "CNT_New1 purification.png"),
    dpi = 300,
    height = 20,
    width = 60,
    units = "cm"
)

ggsave(
    plot_CNT_New2,
    file = paste(saveDir, "CNT_New2 purification.png"),
    dpi = 300,
    height = 20,
    width = 50,
    units = "cm"
)

ggsave(
    plot_CNT_New3,
    file = paste(saveDir, "CNT_New3 purification.png"),
    dpi = 300,
    height = 20,
    width = 40,
    units = "cm"
)


# CD_IC PLOTS ####
ggsave(
    plot_CD_IC_New1,
    file = paste(saveDir, "CD_IC_New1 purification.png"),
    dpi = 300,
    height = 20,
    width = 50,
    units = "cm"
)

ggsave(
    plot_CD_IC_New2,
    file = paste(saveDir, "CD_IC_New2 purification.png"),
    dpi = 300,
    height = 20,
    width = 40,
    units = "cm"
)


# CD_PC PLOTS ####
ggsave(
    plot_CD_PC_New1,
    file = paste(saveDir, "CD_PC_New1 purification.png"),
    dpi = 300,
    height = 20,
    width = 40,
    units = "cm"
)

ggsave(
    plot_CD_PC_New2,
    file = paste(saveDir, "CD_PC_New2 purification.png"),
    dpi = 300,
    height = 20,
    width = 60,
    units = "cm"
)


# EC PLOTS ####
ggsave(
    plot_EC_New1,
    file = paste(saveDir, "EC_New1 purification.png"),
    dpi = 300,
    height = 20,
    width = 40,
    units = "cm"
)

