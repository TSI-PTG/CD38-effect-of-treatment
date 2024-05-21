# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(ggbeeswarm) # install.packages("ggbeeswarm")
library(ggpubr) # install.packages("ggpubr")
library(patchwork) # install.packages("patchwork")
library(ggrepel) # install.packages("ggrepel")
library(gghalves) # pak::pak("erocoar/gghalves")
library(ggh4x) # install.packages("ggh4x")
library(broom) # install.packages("broom") #for tabular model object transformations
library(car) # install.packages("car") #for type-II manovas
library(rstatix) # install.packages("rstatix") #for testing manova assumptions
library(vegan) # install.packages("vegan") #for permanova
library(pairwiseAdonis) # library(devtools) #install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(flextable) # install.packages("flextable") #for table outputs
library(officer) # install.packages("officer")
library(emmeans) # install.packages("emmeans") #for post-hoc testing and CLD
library(multcomp) # install.packages("multcomp") #for for CLD
library(ARTool) # install.packages("ARTool") #for non-parametric anova (aligned-rank test) and post-hoc
library(rcompanion) # install.packages("rcompanion") #for non-parametric anova (aligned-rank test) and post-hoc
# Custom operators, functions, and datasets
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
log10zero <- scales::trans_new(
    name = "log10zero",
    transform = function(x) log10(x + 1),
    inverse = function(x) 10^x - 1
)
# Suppress pesky dplyr reframe info
options(dplyr.reframe.inform = FALSE)
# source plot function
source("C:/R/CD38-effect-of-treatment/code/functions/plot.gg_violin_interaction.r")
# laod reference set
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0Georg Felz CD38 Vienna/G_Rstuff/data/vienna_5086_6Mar24.RData")


# DEFINE SEED ####
seed <- 42

# data_RF3min %>% print(n = "all")


# DEFINE THE SET ####
set <- vienna_5086[, vienna_5086$STUDY_EVALUATION_ID %nin% c(15, 18)]
# set <- vienna_5086[, vienna_5086$STUDY_EVALUATION_ID != 15 & vienna_5086$CEL != "FBN003_NBN010_B2_(PrimeView).CEL"]


set  %>% pData  %>% colnames


# WRANGLE THE PHENOTYPE DATA ####
df00 <- set %>%
    pData() %>%
    tibble() %>%
    dplyr::rename(
        Patient = STUDY_EVALUATION_ID,
        Felz = Felzartamab_presumed
    ) %>%
    mutate(
        Patient = Patient %>% factor(),
        Group = Group %>% factor(levels = c("Index", "FU1", "FU2")),
        Felz = Felz %>% factor(labels = c("NoFelz", "Felz")),
        TxBx = TxBx %>% as.numeric(),
        Group_Felz = paste(Group, Felz, sep = ":") %>%
            factor(levels = c("Index_NoFelz", "FU1_NoFelz", "FU2_NoFelz", "Index_Felz", "FU1_Felz", "FU2_Felz"))
    ) %>%
    arrange(Patient, Group) %>%
    expand_grid(category = c("Survival"),  annotation = "Survival", score = "Survival", variable = "Survival") %>%
    nest(.by = c(category, annotation, variable, score)) %>%
    mutate(
        features = "RF3min"  %>% list,
        data = pmap(
            list(features, data),
            function(features, data) {
                data %>%
                    dplyr::select(
                        CEL, Patient, Center, Group, Felz, Group_Felz, TxBx,
                        all_of(features)
                    )
            }
        )
    )



# UNIVARIATE MEANS AND MEDIANS ####
df_univariate_01 <- df00 %>%
    mutate(
        medians = map(
            data,
            function(data) {
                medians <- data %>%
                    reframe(
                        median = RF3min %>% median(),
                        IQR = RF3min %>% IQR(),
                        .by = c(Group, Felz)
                    ) %>%
                    mutate(
                        Group_Felz = paste(Felz, Group, sep = ":") %>%
                            factor(levels = c(
                                "NoFelz:Index", "NoFelz:FU1", "NoFelz:FU2",
                                "Felz:Index", "Felz:FU1", "Felz:FU2"
                                # "NoFelz:Index", "NoFelz:FU0","NoFelz:FU1", "NoFelz:FU2",
                                # "Felz:Index", "Felz:FU0", "Felz:FU1", "Felz:FU2"
                            ))
                    ) %>%
                    arrange(Group_Felz)
            }
        ),
        medians_delta = map(
            medians,
            function(medians) {
                median <- medians %>%
                    reframe(
                        median_delta = combn(median, 2, diff) %>% as.numeric(),
                        IQR_delta = combn(IQR, 2, mean) %>% as.numeric(),
                        .by = c(Felz)
                    ) %>%
                    mutate(Group_pairwise = rep(c("Index - FU1", "Index - FU2", "FU1 - FU2"), 2)) %>%
                    mutate(
                        median_delta_delta = combn(median_delta, 2, diff) %>% as.numeric(),
                        IQR_delta_delta = combn(IQR_delta, 2, mean) %>% as.numeric(),
                        .by = c(Group_pairwise)
                    )
            }
        )
    )
df_univariate_01$data[[1]]
df_univariate_01$medians[[1]]
df_univariate_01$medians_delta[[1]]


# UNIVARIATE NONPARAMETRIC TESTS ####
df_univariate_02 <- df_univariate_01 %>%
    mutate(
        art = map(data, art, formula = RF3min ~ Group * Felz + (1 | Patient)),
        art_lm = map(
            art,
            function(art) {
                art %>%
                    artlm(term = "Group:Felz", response = "aligned") %>%
                    summary() %>%
                    pluck("coefficients") %>%
                    as.data.frame() %>%
                    rownames_to_column("Term_lm") %>%
                    dplyr::filter(Term_lm != "(Intercept)")
            }
        ),
        art_aov = map(art, anova, type = "II"),
        art_aov_tidy = map(art_aov, tidy) %>% suppressWarnings(),
        art_con_interaction_default = map(
            art,
            art.con,
            formula = "Group:Felz", adjust = "fdr", method = "pairwise", interaction = TRUE, response = "art"
        ),
        art_con_interaction_default_tidy = map(art_con_interaction_default, tidy),
        art_con_cld = map(
            art_con_interaction_default_tidy,
            function(art_con_interaction_default_tidy) {
                art_con_interaction_default_tidy %>%
                    as.data.frame() %>%
                    cldList(adj.p.value ~ Group:Felz, data = .) # %>%
                #         arrange(Group %>%
                #             factor(
                #                 levels = c(
                #                     "Index,NoFelz", "FU1,NoFelz",
                #                     "Index,Felz", "FU1,Felz"
                #                 )
                #             ))
            }
        )
    )

names(df_univariate_02$art_aov_tidy) <- df_univariate_02$variable %>% as.character()
names(df_univariate_02$art_con_interaction_default_tidy) <- df_univariate_02$variable %>% as.character()
names(df_univariate_02$art_con_cld) <- df_univariate_02$variable %>% as.character()

df_univariate_02$art_con_cld

df_univariate_02$medians[[1]]
df_univariate_02$medians_delta[[1]]
df_univariate_02$art_aov_tidy[[1]]
df_univariate_02$art_con_interaction_default[[1]]

df_univariate_02$art_con_interaction_default_tidy[[1]]
df_univariate_02$art_con_cld[[1]]


# CREATE FLEXTABLE OF ART MODELS ####
title_art <- paste("Table i. Non-parametric ANOVA (ART) of molecular scores in biopsies from treated vs untreated patients")
title_art_pairwise <- paste("Table i. Pairwise Non-parametric ANOVA (ART) of molecular scores in biopsies from treated vs untreated patients")

res_art_flextable <- df_univariate_02 %>%
    dplyr::select(score, art_aov) %>%
    unnest(everything()) %>%
    dplyr::rename(p.value = `Pr(>F)`) %>%
    mutate(FDR = p.value) %>%
    dplyr::select(score, Term, `F`, FDR) %>%
    flextable::flextable() %>%
    flextable::add_header_row(values = rep(title_art, ncol_keys(.))) %>%
    flextable::merge_h(part = "header") %>%
    flextable::merge_v(j = 1:2) %>%
    flextable::fontsize(size = 8, part = "all") %>%
    flextable::align(align = "center", part = "all") %>%
    flextable::bg(bg = "white", part = "all") %>%
    flextable::bg(i = ~ FDR < 0.05 & Term == "Group:Felz", j = 3:4, bg = "#fbff00") %>%
    flextable::colformat_double(j = 2:3, digits = 2) %>%
    flextable::colformat_double(j = 4:ncol_keys(.), digits = 3) %>%
    flextable::border_remove() %>%
    flextable::bold(part = "header") %>%
    flextable::padding(padding = 0, part = "all") %>%
    flextable::border(border = fp_border(), part = "all") %>%
    flextable::autofit()

# res_art_flextable %>% print(preview = "pptx")


# FORMAT TABLES OF CONTRASTS ####
data_pairwise_formatted_delta <- df_univariate_02 %>%
    dplyr::select(annotation, score, art_con_interaction_default_tidy, art_con_cld, medians_delta) %>%
    unnest(c(art_con_interaction_default_tidy, art_con_cld), names_repair = tidyr_legacy) %>%
    unnest(medians_delta, names_repair = tidyr_legacy) %>%
    dplyr::rename(p_interaction = adj.p.value) %>%
    mutate(
        FDR_interaction = p_interaction %>% p.adjust(method = "fdr"), .by = annotation
    ) %>%
    mutate_at(
        vars(contains("p_i"), contains("FDR")),
        ~ ifelse(
            . < 0.01,
            formatC(., digits = 0, format = "e"),
            formatC(., digits = 3, format = "f")
        )
    ) %>%
    mutate(
        delta = paste(
            format(round(median_delta, 2), nsmall = 1),
            "\u00B1",
            round(IQR_delta, 2),
            # Letter,
            sep = " "
        ),
        deltadelta = paste(
            format(round(median_delta_delta, 2), nsmall = 1),
            "\u00B1",
            round(IQR_delta_delta, 2),
            # Letter,
            sep = " "
        )
    ) %>%
    dplyr::select(
        annotation, score, Group_pairwise, Group_pairwise1, FDR_interaction,
        Felz, delta, deltadelta
    ) %>%
    distinct(Group_pairwise1, score, Felz, .keep_all = TRUE) %>%
    nest(.by = c(annotation, score, Group_pairwise1)) %>%
    mutate(
        data = map(
            data,
            function(data) {
                data %>%
                    pivot_wider(names_from = Felz, values_from = delta) %>%
                    relocate(deltadelta, .after = last_col())
            }
        )
    )

data_delta_formatted <- df_univariate_02 %>%
    dplyr::select(annotation, score, art_con_interaction_default_tidy, art_con_cld, medians_delta) %>%
    unnest(medians_delta, names_repair = tidyr_legacy) %>%
    mutate_at(
        vars(contains("p_i"), contains("FDR")),
        ~ ifelse(
            . < 0.01,
            formatC(., digits = 0, format = "e"),
            formatC(., digits = 3, format = "f")
        )
    ) %>%
    mutate(
        delta = paste(
            format(round(median_delta, 2), nsmall = 1),
            "\u00B1",
            round(IQR_delta, 2),
            sep = " "
        ),
        deltadelta = paste(
            format(round(median_delta_delta, 2), nsmall = 1),
            "\u00B1",
            round(IQR_delta_delta, 2),
            sep = " "
        )
    ) %>%
    dplyr::select(annotation, score, Group_pairwise, Felz, delta, deltadelta) %>%
    distinct(Group_pairwise, score, Felz, .keep_all = TRUE) %>%
    nest(
        .by = c(annotation, score, Group_pairwise),
        medians = everything()
    ) %>%
    mutate(
        medians = map(
            medians,
            function(medians) {
                medians %>%
                    pivot_wider(names_from = Felz, values_from = delta) %>%
                    relocate(deltadelta, .after = last_col())
            }
        )
    )

data_res_art_formatted <- df_univariate_02 %>%
    dplyr::select(annotation, score, art_con_interaction_default_tidy, art_con_cld, medians_delta) %>%
    unnest(c(art_con_interaction_default_tidy, art_con_cld), names_repair = tidyr_legacy) %>%
    # mutate(FDR_interaction = adj.p.value %>% p.adjust(method = "fdr"), .by = annotation) %>%
    mutate(FDR_interaction = adj.p.value) %>%
    mutate_at(
        vars(contains("p_i"), contains("FDR")),
        ~ ifelse(
            . < 0.01,
            formatC(., digits = 0, format = "e"),
            formatC(., digits = 3, format = "f")
        )
    ) %>%
    dplyr::select(annotation, score, Group_pairwise, FDR_interaction, Letter) %>%
    nest(
        .by = c(annotation, score, Group_pairwise),
        contrasts = everything()
    )

data_pairwise_formatted <- data_delta_formatted %>%
    left_join(data_res_art_formatted, by = c("annotation", "score", "Group_pairwise")) %>%
    mutate(
        data = pmap(
            list(medians, contrasts),
            function(medians, contrasts) {
                medians %>%
                    left_join(contrasts, by = c("annotation", "score", "Group_pairwise")) %>%
                    # mutate(deltadelta = deltadelta %>% paste(Letter)) %>%
                    dplyr::select(-annotation, -score, -Group_pairwise, -Letter)
            }
        )
    ) %>%
    dplyr::select(-medians, -contrasts) %>%
    pivot_wider(names_from = Group_pairwise, values_from = data) %>%
    relocate(`FU1 - FU2`, .before = `Index - FU2`)

data_pairwise_formatted %>%
    dplyr::slice(1) %>%
    pull(4)

data_pairwise_formatted %>%
    dplyr::slice(1) %>%
    unnest(everything(), names_repair = tidyr_legacy)


# FORMAT FLEXTABLE ####
# define sample sizes
# df_n <- df_00 %>%
#     group_by(Group_Felz) %>%
#     tally()

# N <- df_n %>% pull(n)
title_art_pairwise <- paste("Table i. Median \u00B1 IQR molecular scores in biopsies from treated vs untreated patients")

footnoteText <- c(
    paste(
        "Grey shading denotes ANOVA interactive effect FDR < 0.05\n",
        "Bold denotes FDR < 0.05\n",
        "FDR correction was carried out within each annotation grouping",
        sep = ""
    ) # ,
    # "Scores are the mean fold change in expression for all Probes within the set vs the mean expression for all Probes in the NoPGD control set"
)

header1 <- c(
    "Annotation", "Score",
    rep("FU1 - Index", 4),
    rep("FU2 - FU1", 4),
    rep("FU2 - Index", 4)
)

header2 <- c(
    "Annotation", "Score",
    rep(c("\u394 NoFelz\n(N=11)", "\u394 Felz\n(N=11)", "\u394\u394", "\u394\u394 FDR"), 1),
    rep(c("\u394 NoFelz\n(N=10)", "\u394 Felz\n(N=10)", "\u394\u394", "\u394\u394 FDR"), 2)
)

cellWidths <- c(4, 11, rep(c(3, 3, 3, 2), 3))

category_vec1 <- str_remove(data_pairwise_formatted$score, "Prob\\)")
category_vec2 <- str_extract(data_pairwise_formatted$score, "Prob")
category_vec3 <- str_replace(category_vec2, "Prob", ")")


# CREATE FELXTABLE OF PAIRWISE ART MODELS ####
flextable_pairwise <- data_pairwise_formatted %>%
    unnest(everything(), names_repair = tidyr_legacy) %>%
    flextable::flextable() %>%
    flextable::compose(
        j = "score",
        value = as_paragraph(category_vec1, as_sub(category_vec2), category_vec3)
    ) %>%
    flextable::delete_part("header") %>%
    flextable::add_header_row(top = TRUE, values = header2) %>%
    flextable::add_header_row(top = TRUE, values = header1) %>%
    flextable::add_header_row(top = TRUE, values = rep(title_art_pairwise, ncol_keys(.))) %>%
    flextable::add_footer_row(values = footnoteText[[1]], colwidths = ncol_keys(.)) %>%
    flextable::merge_v(j = 1:2) %>%
    flextable::merge_v(part = "header") %>%
    flextable::merge_h(part = "header") %>%
    flextable::border_remove() %>%
    flextable::border(part = "header", border = fp_border()) %>%
    flextable::border(part = "body", border = fp_border()) %>%
    flextable::border(part = "footer", border.left = fp_border(), border.right = fp_border()) %>%
    flextable::border(i = 1, part = "footer", border.bottom = fp_border()) %>%
    flextable::align(align = "center") %>%
    flextable::align(align = "center", part = "header") %>%
    # flextable::valign(i = 3, j = ncol_keys(.), valign = "bottom", part = "header") %>%
    flextable::font(fontname = "Arial", part = "all") %>%
    flextable::fontsize(size = 8, part = "all") %>%
    flextable::fontsize(size = 8, part = "footer") %>%
    flextable::fontsize(i = 1, size = 12, part = "header") %>%
    flextable::bold(part = "header") %>%
    flextable::bold(j = 1, part = "body") %>%
    flextable::bg(bg = "white", part = "all") %>%
    flextable::bg(i = ~ as.numeric(FDR_interaction) < 0.05, j = 2:ncol_keys(.), bg = "grey90", part = "body") %>%
    flextable::bold(i = ~ as.numeric(FDR_interaction) < 0.05, j = ncol_keys(.), part = "body") %>%
    flextable::padding(padding = 0, part = "all") %>%
    flextable::width(width = cellWidths, unit = "cm") %>%
    flextable::width(., width = dim(.)$widths * 30 / (flextable_dim(.)$widths), unit = "cm")


# PRINT THE FLEXTABLES ####
# flextable_pairwise %>% print(preview = "pptx")


# PROCESS DATA FOR PLOTTING ####
df_plot <- df_univariate_02 %>%
    mutate(
        data_plot = pmap(
            list(data, variable, score, medians, medians_delta, art_con_interaction_default_tidy),
            function(data, variable, score, medians, medians_delta, art_con_interaction_default_tidy) {
                delta <- medians_delta %>%
                    mutate(
                        Felz = Felz %>% factor(labels = c("No Felzartamab", "Felzartamab"))
                    )
                delta_delta <- medians_delta %>%
                    distinct(Group_pairwise, .keep_all = TRUE) %>%
                    mutate(Group_pairwise = c("FU1 - Index", "FU2 - Index", "FU2 - FU1"))
                delta_delta_p <- art_con_interaction_default_tidy %>%
                    dplyr::select(Group_pairwise, adj.p.value) %>%
                    mutate(FDR = ifelse(
                        adj.p.value < 0.01,
                        formatC(adj.p.value, digits = 0, format = "e"),
                        formatC(adj.p.value, digits = 3, format = "f")
                    ))
                data <- data %>%
                    mutate(
                        variable = variable,
                        Felz = Felz %>% factor(labels = c("No Felzartamab", "Felzartamab"))
                    )
                data <- data %>%
                    left_join(
                        data %>%
                            reframe(median = median(RF3min), .by = c(Group, Felz, Patient)) %>%
                            spread(Group, median) %>%
                            mutate(delta = FU1 - Index, delta2 = FU2 - FU1) %>%
                            pivot_longer(cols = c(Index, FU1, FU2), names_to = "Group", values_to = "value") %>%
                            dplyr::select(Felz, Patient, Group, delta, delta2),
                        by = c("Felz", "Patient", "Group")
                    ) %>%
                    mutate(Group = Group %>% factor(levels = c("Index", "FU1", "FU2")))
                ymax <- data$RF3min %>% max()
                tibble(
                    data = data %>% list(),
                    delta_delta = delta_delta %>% list(),
                    delta_delta_p = delta_delta_p %>% list(),
                    ymax = ymax %>% list(),
                )
            }
        )
    )

df_plot$data_plot[[1]]$data


# PLOTTING GLOBALS ####
log_ticks <- c(
    seq(0, 1, length.out = 11),
    seq(1, 10, length.out = 10),
    seq(10, 100, length.out = 10),
    seq(100, 1000, length.out = 10)
)
labels <- c(0, 1, 10, 100, 1000)
dodge <- 0.3


# MAKE BIOPSY PAIR PLOTS ####
df_plot$data_plot[[1]]$data[[1]] %>%
    dplyr::filter(Patient == 9)


df_plot$data_plot[[1]]$data[[1]] %>%
    dplyr::select(Patient, Group, RF3min, delta, delta2) %>%
    mutate(
                        delta_prop = case_when(
                            Group == "Index" ~ lead(RF3min) / RF3min,
                            Group == "FU1" ~ RF3min / lag(RF3min),
                            Group == "FU2" ~ RF3min / lag(RF3min),
                        )%>% log2(),
                        delta_prop = case_when(
                            delta_prop == -Inf ~ min(delta_prop[delta_prop != -Inf]),
                            delta_prop == Inf ~ max(delta_prop[delta_prop != Inf]),
                            TRUE ~ delta_prop
                        )
    ) %>%
    print(n = "all")




c(min(df_plot$data_plot[[1]]$data[[1]]$delta), max(df_plot$data_plot[[1]]$data[[1]]$delta))


plot_violin_pairs <- df_plot %>%
    mutate(
        plot_violin = pmap(
            list(data_plot),
            function(data_plot) {
                data <- data_plot$data[[1]] %>%
                    mutate(
                        delta_prop = case_when(
                            Group == "Index" ~ lead(RF3min) / RF3min,
                            Group == "FU1" ~ RF3min / lag(RF3min),
                            Group == "FU2" ~ RF3min / lag(RF3min),
                        )%>% log2(),
                        delta_prop = case_when(
                            delta_prop == -Inf ~ min(delta_prop[delta_prop != -Inf]),
                            delta_prop == Inf ~ max(delta_prop[delta_prop != Inf]),
                            TRUE ~ delta_prop
                        )
                    )
                delta_delta <- data_plot$delta_delta[[1]]
                delta_delta_p <- data_plot$delta_delta_p[[1]]
                ymax <- data_plot$ymax[[1]]
                data %>%
                    ggplot(aes(x = Group, y = RF3min)) +
                    geom_half_violin(
                        inherit.aes = FALSE,
                        data = data %>% dplyr::filter(Group %in% c("Index")),
                        mapping = aes(
                            x = Group,
                            y = RF3min
                        ),
                        # bw = 0.05,
                        side = c("l"),
                        fill = "grey95",
                        col = "grey60",
                        trim = FALSE,
                        scale = "width"
                    ) +
                    geom_violin(
                        inherit.aes = FALSE,
                        data = data %>% dplyr::filter(Group == c("FU1")),
                        mapping = aes(
                            x = Group,
                            y = RF3min
                        ),
                        # bw =0.05,
                        fill = "#f2f2f2ec",
                        col = "#cdcdcdb8",
                        trim = FALSE,
                        scale = "width"
                    ) +
                    geom_half_violin(
                        inherit.aes = FALSE,
                        data = data %>% dplyr::filter(Group %in% c("FU2")),
                        mapping = aes(
                            x = Group,
                            y = RF3min
                        ),
                        # bw =  0.05),
                        side = c("r"),
                        fill = "grey95",
                        col = "grey60",
                        trim = FALSE,
                        scale = "width"
                    ) +
                    geom_point(
                        inherit.aes = FALSE,
                        data = data %>% dplyr::filter(Group %in% c("Index")),
                        mapping = aes(
                            x = Group,
                            y = RF3min,
                            col = delta_prop
                        ),
                        position = position_nudge(x = 0.1),
                        size = 2, alpha = 0.75
                    ) +
                    geom_point(
                        inherit.aes = FALSE,
                        data = data %>% dplyr::filter(Group %in% c("FU1")),
                        mapping = aes(
                            x = Group,
                            y = RF3min,
                            col = delta_prop
                        ),
                        size = 2, alpha = 0.75
                    ) +
                    geom_point(
                        inherit.aes = FALSE,
                        data = data %>% dplyr::filter(Group %in% c("FU2")),
                        mapping = aes(
                            x = Group,
                            y = RF3min,
                            col = delta_prop
                        ),
                        position = position_nudge(x = -0.1),
                        size = 2, alpha = 0.75
                    ) +
                    geom_line(
                        data = data %>% dplyr::filter(Group %in% c("Index", "FU1")),
                        mapping = aes(
                            x = Group,
                            col = delta_prop,
                            group = Patient
                        ),
                        position = position_nudge(
                            x = ifelse(data %>%
                                dplyr::filter(Group %in% c("Index", "FU1")) %>%
                                pull(Group) == "Index", 0.1, 0)
                        ),
                        linewidth = 0.5, alpha = 0.25,
                        show.legend = FALSE
                    ) +
                    geom_line(
                        data = data %>% dplyr::filter(Group %in% c("FU1", "FU2")),
                        mapping = aes(
                            x = Group,
                            col = delta_prop  %>% lead,
                            group = Patient
                        ),
                        position = position_nudge(
                            x = ifelse(data %>%
                                dplyr::filter(Group %in% c("FU1", "FU2")) %>%
                                pull(Group) == "FU1", 0, -0.1)
                        ),
                        linewidth = 0.5, alpha = 0.25,
                        show.legend = FALSE
                    ) +
                    stat_summary(
                        fun = median, geom = "point", size = 4
                    ) +
                    stat_summary(
                        aes(group = Felz),
                        fun = median,
                        geom = "line",
                        linewidth = 1,
                        linetype = "dashed"
                    ) +
                    labs(
                        title = paste(
                            paste(
                                "\u394\u394", delta_delta$Group_pairwise[[1]], "=",
                                delta_delta$median_delta_delta[[1]] %>% round(2),
                                "| FDR =",
                                delta_delta_p$FDR[[1]]
                            ),
                            paste(
                                "\u394\u394", delta_delta$Group_pairwise[[3]], "=",
                                delta_delta$median_delta_delta[[3]] %>% round(2),
                                "| FDR =",
                                delta_delta_p$FDR[[3]]
                            ),
                            paste(
                                "\u394\u394", delta_delta$Group_pairwise[[2]], "=",
                                delta_delta$median_delta_delta[[2]] %>% round(2),
                                "| FDR =",
                                delta_delta_p$FDR[[2]]
                            ),
                            sep = "\n"
                        ),
                        x = NULL,
                        y = "Survival probability",
                        col = "Individual patient response     ",
                        parse = TRUE
                    ) +
                    scale_color_gradient2(
                        low = "red",
                        mid = "grey60",
                        high = "#00ff00bc",
                        midpoint = 0,
                        breaks = c(min(data$delta_prop), max(data$delta_prop)),
                        labels = c("worsened", "improved"),
                        guide = guide_colorbar(
                            title.position = "top",
                            barwidth = 20,
                            ticks = FALSE,
                            label.hjust = c(1.1, -0.1),
                            label.vjust = 8,
                            reverse = TRUE
                        ),
                    ) +
                    # scale_x_discrete(labels = c("Index\n(n = 11)", "FU1\n(n = 11)", "FU2\n(n = 11)")) +
                    scale_x_discrete(labels = c("Index", "FU1", "FU2")) +
                    # scale_y_continuous(
                    #     expand = c(0, 0),
                    #     breaks = c(0, 1, 10, 100, 1000),
                    #     labels = c(0, 1, 10, 100, 1000),
                    #     minor_breaks = log_ticks,
                    #     trans = log10zero
                    # ) +
                    coord_cartesian(
                        xlim = c(NA, NA),
                        ylim = c(0, 1)
                    ) +
                    theme_bw() +
                    theme(
                        axis.ticks.length.y = unit(0.25, "cm"),
                        axis.title = element_text(size = 15),
                        axis.text = element_text(size = 12, colour = "black"),
                        panel.grid = element_blank(),
                        strip.text = element_text(size = 12, colour = "black"),
                        legend.position = "top",
                        legend.box.spacing = unit(-10, "pt"),
                        legend.background = element_blank(),
                        legend.title = element_text(size = 15, hjust = 0.66, vjust = 1, face = "bold"),
                        legend.text = element_text(size = 15),
                        plot.title = element_text(
                            colour = "black",
                            hjust = 0,
                            size = 12,
                            face = case_when(
                                delta_delta_p$adj.p.value[[1]] < 0.05 ~ "bold.italic",
                                delta_delta_p$adj.p.value[[2]] < 0.05 ~ "bold.italic",
                                TRUE ~ "italic"
                            )
                        )
                    ) +
                    facet_wrap(~Felz)
            }
        )
    )

plot_ddRF3min <- plot_violin_pairs$plot_violin[[1]] %>% ggarrange(common.legend = TRUE)



# MAKE JOINT PATIENT PAIR PLOTS ####
plot_patient_pairs <- df_univariate_02 %>%
    mutate(
        gg_line = pmap(
            list(data, variable, score),
            function(data, variable, score) {
                data <- data %>%
                    mutate(
                        variable = variable,
                        Felz = Felz %>% factor(labels = c("No Felzartamab", "Felzartamab"))
                    )
                data %>%
                    ggplot(
                        aes(
                            x = Group,
                            y = RF3min,
                            col = Patient,
                            group = Patient
                        )
                    ) +
                    geom_point(size = 5, position = position_dodge(width = dodge)) +
                    geom_line(
                        linewidth = 0.75,
                        linetype = "dashed",
                        position = position_dodge(width = dodge),
                        show.legend = FALSE
                    ) +
                    geom_text(
                        aes(label = Patient),
                        size = 4,
                        col = "black",
                        show.legend = FALSE,
                        position = position_dodge(width = dodge)
                    ) +
                    labs(
                        x = NULL,
                        y = "dd-RF3min (cp/mL)",
                        col = "STUDY_EVALUATION_ID     ",
                        parse = TRUE
                    ) +
                    scale_y_continuous(
                        breaks = c(0, 1, 10, 100, 1000),
                        labels = c(0, 1, 10, 100, 1000),
                        minor_breaks = log_ticks,
                        trans = log10zero
                    ) +
                    coord_cartesian(xlim = c(1.3, 2.7), ylim = c(0, 1000)) +
                    theme_bw() +
                    theme(
                        axis.ticks.length.y = unit(0.25, "cm"),
                        axis.title = element_text(size = 15),
                        axis.text = element_text(size = 12, colour = "black"),
                        panel.grid = element_blank(),
                        strip.text = element_text(size = 12, colour = "black"),
                        legend.position = "none",
                        legend.title = element_text(size = 15, hjust = 0.66, vjust = 1, face = "bold"),
                        legend.text = element_text(size = 15),
                    ) +
                    facet_wrap(data$Felz) +
                    guides(
                        col = guide_legend(nrow = 1),
                        y = ggprism::guide_prism_minor()
                    )
            }
        )
    )



# SAVE THE PLOTS ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0Georg Felz CD38 Vienna/G_Rstuff/output/"
ggsave(
    filename = paste(saveDir, "Felzartamab RF3min.png"),
    plot = plot_ddRF3min,
    dpi = 600,
    width = 18,
    height = 14,
    units = "cm",
    bg = "white"
)
ggsave(
    filename = paste(saveDir, "Felzartamab ddRF3min patient pairs.png"),
    plot = plot_patient_pairs$gg_line[[1]],
    dpi = 600,
    width = 18,
    height = 12,
    units = "cm",
    bg = "white"
)



# END ####
