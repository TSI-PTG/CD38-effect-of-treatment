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
library(flextable) # install.packages("flextable")  
library(officer) # install.packages("officer")
library(emmeans) # install.packages("emmeans") #for post-hoc testing and CLD
library(multcomp) # install.packages("multcomp") #for for CLD
library(ARTool) # install.packages("ARTool") #for non-parametric anova (aligned-rank test) and post-hoc
library(rcompanion) # install.packages("rcompanion") #for non-parametric anova (aligned-rank test) and post-hoc
# Custom operators, functions, and datasets
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
log10zero <- scales::trans_new(
    name = "log10zero",
    transform = function(x) log10(x + 0.1),
    inverse = function(x) 10^x - 0.1
)
# Suppress pesky dplyr reframe info
options(dplyr.reframe.inform = FALSE)
# source plot function
source("C:/R/CD38-effect-of-treatment/code/functions/plot.gg_violin_interaction.r")
source("C:/R/CD38-effect-of-treatment/code/functions/plot.gg_patient_pairs_interaction.r")
# load data
# load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0Georg Felzartamab CD38 Vienna/G_Rstuff/data/data_cfdna.RData")
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/data_scores_k1208.RData")


# DEFINE SEED ####
seed <- 42


# DEFINE DATA ####
data_cfdna <- data_scores_k1208 %>%
    dplyr::select(Center, Patient, Felzartamab, Group, Followup, Felzartamab_Group, Felzartamab_Followup, cfDNA) %>%
    dplyr::filter(Group != "FU1b", Patient %nin% c(15, 18)) %>%
    drop_na(cfDNA) %>%
    left_join(., summarise(., sample_pairs = n(), .by = Felzartamab_Followup), by = "Felzartamab_Followup") %>%
    arrange(Felzartamab, Patient, Group)


# WRANGLE THE PHENOTYPE DATA ####
df00 <- data_cfdna %>%
    expand_grid(
        category = c("cfDNA"),
        score = "Donor-derived cell-free DNA (dd-cfDNA, cp/mL)",
        variable = "cfDNA"
    ) %>%
    dplyr::rename(value = cfDNA) %>%
    mutate(annotation = "cfDNA") %>%
    nest(.by = c(category, annotation, variable, score))


# UNIVARIATE MEANS AND MEDIANS ####
# TODO: find a way to determine median_deltas without using combn while also being nonspecific to new data
df_univariate_01 <- df00 %>%
    mutate(
        medians = map(
            data,
            function(data) {
                data %>%
                    reframe(
                        median = value %>% median(),
                        IQR = value %>% IQR(),
                        sample_pairs = n(),
                        .by = c(Followup, Felzartamab, Felzartamab_Followup)
                    ) %>%
                    arrange(Felzartamab_Followup)
            }
        ),
        medians_delta = map(
            medians,
            function(medians) {
                medians %>%
                    reframe(
                        median_delta = combn(median, 2, diff) %>% as.numeric(),
                        IQR_delta = combn(IQR, 2, mean) %>% as.numeric(),
                        sample_pairs = sample_pairs,
                        .by = c(Felzartamab)
                    ) %>%
                    mutate(
                        Followup_pairwise = rep(c("Day0 - Week24", "Day0 - Week52", "Week24 - Week52"), 2) %>%
                            factor(levels = c("Day0 - Week24", "Day0 - Week52", "Week24 - Week52")),
                        .before = sample_pairs
                    ) %>%
                    mutate(
                        median_delta_delta = combn(median_delta, 2, diff) %>% as.numeric(),
                        IQR_delta_delta = combn(IQR_delta, 2, mean) %>% as.numeric(),
                        .by = c(Followup_pairwise),
                        .before = sample_pairs
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
        art = map(data, art, formula = value ~ Followup * Felzartamab + (1 | Patient)),
        art_lm = map(
            art,
            function(art) {
                art %>%
                    artlm(term = "Followup:Felzartamab", response = "aligned") %>%
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
            formula = "Followup:Felzartamab", adjust = "fdr", method = "pairwise", interaction = TRUE, response = "art"
        ),
        art_con_interaction_default_tidy = map(art_con_interaction_default, tidy),
        art_con_cld = map(
            art_con_interaction_default_tidy,
            function(art_con_interaction_default_tidy) {
                art_con_interaction_default_tidy %>%
                    as.data.frame() %>%
                    cldList(adj.p.value ~ Followup:Felzartamab, data = .) 
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
    dplyr::select(annotation, score, art_aov) %>%
    unnest(everything()) %>%
    dplyr::rename(p.value = `Pr(>F)`) %>%
    mutate(FDR = p.value) %>%
    dplyr::select(annotation, score, Term, `F`, FDR) %>%
    flextable::flextable() %>%
    flextable::add_header_row(values = rep(title_art, ncol_keys(.))) %>%
    flextable::merge_h(part = "header") %>%
    flextable::merge_v(j = 1:2) %>%
    flextable::fontsize(size = 8, part = "all") %>%
    flextable::align(align = "center", part = "all") %>%
    flextable::bg(bg = "white", part = "all") %>%
    flextable::bg(i = ~ FDR < 0.05 & Term == "Group:Felzartamab", j = 3:4, bg = "#fbff00") %>%
    flextable::colformat_double(j = 2:3, digits = 2) %>%
    flextable::colformat_double(j = 4:ncol_keys(.), digits = 3) %>%
    flextable::border_remove() %>%
    flextable::bold(part = "header") %>%
    flextable::padding(padding = 0, part = "all") %>%
    flextable::border(border = officer::fp_border(), part = "all") %>%
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
        annotation, score, Followup_pairwise, Followup_pairwise1, FDR_interaction,
        Felzartamab, delta, deltadelta
    ) %>%
    distinct(score, Felzartamab, .keep_all = TRUE) %>%
    nest(.by = c(annotation, score)) %>%
    mutate(
        data = map(
            data,
            function(data) {
                data %>%
                    pivot_wider(names_from = Felzartamab, values_from = delta) %>%
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
    dplyr::select(annotation, score, Followup_pairwise, Felzartamab, delta, deltadelta) %>%
    distinct(Followup_pairwise, score, Felzartamab, .keep_all = TRUE) %>%
    nest(
        .by = c(annotation, score, Followup_pairwise),
        medians = everything()
    ) %>%
    mutate(
        medians = map(
            medians,
            function(medians) {
                medians %>%
                    pivot_wider(names_from = Felzartamab, values_from = delta) %>%
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
    dplyr::select(annotation, score, Followup_pairwise, FDR_interaction, Letter) %>%
    nest(
        .by = c(annotation, score, Followup_pairwise),
        contrasts = everything()
    )

data_pairwise_formatted <- data_delta_formatted %>%
    left_join(data_res_art_formatted, by = c("annotation", "score", "Followup_pairwise")) %>%
    mutate(
        data = pmap(
            list(medians, contrasts),
            function(medians, contrasts) {
                medians %>%
                    left_join(contrasts, by = c("annotation", "score", "Followup_pairwise")) %>%
                    dplyr::select(-annotation, -score, -Followup_pairwise, -Letter)
            }
        )
    ) %>%
    dplyr::select(-medians, -contrasts) %>%
    pivot_wider(names_from = Followup_pairwise, values_from = data) %>%
    relocate(`Week24 - Week52`, .before = `Day0 - Week52`)

# data_pairwise_formatted %>%
#     dplyr::slice(1) %>%
#     pull(4)

# data_pairwise_formatted %>%
#     dplyr::slice(1) %>%
#     unnest(everything(), names_repair = tidyr_legacy)


# UNIVERSAL VARIABLES FOR FLEXTABLE ####
# define sample sizes
# df_n <- df_00 %>%
#     group_by(Group_Felz) %>%
#     tally()

# N <- df_n %>% pull(n)
title_art_pairwise <- paste("Table i. Median \u00B1 IQR molecular scores in biopsies from placebo vs felzartamab patients")

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
    rep("Week24 - Day0", 4),
    rep("Week52 - Week24", 4),
    rep("Week52 - Day0", 4)
)

header2 <- c(
    "Annotation", "Score",
    rep(c("\u394 Placebo\n(N=10)", "\u394 Felzartamab\n(N=10)", "\u394\u394", "\u394\u394 FDR"), 1),
    rep(c("\u394 Placebo\n(N=10)", "\u394 Felzartamab\n(N=10)", "\u394\u394", "\u394\u394 FDR"), 2)
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
    flextable::border(part = "header", border = officer::fp_border()) %>%
    flextable::border(part = "body", border = officer::fp_border()) %>%
    flextable::border(part = "footer", border.left = officer::fp_border(), border.right = officer::fp_border()) %>%
    flextable::border(i = 1, part = "footer", border.bottom = officer::fp_border()) %>%
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



# MAKE BIOPSY PAIR PLOTS ####
plot_cfDNA <- df_univariate_02 %>%
    mutate(
        plot_violin = pmap(
            list(data, variable, score, medians_delta, art_con_interaction_default_tidy),
            gg_violin_interaction
        ),
        plot_patient_pairs = pmap(
            list(data, variable, score),
            gg_patient_pairs_interaction
        )
    )

plot_cfDNA$plot_violin[[1]]
plot_cfDNA$plot_patient_pairs[[1]]


# # SAVE THE PLOTS ####
# saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/output/"
# ggsave(
#     filename = paste(saveDir, "Felzartamab ddcfDNA.png"),
#     plot = plot_ddcfDNA,
#     dpi = 600,
#     width = 18,
#     height = 14,
#     units = "cm",
#     bg = "white"
# )
# ggsave(
#     filename = paste(saveDir, "Felzartamab ddcfDNA patient pairs.png"),
#     plot = plot_patient_pairs$gg_line[[1]],
#     dpi = 600,
#     width = 18,
#     height = 12,
#     units = "cm",
#     bg = "white"
# )


# SAVE DATA NEST ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/"
save(plot_cfDNA, file = paste(saveDir, "Felzartamab_cfDNA_results.RData", sep = ""))



# END ####
