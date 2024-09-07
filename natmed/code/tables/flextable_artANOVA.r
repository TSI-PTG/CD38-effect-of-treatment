# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(flextable) # install.packages("flextable")
library(officer) # install.packages("officer")
# Custom operators and functions
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
# load anova results
load("natmed/results/ARTanova.RData")


# CREATE FLEXTABLE OF ART MODELS ####
title_art <- paste("Table i. Non-parametric ANOVA (ART) of molecular scores in biopsies from treated vs untreated patients")
title_art_pairwise <- paste("Table i. Pairwise Non-parametric ANOVA (ART) of molecular scores in biopsies from treated vs untreated patients")

res_art_flextable <- artANOVA %>%
    dplyr::select(annotation, n_cat, score, art_aov) %>%
    unnest(everything()) %>%
    dplyr::rename(p.value = `Pr(>F)`) %>%
    mutate(FDR = p.value %>% p.adjust(method = "fdr"), .by = annotation) %>%
    dplyr::select(annotation, score, Term, `F`, p.value, FDR) %>%
    flextable::flextable() %>%
    flextable::add_header_row(values = rep(title_art, ncol_keys(.))) %>%
    flextable::merge_h(part = "header") %>%
    flextable::merge_v(j = 1:2) %>%
    flextable::fontsize(size = 8, part = "all") %>%
    flextable::align(align = "center", part = "all") %>%
    flextable::bg(bg = "white", part = "all") %>%
    flextable::bg(i = ~ FDR < 0.05 & Term == "Followup:Felzartamab", j = 3:6, bg = "#fbff00") %>%
    flextable::colformat_double(j = 2:3, digits = 2) %>%
    flextable::colformat_double(j = 4:ncol_keys(.), digits = 3) %>%
    flextable::border_remove() %>%
    flextable::bold(part = "header") %>%
    flextable::padding(padding = 0, part = "all") %>%
    flextable::border(border = officer::fp_border(), part = "all") %>%
    flextable::autofit()
# res_art_flextable %>% print(preview = "pptx")


# FORMAT TABLES OF CONTRASTS ####
# TODO: can this be simplified?
data_pairwise_formatted_delta <- artANOVA %>%
    dplyr::select(annotation, score, art_aov_contrast_tidy, art_aov_contrast_cld, medians_delta) %>%
    unnest(c(art_aov_contrast_tidy, art_aov_contrast_cld), names_repair = tidyr_legacy) %>%
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
            " (",
            round(IQR_delta, 2),
            ")",
            sep = ""
        ),
        deltadelta = paste(
            format(round(median_delta_delta, 2), nsmall = 1),
            " (",
            round(IQR_delta_delta, 2),
            ")",
            sep = ""
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

data_delta_formatted <- artANOVA %>%
    dplyr::select(annotation, score, art_aov_contrast_tidy, art_aov_contrast_cld, medians_delta) %>%
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
            " (",
            round(IQR_delta, 2),
            ")",
            sep = ""
        ),
        deltadelta = paste(
            format(round(median_delta_delta, 2), nsmall = 1),
            " (",
            round(IQR_delta_delta, 2),
            ")",
            sep = ""
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

data_res_art_formatted <- artANOVA %>%
    dplyr::select(annotation, score, art_aov_contrast_tidy, art_aov_contrast_cld, medians_delta) %>%
    unnest(c(art_aov_contrast_tidy, art_aov_contrast_cld), names_repair = tidyr_legacy) %>%
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
    relocate(`Week24 - Week52`, .before = `Baseline - Week52`)

# data_pairwise_formatted %>%
#     dplyr::slice(1) %>%
#     pull(4)

# data_pairwise_formatted %>%
#     dplyr::slice(1) %>%
#     unnest(everything(), names_repair = tidyr_legacy)


# UNIVERSAL VARIABLES FOR FLEXTABLE ####
title_art_pairwise <- paste("Table 1. Effect of felzartamab on molecular scores")

footnoteText <- c(
    paste(
        "Grey shading denotes ANOVA interactive effect FDR < 0.05\n",
        "FDR correction was carried out within each annotation grouping",
        sep = ""
    ) # ,
    # "Scores are the mean fold change in expression for all Probes within the data vs the mean expression for all Probes in the NoPGD control data"
)

header1 <- c(
    "Annotation", "Score",
    rep("Baseline - Week 24", 4),
    rep("Week 24 - Week 52", 4),
    rep("Baseline - Week 52", 4)
)

header2 <- c(
    "Annotation", "Score",
    rep(c("\u394 Placebo\n(N=10)", "\u394 Felzartamab\n(N=10)", "\u394\u394", "\u394\u394 FDR"), 1),
    rep(c("\u394 Placebo\n(N=10)", "\u394 Felzartamab\n(N=10)", "\u394\u394", "\u394\u394 FDR"), 2)
)

cellWidths <- c(4, 11, rep(c(3.5, 3.5, 3.5, 2), 3))

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
    flextable::font(fontname = "Arial", part = "all") %>%
    flextable::fontsize(size = 8, part = "all") %>%
    flextable::fontsize(size = 8, part = "footer") %>%
    flextable::fontsize(i = 1, size = 12, part = "header") %>%
    flextable::bold(part = "header") %>%
    flextable::bold(j = 1, part = "body") %>%
    flextable::bg(bg = "white", part = "all") %>%
    flextable::bg(i = ~ as.numeric(FDR_interaction) < 0.05, j = 2:ncol_keys(.), bg = "grey90", part = "body") %>%
    # flextable::bold(i = ~ as.numeric(FDR_interaction) < 0.05, j = ncol_keys(.), part = "body") %>%
    flextable::padding(padding = 0, part = "all") %>%
    flextable::width(width = cellWidths, unit = "cm") %>%
    flextable::width(., width = dim(.)$widths * 30 / (flextable_dim(.)$widths), unit = "cm")


# SAVE THE FLEXTABLE
saveDir <- "natmed/results/"
save(flextable_pairwise, file = paste0(saveDir, "flextable_artANOVA.RData"))


# PRINT THE FLEXTABLES ####
# flextable_pairwise %>% print(preview = "pptx")
