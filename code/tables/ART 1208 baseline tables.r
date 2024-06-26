# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(flextable) # install.packages("flextable") #for table outputs
library(officer) # install.packages("officer")
# Custom operators, functions, and datasets
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
# load anova results
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/felzartamab_ARTanova_baseline.RData")


# CREATE FLEXTABLE OF ART MODELS ####
title_art <- paste("Table i. Non-parametric ANOVA (ART) of molecular scores in baseline biopsies from treated vs untreated patients")

res_art_flextable <- felzartamab_ARTanova_baseline %>%
    dplyr::select(annotation, n_cat, score, art_aov) %>%
    unnest(everything()) %>%
    dplyr::rename(p.value = `Pr(>F)`) %>%
    mutate(FDR = p.value %>% p.adjust(method = "fdr"), .by = annotation) %>%
    dplyr::select(annotation, score, Term, p.value, FDR) %>%
    flextable::flextable() %>%
    flextable::add_header_row(values = rep(title_art, ncol_keys(.))) %>%
    flextable::merge_h(part = "header") %>%
    flextable::merge_v(j = 1:2) %>%
    flextable::fontsize(size = 8, part = "all") %>%
    flextable::align(align = "center", part = "all") %>%
    flextable::bg(bg = "white", part = "all") %>%
    flextable::bg(i = ~ FDR < 0.05 & Term == "Felzartamab", j = 3:5, bg = "#fbff00") %>%
    flextable::colformat_double(j = 2:3, digits = 2) %>%
    flextable::colformat_double(j = 4:ncol_keys(.), digits = 3) %>%
    flextable::border_remove() %>%
    flextable::bold(part = "header") %>%
    flextable::padding(padding = 0, part = "all") %>%
    flextable::border(border = fp_border(), part = "all") %>%
    flextable::autofit()

# res_art_flextable %>% print(preview = "pptx")


# FORMAT TABLES OF CONTRASTS ####
# TODO: can this be simplified?
data_formatted_medians <- felzartamab_ARTanova_baseline %>%
    dplyr::select(annotation, score, art_aov_tidy, medians) %>%
    unnest(c(art_aov_tidy), names_repair = tidyr_legacy) %>%
    unnest(medians, names_repair = tidyr_legacy) %>%
    dplyr::rename(p = p.value)  %>% 
    mutate(FDR = p %>% p.adjust(method = "fdr"), .by = annotation) %>%
    mutate_at(
        vars("p", "FDR"),
        ~ ifelse(
            . < 0.01,
            formatC(., digits = 0, format = "e"),
            formatC(., digits = 3, format = "f")
        )
    ) %>%
    mutate(
        median = paste(
            format(round(median, 2), nsmall = 1),
            "\u00B1",
            round(IQR, 2),
            # Letter,
            sep = " "
        )
    ) %>%
    dplyr::select(
        annotation, score, Felzartamab, median, p, FDR
    ) %>%
    distinct(score, Felzartamab, .keep_all = TRUE) %>%
    nest(.by = c(annotation, score)) %>%
    mutate(
        data = map(
            data,
            function(data) {
                data %>%
                    pivot_wider(names_from = Felzartamab, values_from = median)  %>% 
                    relocate(c("p", "FDR"), .after = last_col())
            }
        )
    )


# UNIVERSAL VARIABLES FOR FLEXTABLE ####
title_art_baseline <- paste("Table i. Median \u00B1 IQR molecular scores in baseline biopsies from Felzartamab-treated vs placebo-treated patients")

footnoteText <- c(
    paste(
        "Grey shading denotes ANOVA effect FDR < 0.05\n",
        "FDR correction was carried out within each annotation grouping",
        sep = ""
    ) # ,
    # "Scores are the mean fold change in expression for all Probes within the data vs the mean expression for all Probes in the NoPGD control data"
)

header1 <- c(
    "Annotation", "Score",
    rep("Baseline", 4)
)

header2 <- c(
    "Annotation", "Score",
    rep(c("Placebo\n(N=10)", "Felzartamab\n(N=10)", "p", "FDR"))
)

cellWidths <- c(4, 11, rep(c(3.5), 4))

category_vec1 <- str_remove(data_formatted_medians$score, "Prob\\)")
category_vec2 <- str_extract(data_formatted_medians$score, "Prob")
category_vec3 <- str_replace(category_vec2, "Prob", ")")


# CREATE FELXTABLE OF BASELINE ART MODELS ####
flextable_baseline <- data_formatted_medians %>%
    unnest(everything(), names_repair = tidyr_legacy) %>%
    flextable::flextable() %>%
    flextable::compose(
        j = "score",
        value = as_paragraph(category_vec1, as_sub(category_vec2), category_vec3)
    ) %>%
    flextable::delete_part("header") %>%
    flextable::add_header_row(top = TRUE, values = header2) %>%
    # flextable::add_header_row(top = TRUE, values = header2) %>%
    flextable::add_header_row(top = TRUE, values = rep(title_art_baseline, ncol_keys(.))) %>%
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
    flextable::font(fontname = "Arial", part = "all") %>%
    flextable::fontsize(size = 8, part = "all") %>%
    flextable::fontsize(size = 8, part = "footer") %>%
    flextable::fontsize(i = 1, size = 12, part = "header") %>%
    flextable::bold(part = "header") %>%
    flextable::bold(j = 1, part = "body") %>%
    flextable::bg(bg = "white", part = "all") %>%
    flextable::bg(i = ~ as.numeric(FDR) < 0.05, j = 2:ncol_keys(.), bg = "grey90", part = "body") %>%
    # flextable::bold(i = ~ as.numeric(FDR_interaction) < 0.05, j = ncol_keys(.), part = "body") %>%
    flextable::padding(padding = 0, part = "all") %>%
    flextable::width(width = cellWidths, unit = "cm") %>%
    flextable::width(., width = dim(.)$widths * 30 / (flextable_dim(.)$widths), unit = "cm")


# PRINT THE FLEXTABLES ####
flextable_baseline %>% print(preview = "pptx")
