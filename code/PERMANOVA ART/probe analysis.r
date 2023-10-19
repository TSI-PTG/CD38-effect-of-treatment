# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(ggbeeswarm) # install.packages("ggbeeswarm")
library(ggpubr) # install.packages("ggpubr")
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
library(furrr) # install.packages("furrr")
library(progressr) # install.packages("progressr")
library(parallel) # install.packages("parallel")
library(openxlsx) # install.packages("openxlsx")
# Bioconductor libraries
library(Biobase) # BiocManager::install("Biobase")
library(biobroom) # BiocManager::install("biobroom")
# Custom operators, functions, and datasets
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
se <- function(x) sd(x) / sqrt(length((x)))
log10zero <- scales::trans_new(
    name = "log10zero",
    transform = function(x) log10(x + 0.001),
    inverse = function(x) 10^x - 0.001
)
# SET UP PROGRESS BAR HANDLER
plan(sequential)
handlers(global = TRUE)
handlers(list(handler_progress(
    format = ":message [:current/:total] [:bar] [:percent in :elapsed] [:eta remaining]",
    clear = FALSE,
    intrusiveness = 1,
    complete = "+"
)))
# Suppress pesky dplyr summarise info
options(dplyr.summarise.inform = FALSE)
# laod reference set
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0Georg CD38 Vienna/G_Rstuff/data/Vienna44_18Oct23.RData")
# load affymap
load("Z:/DATA/Datalocks/Other data/affymap219_21Oct2019_1306_JR.RData")


# DEFINE SEED ####
seed <- 42


# DEFINE FEATURES ####
features <- Vienna44 %>% featureNames()


# DEFINE THE SET ####
set <- Vienna44 %>% tidy(addPheno = TRUE)


# WRANGLE THE PHENOTYPE DATA ####
df00 <- set %>%
    dplyr::rename(
        Patient = STUDY_EVALUATION_ID,
        Felz = Felzartamab_presumed
    ) %>%
    mutate(
        Patient = Patient %>% factor(),
        Group = Group %>% factor(levels = c("Index", "FU1")),
        Felz = Felz %>% factor(labels = c("NoFelz", "Felz")),
        TxBx = TxBx %>% as.numeric(),
        Group_Felz = paste(Group, Felz, sep = ":") %>%
            factor(levels = c("Index:NoFelz", "FU1:NoFelz", "Index:Felz", "FU1:Felz"))
    ) %>%
    dplyr::select(
        CEL, Patient, Center, Group, Felz, Group_Felz, TxBx, gene, value
    ) %>%
    tibble() %>%
    pivot_wider(names_from = gene, values_from = value)


# PERMANOVA ####
# permanova assumptions
# df_PBT_popmean %>%
#     dplyr::select(all_of(wanted_pbts)) %>%
#     boxplot()
df_features <- df00 %>% dplyr::select(all_of(features))
Group <- df00$Group
Felz <- df00$Felz
Group_Felz <- df00$Group_Felz

# res <- betadisper(dist(df_features), df00$Group_Felz)
# permutest(res, pairwise = T) # less bad
# pstat <- permutest(res, pairwise = TRUE) %>% permustats()
# densityplot(pstat, scales = list(x = list(relation = "free")))
# qqmath(pstat, scales = list(relation = "free"))
# looks like there are some outliers

# adonis
set.seed(seed)
res_adonis_interaction <- adonis2(
    df_features ~ Group * Felz,
    data = df00,
    method = "euclidean",
    # by = "margin", only specific margin if the sample sizes are unequal
    permutations = 100000
)

# set.seed(seed)
# res_adonis <- adonis2(
#     df_features ~ Group + Felz,
#     data = df00,
#     method = "euclidean",
#     by = "margin",
#     permutations = 10000
# )

# pairwise adonis
set.seed(seed)
res_pairwise_adonis <- pairwise.adonis(
    df_features,
    Group_Felz,
    reduce = "Group|Felz",
    sim.method = "euclidean",
    p.adjust.m = "fdr",
    perm = 10000
)


# PRODUCE TABLE OF ADONIS RESULTS ####
title_adonis <- paste("Table i. PERMANOVA of transcript expression in biopsies from treated vs untreated patients")

res_adonis_flextable <- res_adonis_interaction %>%
    tidy() %>%
    suppressWarnings() %>%
    dplyr::filter(term %nin% c("Residual", "Total")) %>%
    dplyr::select(-df, -SumOfSqs) %>%
    mutate(p.value = p.value %>% formatC(digits = 3, format = "g")) %>%
    flextable::flextable() %>%
    flextable::add_header_row(values = rep(title_adonis, ncol_keys(.))) %>%
    flextable::merge_h(part = "header") %>%
    flextable::fontsize(size = 8, part = "all") %>%
    flextable::align(align = "center", part = "all") %>%
    flextable::bg(bg = "white", part = "all") %>%
    flextable::colformat_double(j = 2:3, digits = 2) %>%
    flextable::bg(i = ~ p.value %>% as.numeric() < 0.05, bg = "#fbff00", part = "body") %>%
    flextable::border_remove() %>%
    flextable::bold(part = "header") %>%
    flextable::padding(padding = 0, part = "all") %>%
    flextable::border(border = fp_border(), part = "all") %>%
    flextable::autofit()

res_adonis_flextable %>% print(preview = "docx")


# PRODUCE TABLE OF PAIRWISE ADONIS RESULTS ####
title_adonis_pairwise <- paste("Table i. Pairwise PERMANOVA of transcript expression in biopsies from treated vs untreated patients")

res_pairwise_adonis_flextable <- res_pairwise_adonis %>%
    dplyr::select(-Df, -sig, -SumsOfSqs) %>%
    dplyr::rename("F-value" = F.Model, FDR = p.adjusted) %>%
    mutate(FDR = FDR %>% formatC(digits = 3, format = "g")) %>%
    arrange(p.value) %>%
    flextable::flextable() %>%
    flextable::add_header_row(values = rep(title_adonis_pairwise, ncol_keys(.))) %>%
    flextable::merge_h(part = "header") %>%
    flextable::fontsize(size = 8, part = "all") %>%
    flextable::align(align = "center", part = "all") %>%
    flextable::bg(bg = "white", part = "all") %>%
    flextable::bg(i = ~ FDR %>% as.numeric() < 0.05, bg = "#fbff00") %>%
    flextable::colformat_double(j = 2:3, digits = 2) %>%
    flextable::border_remove() %>%
    flextable::bold(part = "header") %>%
    flextable::padding(padding = 0, part = "all") %>%
    flextable::border(border = fp_border(), part = "all") %>%
    flextable::autofit()

res_pairwise_adonis_flextable %>% print(preview = "docx")


# WRANGLE THE DATA FOR UNIVARIATE TESTS ####
df01 <- df00 %>%
    pivot_longer(
        cols = c(
            all_of(features)
        ),
        names_to = "AffyID",
        values_to = "value"
    ) %>%
    group_by(AffyID) %>%
    nest() %>%
    tibble() %>%
    right_join(
        affymap219 %>%
            dplyr::select(AffyID, Symb, Gene, PBT),
        .,
        by = "AffyID"
    ) %>%
    tibble() %>%
    mutate(id = row_number())


# FUNCTION TO ITERATE ALIGNED RANK TRANSFORM ANOVA ####
artF <- function(data, AffyID, p) {
    p <- progressor(along = AffyID)
    future_pmap(
        list(data, AffyID),
        function(data, AffyID) {
            p(sprintf("probe=%s", AffyID))
            set.seed(seed)
            art(value ~ Group * Felz + (1 | Patient), data = data)
        },
        .options = furrr_options(seed = TRUE)
    )
}


# FUNCTION TO ITERATE ANOVA TABLE ON ART MODEL ####
art_aovF <- function(art, AffyID, p) {
    p <- progressor(along = AffyID)
    future_pmap(
        list(art, AffyID),
        function(art, AffyID) {
            p(sprintf("probe=%s", AffyID))
            set.seed(seed)
            art_aov <- art %>% anova(type = "II")
            art_aov$part_eta_sq <- with(art_aov, `F` * `Df` / (`F` * `Df` + `Df.res`))
            art_aov %>% tibble
        },
        .options = furrr_options(seed = TRUE)
    )
}


# FUNCTION TO ITERATE ESTIMATED EFFECTS OF INTERACTION TERM ####
art_lmF <- function(art, AffyID, p) {
    p <- progressor(along = AffyID)
    future_pmap(
        list(art, AffyID),
        function(art, AffyID) {
            p(sprintf("probe=%s", AffyID))
            set.seed(seed)
            art %>%
                artlm(term = "Group:Felz", response = "aligned") %>%
                summary() %>%
                pluck("coefficients") %>%
                as.data.frame() %>%
                rownames_to_column("Term_lm") %>%
                dplyr::filter(Term_lm != "(Intercept)")
        },
        .options = furrr_options(seed = TRUE)
    )
}



# RUN ART #####
# SET UP PARALLEL PROCESSING
# n_workers <- detectCores()
plan(multisession, workers = 10) # select the number of workers/cores

message("can take some time to initialize\n")
res_art00 <- df01 %>% mutate(art = artF(data, AffyID))
res_art01 <- res_art00 %>% mutate(art_aov_tidy = art_aovF(art, AffyID))
res_art02 <- res_art01 %>% mutate(art_lm = art_lmF(art, AffyID))

# FREE WORKERS FROM PARALLEL PROCESSING
plan(sequential)




# FORMAT THE UNIVARIATE RESULTS ####
res_art_table <- res_art02 %>%
    dplyr::select(AffyID, Symb, Gene, PBT, art_aov_tidy, art_lm) %>%
    unnest(everything()) %>%
    dplyr::filter(Term %>% str_detect(":")) %>%
    dplyr::rename(p.value = `Pr(>F)`) %>%
    dplyr::select(AffyID, Symb, Gene, PBT, Term, Estimate, F, p.value, part_eta_sq) %>%
    arrange(p.value) %>%
    mutate(
        Term = "Interaction",
        p.value = ifelse(
        p.value < 0.001,
        formatC(p.value, digits = 1, format = "e"),
        round(p.value, 4)
    )) %>%
    mutate_if(
        is.numeric, ~ round(., 2)
    )



# MAKE FLEXTABLE OF UNIVARIATE RESULTS ####
title <- paste("Table i. Top 20 probesets affected by Felzartamab treatment (by non-parametric ANOVA p-value)")
cellWidths <- c(2, 2, 8, 4, rep(1.5, 4))


res_art_flextable <- res_art_table %>%
    dplyr::slice(1:20) %>%
    flextable::flextable() %>%
    flextable::add_header_row(top = TRUE, values = rep(title, ncol_keys(.))) %>%
    flextable::merge_h(part = "header") %>%
    flextable::border_remove() %>%
    flextable::border(part = "header", border = fp_border()) %>%
    flextable::border(part = "body", border = fp_border()) %>%
    flextable::align(align = "center", part = "all") %>%
    flextable::font(fontname = "Arial", part = "all") %>%
    flextable::fontsize(size = 8, part = "all") %>%
    flextable::fontsize(i = 1, size = 12, part = "header") %>%
    flextable::bold(part = "header") %>%
    flextable::bg(bg = "white", part = "all") %>%
    flextable::padding(padding = 0, part = "all") %>%
    flextable::width(width = cellWidths, unit = "cm") %>%
    flextable::width(., width = dim(.)$widths * 26.25 / (flextable_dim(.)$widths), unit = "cm")


res_art_flextable %>% print(preview = "pptx")


# SAVE THE EXCEL FILE ####
write.xlsx(res_art_table, file = "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0Georg CD38 Vienna/G_Rstuff/output/all_probes_ANOVAs_Vienna44_18Oct23.xlsx")
