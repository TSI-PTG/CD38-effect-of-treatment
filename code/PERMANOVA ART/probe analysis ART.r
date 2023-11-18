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
library(genefilter)
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
# Suppress pesky dplyr reframe info
options(dplyr.reframe.inform = FALSE)
# laod reference set
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0Georg CD38 Vienna/G_Rstuff/data/Vienna44_18Oct23.RData")
# load affymap
load("Z:/DATA/Datalocks/Other data/affymap219_21Oct2019_1306_JR.RData")
# load reference expression data
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0Georg CD38 Vienna/G_Rstuff/data/mean_expression_K5086_MMDx.RData")


# IQR FILTER THE DATA ####
f1 <- function(x) (IQR(x) > 0.5)
ff <- filterfun(f1)
if (!exists("selected")) {
    selected <- genefilter(Vienna44, ff)
}
set00 <- Vienna44[selected, ]



# DEFINE SEED ####
seed <- 42


# DEFINE FEATURES ####
features <- set00 %>% featureNames()


# DEFINE THE SET ####
set <- set00 %>% tidy(addPheno = TRUE)


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
    dplyr::select(Patient, Group, Felz, Group_Felz, AffyID, value) %>%
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
art_aovF <- function(data, AffyID, p) {
    p <- progressor(along = AffyID)
    future_pmap(
        list(data, AffyID),
        function(data, AffyID) {
            p(sprintf("probe=%s", AffyID))
            set.seed(seed)
            art(value ~ Group * Felz + (1 | Patient), data = data) %>%
                anova(type = "II") %>%
                tibble() %>%
                mutate(part_eta_sq = F * Df / (F * Df + Df.res))
        },
        .options = furrr_options(seed = TRUE)
    )
}


# FUNCTION TO ITERATE ESTIMATED EFFECTS OF INTERACTION TERM ####
art_lmF <- function(data, AffyID, p) {
    p <- progressor(along = AffyID)
    future_pmap(
        list(data, AffyID),
        function(data, AffyID) {
            p(sprintf("probe=%s", AffyID))
            set.seed(seed)
            art(value ~ Group * Felz + (1 | Patient), data = data) %>%
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


# FUNCTION TO ITERATE GROUPWISE MEANS ####
meansF <- function(data, AffyID, p) {
    p <- progressor(along = AffyID)
    future_pmap(
        list(data, AffyID),
        function(data, AffyID) {
            p(sprintf("probe=%s", AffyID))
            set.seed(seed)
            groupwiseMean(
                value ~ Group + Felz,
                data = data
            ) %>%
                mutate(Group_Felz = paste(Group, Felz, sep = "_")) %>%
                dplyr::select(Group_Felz, Mean) %>%
                pivot_wider(names_from = Group_Felz, values_from = Mean) %>%
                mutate_all(~ 2^. %>% round(0))
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
res_art01 <- res_art00 %>% mutate(art_aov = art_aovF(data, AffyID))
res_art02 <- res_art01 %>% mutate(art_lm = art_lmF(data, AffyID))
res_art03 <- res_art02 %>% mutate(means = meansF(data, AffyID))

# FREE WORKERS FROM PARALLEL PROCESSING
plan(sequential)


# FORMAT THE UNIVARIATE RESULTS ####
res_art_table <- res_art03 %>%
    left_join(
        .,
        means_K5086 %>% dplyr::select(-Symb, -Gene, -PBT),
        by = "AffyID"
    ) %>%
    dplyr::select(AffyID, Symb, Gene, PBT, means, art_aov, art_lm) %>%
    unnest(everything()) %>%
    dplyr::filter(Term %>% str_detect(":")) %>%
    dplyr::rename(
        p.value = `Pr(>F)`
    ) %>%
    mutate(
        Term = "Interaction",
        FDR = p.value %>% p.adjust("fdr")
    ) %>%
    mutate_at(
        vars(p.value, FDR),
        ~ ifelse(
            . < 0.01,
            formatC(., digits = 0, format = "e"),
            formatC(., digits = 3, format = "f")
        )
    ) %>%
    mutate_if(
        is.numeric, ~ round(., 2)
    ) %>%
    dplyr::select(
        AffyID, Symb, Gene, PBT,
        Index_NoFelz, Index_Felz, FU1_NoFelz, FU1_Felz,
        Estimate, p.value, FDR, contains("MMDx")
    ) %>%
    arrange(p.value) %>%
    expand_grid(direction = c("all", "negative", "positive")) %>%
    group_by(direction) %>%
    nest() %>%
    tibble() %>%
    mutate(data = pmap(
        list(direction, data),
        function(direction, data) {
            if (direction == "all") {
                data
            } else if (direction == "negative") {
                data %>% dplyr::filter(Estimate < 0)
            } else if (direction == "positive") {
                data %>% dplyr::filter(Estimate > 0)
            }
        }
    ))


# GLOBALS VARIABLES FOR FLEXTABLES ####
title <- paste("Table i. Top 20 probesets affected by Felzartamab treatment (by non-parametric ANOVA p-value)")
cellWidths <- c(2, 2, 8, 4, rep(1.5, 9))
cellWidths %>% length()

header1 <- c(
    # "AffyID",
    "Symb", "Gene", "PBT", "Estimate", "P", "FDR",
    rep("Felzartamab study", 4),
    rep("K5086 reference set", 6)
)
header2 <- c(
    # "AffyID",
    "Symb", "Gene", "PBT", "Estimate", "P", "FDR",
    rep("Mean expression by group", 4),
    rep("Mean expression by MMDx", 6)
)

header3 <- c(
    # "AffyID",
    "Symb", "Gene", "PBT", 
    "Estimate", "P", "FDR",
    "Index\n(N=11)", "FU1\n(N=11)", "Index\n(N=11)", "FU1\n(N=11)",
    "Term", "Estimate", "F", "p.value", "FDR", "part_eta_sq"
)


# MAKE FLEXTABLE OF UNIVARIATE RESULTS ####
res_art_flextable <- res_art_table %>%
    mutate(flextable = pmap(
        list(direction, data),
        function(direction, data) {
            data %>%
                dplyr::slice(1:20) %>%
                flextable::flextable() %>%
                flextable::delete_part("header") %>%
                flextable::add_header_row(top = TRUE, values = header3) %>%
                flextable::add_header_row(top = TRUE, values = header2) %>%
                flextable::add_header_row(top = TRUE, values = header1) %>%
                flextable::add_header_row(top = TRUE, values = rep(title, ncol_keys(.))) %>%
                flextable::merge_h(part = "header") %>%
                flextable::merge_v(part = "header") %>%
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
                flextable::width(., width = dim(.)$widths * 33 / (flextable_dim(.)$widths), unit = "cm")
        }
    ))

res_art_flextable %>%
    dplyr::filter(direction == "negative") %>%
    pull(flextable) %>%
    print(preview = "pptx")


# EXPORT THE DATA AS AN EXCEL SHEET ####
savedir1 <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0Georg CD38 Vienna/G_Rstuff/output/"

names(res_art_table$data) <- res_art_table$direction
openxlsx::write.xlsx(res_art_table$data,
    asTable = TRUE,
    file = paste(savedir1, "all_probes_ANOVAs_Vienna44_18Oct23 ",
        # Sys.Date(),
        # format(Sys.time(), "_%I%M%p"),
        ".xlsx",
        sep = ""
    )
)
