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
# Bioconductor libraries
library(Biobase) # BiocManager::install("Biobase")
# Custom operators, functions, and datasets
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
se <- function(x) sd(x) / sqrt(length((x)))
log10zero <- scales::trans_new(
    name = "log10zero",
    transform = function(x) log10(x + 0.001),
    inverse = function(x) 10^x - 0.001
)
# Suppress pesky dplyr summarise info
options(dplyr.summarise.inform = FALSE)
# laod reference set
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0Georg CD38 Vienna/G_Rstuff/data/Vienna44_18Oct23.RData")



# DEFINE CATEGORIES FOR FEATURES ####
Rejectionrelated <- c("GRIT3", "Rej-RAT")
ABMRrelated <- c("DSAST", "NKB", "ABMRpm", "ggt0", "cggt0", "ptcgt0")
TCMRrelated <- c("TCMR-RAT", "QCAT", "TCB", "TCMRt", "tgt1", "igt1")
Endothelium <- c("ENDAT")
Parenchyma <- c("KT1", "KT2")
Macrophage <- c("AMAT1", "QCMAT")
Injuryrecent <- c("FICOL", "IRRAT30", "IRITD3", "IRITD5")
Injurylate <- c("IGT", "MCAT", "BAT", "cigt1", "ctgt1")


# DEFINE FEATURES ####
features <- c(Rejectionrelated, ABMRrelated)


# DEFINE THE SET ####
set <- Vienna44


# WRANGLE THE PHENOTYPE DATA ####
df00 <- set %>%
    pData() %>%
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
        CEL, Patient, Center, Group, Felz, Group_Felz, TxBx,
        all_of(features)
    ) %>%
    tibble()


# COMPARE THE MEAN variable SCORES ####
df00 %>%
    group_by(Group_Felz) %>%
    summarise_at(vars(all_of(features)), mean) %>%
    column_to_rownames("Group_Felz") %>%
    t() %>%
    round(3)


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
set.seed(42)
res_adonis_interaction <- adonis2(
    df_features ~ Group * Felz,
    data = df00,
    method = "euclidean",
    # by = "margin", only specific margin if the sample sizes are unequal
    permutations = 100000
)

# set.seed(42)
# res_adonis <- adonis2(
#     df_features ~ Group + Felz,
#     data = df00,
#     method = "euclidean",
#     # by = "margin",
#     permutations = 10000
# )

# pairwise adonis
set.seed(42)
res_pairwise_adonis <- pairwise.adonis(
    df_features,
    Group_Felz,
    reduce = "Group|Felz",
    sim.method = "euclidean",
    p.adjust.m = "fdr",
    perm = 10000
)


# PRODUCE TABLE OF ADONIS RESULTS ####
title_adonis <- paste("Table i. PERMANOVA of molecular scores in biopsies from treated vs untreated patients")

# res_adonis_flextable <- bind_rows(
#     res_adonis %>%
#         tidy() %>%
#         suppressWarnings(),
#     res_adonis_interaction %>%
#         tidy() %>%
#         suppressWarnings()
# ) %>%
res_adonis_flextable <- res_adonis_interaction %>%
        tidy() %>%
        suppressWarnings() %>%
    dplyr::filter(term %nin% c("Residual", "Total")) %>%
    dplyr::select(-df, -SumOfSqs) %>%
    mutate(p.value = p.value %>% formatC(digits = 3, format = "fg")) %>%
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

res_adonis_flextable  %>% print(preview = "docx")


# PRODUCE TABLE OF PAIRWISE ADONIS RESULTS ####
title_adonis_pairwise <- paste("Table i. Pairwise PERMANOVA of molecular scores in biopsies from treated vs untreated patients")

res_pairwise_adonis_flextable <- res_pairwise_adonis %>%
    dplyr::select(-Df, -sig, -SumsOfSqs) %>%
    dplyr::rename("F-value" = F.Model, FDR = p.adjusted) %>%
    mutate(
        p.value = p.value %>% formatC(digits = 0, format = "e"),
        FDR = FDR %>% formatC(digits = 0, format = "e"),
    ) %>%
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
        names_to = "variable",
        values_to = "value"
    ) %>%
    group_by(variable) %>%
    nest() %>%
    tibble() %>%
    mutate(
        variable = variable %>%
            factor(
                levels = c(
                    "Rej-RAT", "GRIT3", "ptcgt0", "BAT",
                    "ABMRpm", "ggt0", "cggt0", "NKB", "DSAST",
                    "TCMRt", "tgt1", "igt1", "TCB", "TCMR-RAT", "QCAT",
                    "AMAT1", "QCMAT",
                    "FICOL", "IRRAT30", "IRITD3", "IRITD5",
                    "cigt1", "ctgt1", "IGT", "MCAT",
                    "KT1", "KT2"
                )
            ),
        annotation = case_when(
            variable %in% Rejectionrelated ~ "Rejection-related",
            variable %in% TCMRrelated ~ "TCMR-related",
            variable %in% ABMRrelated ~ "ABMR-related",
            # variable %in% Endothelium ~ "Endothelium-related",
            variable %in% Parenchyma ~ "Parenchyma-related",
            variable %in% Macrophage ~ "Macrophage-related",
            variable %in% Injurylate ~ "Atrophy-fibrosis-related",
            variable %in% Injuryrecent ~ "Recent injury-related",
            TRUE ~ " "
        ) %>%
            factor(
                levels = c(
                    "Rejection-related", "ABMR-related", "TCMR-related",
                    "Macrophage-related", "Recent injury-related", "Atrophy-fibrosis-related",
                    "Parenchyma-related"
                )
            ),
        score = case_when(
            variable == "TCMRt" ~ "TCMR classifier (TCMRProb)",
            variable == "TCB" ~ "T-cell burden (TCB)",
            variable == "TCMR-RAT" ~ "TCMR-associated RATs (TCMR-RAT)",
            variable == "tgt1" ~ "Tubulitis classifier (t>1Prob)",
            variable == "igt1" ~ "Interstitial infiltrate classifier (i>1Prob)",
            variable == "QCAT" ~ "Cytotoxic T cell-associated transcripts (QCAT)",
            variable == "Rej-RAT" ~ "Rejection-associated transcripts (Rej-RAT)",
            variable == "GRIT3" ~ "Interferon gamma-inducible transcripts (GRIT3)",
            variable == "ABMR-RAT" ~ "ABMR-associated RATs (ABMR-RAT)",
            variable == "ABMRpm" ~ "ABMR classifier (ABMRProb)",
            variable == "DSAST" ~ "DSA-selective transcripts (DSAST)",
            variable == "NKB" ~ "NK cell burden (NKB)",
            variable == "ggt0" ~ "Glomerulitis classifier (g>0Prob)",
            variable == "cggt0" ~ "Glomerular double contours classifier (cg>0Prob) ",
            variable == "ptcgt0" ~ "Peritubular capillaritis classifier (ptc>0Prob)",
            variable == "AMAT1" ~ "Alternatively activated macrophage (AMAT1)",
            variable == "QCMAT" ~ "Constitutive macrophage (QCMAT)",
            variable == "FICOL" ~ "Fibrillar collagen (FICOL)",
            variable == "DAMP" ~ "Damage-associated molecular pattern transcripts (DAMP)",
            variable == "cIRIT" ~ "Cardiac injury and repair–induced transcripts (cIRIT)",
            variable == "IRITD3" ~ "Injury-repair induced, day 3 (IRITD3)",
            variable == "IRITD5" ~ "Injury-repair induced, day 5 (IRITD5)",
            variable == "IRRAT30" ~ "Injury-repair associated (IRRAT30)",
            variable == "cigt1" ~ "Fibrosis classifier (ci>1Prob)",
            variable == "ctgt1" ~ "Atrophy classifier (ct>1Prob)",
            variable == "txbxCorrCLAGdn_tbb" ~ "Transcripts downregulated in CLAD",
            variable == "txbxCorrCLAGup_tbb" ~ "Transcripts upregulated in CLAD",
            variable == "SFT" ~ "Surfactant-associated transcripts (SFT)",
            variable == "ENDAT" ~ "Endothelial cell-associated transcripts (ENDAT)",
            variable == "eDSAST" ~ "Endothelium-expressed DSA-selective transcripts (eDSAST)",
            variable == "IGT" ~ "Immunoglobulin transcripts (IGT)",
            variable == "BAT" ~ "B cell–associated transcripts (BAT)",
            variable == "MCAT" ~ "Mast cell-associated transcripts (MCAT)",
            variable == "KT1" ~ "Kidney parenchymal transcripts (KT1)",
            variable == "KT2" ~ "Kindey parenchymal transcripts - no solute carriers (KT2)"
        ), .before = 1
    ) %>%
    arrange(annotation, variable)
df01 %>% print(n = "all")


# UNIVARIATE NONPARAMETRIC TESTS ####
# univariate tests
df02 <- df01 %>%
    mutate(
        means = map(
            data,
            function(data) {
                data %>%
                    group_by(Group_Felz) %>%
                    summarise(mean = mean(value), se = se(value)) %>%
                    arrange(
                        Group_Felz %>%
                            factor(levels = c(
                                "Index:NoFelz", "FU1:NoFelz",
                                "Index:Felz", "FU1:Felz"
                            ))
                    )
            }
        ),
        medians = map(
            data,
            function(x) {
                x %>%
                    group_by(Group_Felz) %>%
                    summarise(median = median(value), IQR = IQR(value)) %>%
                    arrange(
                        Group_Felz %>%
                            factor(levels = c(
                                "Index:NoFelz", "FU1:NoFelz",
                                "Index:Felz", "FU1:Felz"
                            ))
                    )
            }
        ),
        art = map(data, ~ art(value ~ Group*Felz + (1 | Patient), data = .)),
        art_aov = map(art, ~ anova(., type = "II")),
        art_aov_tidy = map(art_aov, . %>% tidy()) %>% suppressWarnings(),
        art_con = map(art, ~ art.con(.x, "Group:Felz", adjust = "none")),
        art_con_tidy = map(art_con, tidy),
        art_con_cld = map(art_con, . %>% as.data.frame() %>%
            cldList(p.value ~ contrast, data = .) %>%
            arrange(Group %>%
                factor(
                    levels = c(
                        "Index,NoFelz", "FU1,NoFelz",
                        "Index,Felz", "FU1,Felz"
                    )
                )))
    )
names(df02$art_aov_tidy) <- df02$variable %>% as.character()
names(df02$art_con) <- df02$variable %>% as.character()
names(df02$art_con_tidy) <- df02$variable %>% as.character()
names(df02$art_con_cld) <- df02$variable %>% as.character()

df02$art_aov
df02$art_aov_tidy
df02$art_con_tidy


# CREATE FLEXTABLE OF ART MODELS ####
title_art <- paste("Table i. Non-parametric ANOVA (ART) of molecular scores in biopsies from treated vs untreated patients")
title_art_pairwise <- paste("Table i. Pairwise Non-parametric ANOVA (ART) of molecular scores in biopsies from treated vs untreated patients")

res_art_flextable <- df02 %>%
    dplyr::select(annotation, score, art_aov) %>%
    unnest(everything()) %>%
    dplyr::select(annotation, score, Term, `F`, `Pr(>F)`) %>%
    dplyr::rename(p.value = `Pr(>F)`) %>%
    mutate(FDR = p.value %>% p.adjust(method = "fdr")) %>%
    flextable::flextable() %>%
    flextable::add_header_row(values = rep(title_art, ncol_keys(.))) %>%
    flextable::merge_h(part = "header") %>%
    flextable::merge_v(j = 1:2) %>%
    flextable::fontsize(size = 8, part = "all") %>%
    flextable::align(align = "center", part = "all") %>%
        flextable::bg(bg = "white", part = "all") %>%
    flextable::bg(i = ~ FDR < 0.05 & Term == "Group:Felz", j = 3:6, bg = "#fbff00") %>%
    flextable::colformat_double(j = 2:3, digits = 2) %>%
    flextable::colformat_double(j = 4:ncol_keys(.), digits = 3) %>%
    flextable::border_remove() %>%
    flextable::bold(part = "header") %>%
    flextable::padding(padding = 0, part = "all") %>%
    flextable::border(border = fp_border(), part = "all") %>%
    flextable::autofit()

res_art_flextable %>% print(preview = "pptx")


# FORMAT TABLES OF PAIRWISE ART MODELS ####
res_art_pairwise_formatted <- df02 %>%
    dplyr::select(annotation, score, art_con_cld, medians) %>%
    unnest(c(art_con_cld, medians), names_repair = tidyr_legacy) %>%
    mutate(
        medians =
            paste(
                format(round(median, 2), nsmall = 1),
                "\u00B1",
                round(IQR, 2),
                Letter,
                sep = " "
            )
    ) %>%
    dplyr::select(annotation, score, Group, medians) %>%
    pivot_wider(names_from = c(Group), values_from = c(medians)) %>%
    as.data.frame()

res_art_pairwise_formatted


# FORMAT variable TABLE ####
# define sample sizes
df_n <- df00 %>%
    group_by(Group_Felz) %>%
    tally()

N <- df_n %>% pull(n)
cellWidths <- c(5, 9.25, rep(3, 4))
title_art_pairwise <- paste("Table i. Pairwise Non-parametric ANOVA (ART) of molecular scores in biopsies from treated vs untreated patients")

footnoteText <- c(
    paste(
        "Yellow shading indicates significant interactive effect of treatment in follow-up biopsy.\n",
        "NOTE: Lettering represents contrasts within each row (i.e., multigene score).\n",
        "Groups sharing the same letters are not different from one another for that molecular score",
        sep = ""
    ) # ,
    # "Scores are the mean fold change in expression for all Probes within the set vs the mean expression for all Probes in the NoPGD control set"
)

header1 <- c("Annotation", "Score", "No Felz", "No Felz", "Felz", "Felz")
header2 <- c("Annotation", "Score", "Index", "FU1", "Index", "FU1")
header3 <- c("Annotation", "Score", paste("N = ", N, sep = ""))



# CREATE FELXTABLE OF PAIRWISE ART MODELS ####
res_art_pairwise_flextable <- res_art_pairwise_formatted %>%
    flextable::flextable() %>%
    flextable::delete_part("header") %>%
    flextable::add_header_row(top = TRUE, values = header3) %>%
    flextable::add_header_row(top = TRUE, values = header2) %>%
    flextable::add_header_row(top = TRUE, values = header1) %>%
    flextable::add_header_row(top = TRUE, values = rep(title_art_pairwise, ncol_keys(.))) %>%
    flextable::add_footer_row(values = footnoteText[[1]], colwidths = ncol_keys(.)) %>%
    flextable::merge_v(j = 1:2) %>%
    flextable::merge_v(part = "header") %>%
    flextable::merge_h(i = 1:3, part = "header") %>%
    flextable::border_remove() %>%
    flextable::border(part = "header", border = fp_border()) %>%
    flextable::border(part = "body", border = fp_border()) %>%
    flextable::border(part = "footer", border.left = fp_border(), border.right = fp_border()) %>%
    flextable::border(i = 1, part = "footer", border.bottom = fp_border()) %>%
    flextable::align(align = "center") %>%
    flextable::align(align = "center", part = "header") %>%
    flextable::font(fontname = "Arial", part = "all") %>%
    flextable::fontsize(size = 10, part = "all") %>%
    flextable::fontsize(size = 10, part = "footer") %>%
    flextable::fontsize(i = 1, size = 12, part = "header") %>%
    flextable::bold(part = "header") %>%
    flextable::bold(j = 1, part = "body") %>%
    flextable::bg(bg = "white", part = "all") %>%
    flextable::bg(i = ~ score %>% str_detect("NKB|DSAST|Glomerulitis classifier"), j = 2:6, bg = "yellow", part = "body") %>%
    # flextable::bg(i = ~ score %>% str_detect(""), j = 2:6, bg = "yellow", part = "body") %>%
    flextable::padding(padding = 0, part = "all") %>%
    flextable::width(width = cellWidths, unit = "cm") %>%
    flextable::width(., width = dim(.)$widths * 26.25 / (flextable_dim(.)$widths), unit = "cm")

res_art_pairwise_flextable


# PRINT THE FLEXTABLES ####
res_art_pairwise_flextable %>% print(preview = "pptx")



