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
# LOAD REFERENCE SET(S)
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0Georg CD38 Vienna/G_Rstuff/data/Vienna44_18Oct23.RData")
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0Georg CD38 Vienna/G_Rstuff/data/Vienna.RData") # N=52


# DEFINE FEATURES ####
wanted_pbts <- c(
    "TCB", "DSAST", "NKB", "Rej-RAT", "GRIT3",
    "QCAT", "IGT", "AMAT1", "BAT", "IRRAT30",
    "IRITD3", "IRITD5", "MCAT", "QCMAT", "KT1", "KT2"
)
wanted_classifiers <- c(
    "TCMRt", "igt1", "tgt1", "cggt0", "ABMRpm", "ggt0", "ptcgt0"
)
features <- c(wanted_pbts, wanted_classifiers)


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


# COMPARE THE MEAN PBT SCORES ####
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

# res <- betadisper(dist(PBTs), df_PBT_popmean$BLADgrade_1yr)
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
    by = "margin",
    permutations = 10000
)
set.seed(42)
res_adonis <- adonis2(
    df_features ~ Group + Felz,
    data = df00,
    method = "euclidean",
    by = "margin",
    permutations = 10000
)

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

res_adonis_flextable <- bind_rows(
    res_adonis %>%
        tidy() %>%
        suppressWarnings(),
    res_adonis_interaction %>%
        tidy() %>%
        suppressWarnings()
) %>%
    dplyr::filter(term %nin% c("Residual", "Total")) %>%
    dplyr::select(-df, -SumOfSqs) %>%
    mutate(p.value = p.value %>% formatC(digits = 0, format = "e")) %>%
    flextable::flextable() %>%
    flextable::add_header_row(values = rep(title_adonis, ncol_keys(.))) %>%
    flextable::merge_h(part = "header") %>%
    flextable::fontsize(size = 8, part = "all") %>%
    flextable::align(align = "center", part = "all") %>%
    flextable::bg(i = 1, bg = "grey80", part = "header") %>%
    flextable::colformat_double(j = 2:3, digits = 2) %>%
    flextable::bg(i = ~ p.value %>% as.numeric() < 0.05, bg = "#fbff00", part = "body") %>%
    flextable::border_remove() %>%
    flextable::bold(part = "header") %>%
    flextable::padding(padding = 0, part = "all") %>%
    flextable::border(border = fp_border(), part = "all") %>%
    flextable::autofit()

res_adonis_flextable


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
    flextable::bg(i = ~ FDR %>% as.numeric() < 0.05, bg = "#fbff00") %>%
    flextable::bg(i = 1, bg = "grey80", part = "header") %>%
    flextable::colformat_double(j = 2:3, digits = 2) %>%
    flextable::border_remove() %>%
    flextable::bold(part = "header") %>%
    flextable::padding(padding = 0, part = "all") %>%
    flextable::border(border = fp_border(), part = "all") %>%
    flextable::autofit()

res_pairwise_adonis_flextable


# UNIVARIATE NONPARAMETRIC TESTS ####
# univariate tests
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
        art = map(data, ~ art(value ~ Group * Felz + (1 | Patient), data = .)),
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
names(df01$art_aov_tidy) <- df01$variable
names(df01$art_con) <- df01$variable
names(df01$art_con_tidy) <- df01$variable

names(df01$art_con_cld) <- df01$variable

df01$means[[1]]
df01$medians[[1]]
df01$art_aov
df01$art_aov_tidy
df01$art_con_tidy
df01$art_con_cld[[1]]


# CREATE FLEXTABLE OF ART MODELS ####
title_art <- paste("Table i. Non-parametric ANOVA (ART) of molecular scores in biopsies from treated vs untreated patients")
title_art_pairwise <- paste("Table i. Pairwise Non-parametric ANOVA (ART) of of molecular scores in biopsies from treated vs untreated patients")

res_art_flextable <- df01 %>%
    dplyr::select(variable, art_aov) %>%
    unnest(everything()) %>%
    dplyr::select(variable, Term, `F`, `Pr(>F)`) %>%
    dplyr::rename(p.value = `Pr(>F)`) %>%
    mutate(FDR = p.value %>% p.adjust(method = "fdr")) %>%
    flextable::flextable() %>%
    flextable::add_header_row(values = rep(title_art, ncol_keys(.))) %>%
    flextable::merge_h(part = "header") %>%
    flextable::fontsize(size = 8, part = "all") %>%
    flextable::align(align = "center", part = "all") %>%
    flextable::bg(i = ~ FDR < 0.05, bg = "#fbff00") %>%
    flextable::bg(i = 1, bg = "grey80", part = "header") %>%
    flextable::colformat_double(j = 2:3, digits = 2) %>%
    flextable::colformat_double(j = 4:ncol_keys(.), digits = 3) %>%
    flextable::border_remove() %>%
    flextable::bold(part = "header") %>%
    flextable::padding(padding = 0, part = "all") %>%
    flextable::border(border = fp_border(), part = "all") %>%
    flextable::autofit()

res_art_flextable %>% print(preview = "docx")


# CREATE TABLES OF PAIRWISE ART MODELS ####

res_art_pairwise_flextable <- df01 %>%
                    dplyr::select(variable, art_con_cld, medians) %>%
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
                    dplyr::select(variable, Group, medians) %>%
                    pivot_wider(names_from = c(Group), values_from = c(medians)) %>%
                    as.data.frame() %>%
                    # dplyr::rename("NoPGD" = "0", "PGD"  = "1")  %>%
                    relocate(c("p.value", "FDR"), .after = PGD)

df_res_art$art_pairwise_table_medianIQR[[1]]


# FORMAT PBT TABLE ####
# define sample sizes
df_n <- PBTnest %>%
    group_by(PGDstatus, txBxBin) %>%
    tally() %>%
    arrange(txBxBin)

wanted_pbts
# define categories for the pbts
Rejectionrelated <- c("GRIT1", "Rej-RAT")
TCMRrelated <- c("TCMR-RAT", "QCAT", "TCB", "TCMRtPBT", "tgt1PBT", "igt1PBT")
ABMRrelated <- c("DSAST", "NKB", "ABMR-RAT", "ABMRpmPBT", "ggt0PBT", "ptcgt0PBT")
SFT <- c("SFT")
Endothelium <- c("ENDAT")
Parenchyma <- c("LT1", "LT2")
Macrophage <- c("AMAT1", "QCMAT")
Injuryrecent <- c("FICOL", "IRRAT30", "IRITD3", "IRITD5")
Injurylate <- c("IGT", "MCAT", "cigt1PBT", "ctgt1PBT")
CLAD <- c("txbxCorrCLAGdn_tbb", "txbxCorrCLAGup_tbb")

cellWidths <- c(5, 9.25, rep(3, 4))
title_art_pairwise <- paste("Table x. Median \u00B1 IQR PBT scores in PGD vs No PGD ")

footnoteText <- c(
    paste("NOTE: Lettering represents contrasts within each row (i.e., transcript set). ",
        "BLAD grades sharing the same letters are not different from one another for that molecular score",
        sep = ""
    ),
    "Scores are the mean fold change in expression for all Probes within the set vs the mean expression for all Probes in the NoPGD control set"
)

# demon table
df_res_art <- df_res_art %>%
    mutate(
        art_pairwise_table_medianIQR_formatted = map(
            art_pairwise_table_medianIQR,
            function(x) {
                x %>%
                    mutate(
                        category = case_when(
                            PBT %in% Rejectionrelated ~ "Rejection-related",
                            PBT %in% TCMRrelated ~ "TCMR-related transcripts",
                            PBT %in% ABMRrelated ~ "ABMR-related transcripts",
                            PBT %in% SFT ~ "Surfactant-related",
                            PBT %in% Endothelium ~ "Endothelium-related transcripts",
                            PBT %in% Parenchyma ~ "Parenchyma-related transcripts",
                            PBT %in% Macrophage ~ "Macrophage-related transcripts",
                            PBT %in% Injurylate ~ "Late injury-related transcripts (atrophy-fibrosis)",
                            PBT %in% Injuryrecent ~ "Recent injury-related transcripts",
                            TRUE ~ " "
                        ),
                        order = case_when(
                            PBT %in% Rejectionrelated ~ 7,
                            PBT %in% TCMRrelated ~ 8,
                            PBT %in% ABMRrelated ~ 9,
                            PBT %in% SFT ~ 5,
                            PBT %in% Endothelium ~ 4,
                            PBT %in% Parenchyma ~ 3,
                            PBT %in% Macrophage ~ 6,
                            PBT %in% Injurylate ~ 2,
                            PBT %in% Injuryrecent ~ 1,
                            TRUE ~ 0
                        ),
                        order2 = case_when(
                            PBT == "FICOL" ~ 1,
                            PBT == "DAMP" ~ 2,
                            PBT == "IRRAT30" ~ 3,
                            PBT == "IRITD3" ~ 4,
                            PBT == "IRITD5" ~ 5,
                            PBT == "cIRIT" ~ 6,
                            PBT == "cigt1PBT" ~ 7,
                            PBT == "ctgt1PBT" ~ 8,
                            PBT == "IGT" ~ 9,
                            PBT == "BAT" ~ 10,
                            PBT == "MCAT" ~ 11,
                            PBT == "LT1" ~ 12,
                            PBT == "LT2" ~ 13,
                            # PBT == 'LT3' ~ 14,
                            # PBT == 'LT4' ~ 15,
                            PBT == "eDSAST" ~ 14,
                            PBT == "ENDAT" ~ 15,
                            PBT == "SFT" ~ 16,
                            PBT == "AMAT1" ~ 17,
                            PBT == "QCMAT" ~ 18,
                            PBT == "GRIT1" ~ 19,
                            PBT == "Rej-RAT" ~ 20,
                            PBT == "TCMR-RAT" ~ 21,
                            PBT == "QCAT" ~ 22,
                            PBT == "TCB" ~ 23,
                            PBT == "TCMRtPBT" ~ 24,
                            PBT == "tgt1PBT" ~ 25,
                            PBT == "igt1PBT" ~ 26,
                            PBT == "ABMR-RAT" ~ 27,
                            PBT == "NKB" ~ 28,
                            PBT == "DSAST" ~ 29,
                            PBT == "ABMRpmPBT" ~ 30,
                            PBT == "ggt0PBT" ~ 31,
                            PBT == "ptcgt0PBT" ~ 32
                        ),
                        score = case_when(
                            PBT == "TCMRtPBT" ~ "TCMR classifier transcripts (TCMRProb)",
                            PBT == "TCB" ~ "T-cell burden (TCB)",
                            PBT == "TCMR-RAT" ~ "TCMR-associated RATs (TCMR-RAT)",
                            PBT == "tgt1PBT" ~ "Tubulitis classifier transcripts (t>1Prob)",
                            PBT == "igt1PBT" ~ "Interstitial infiltrate classifier transcripts (i>1Prob)",
                            PBT == "QCAT" ~ "Cytotoxic T cell-associated transcripts (QCAT)",
                            PBT == "Rej-RAT" ~ "Rejection-associated transcripts (Rej-RAT)",
                            PBT == "GRIT1" ~ "Interferon gamma-inducible transcripts (GRIT1)",
                            PBT == "ABMR-RAT" ~ "ABMR-associated RATs (ABMR-RAT)",
                            PBT == "ABMRpmPBT" ~ "ABMR classifier transcripts (ABMRProb)",
                            PBT == "DSAST" ~ "DSA-selective transcripts (DSAST)",
                            PBT == "NKB" ~ "NK cell burden (NKB)",
                            PBT == "ggt0PBT" ~ "Glomerulitis classifier transcripts (g>0Prob)",
                            PBT == "ptcgt0PBT" ~ "Capillaritis classifier transcripts (ptc>0Prob)",
                            PBT == "AMAT1" ~ "Alternatively activated macrophage (AMAT1)",
                            PBT == "QCMAT" ~ "Constitutive macrophage (QCMAT)",
                            PBT == "FICOL" ~ "Fibrillar collagen (FICOL)",
                            PBT == "DAMP" ~ "Damage-associated molecular pattern transcripts (DAMP)",
                            PBT == "cIRIT" ~ "Cardiac injury and repair–induced transcripts (cIRIT)",
                            PBT == "IRITD3" ~ "Injury-repair induced, day 3 (IRITD3)",
                            PBT == "IRITD5" ~ "Injury-repair induced, day 5 (IRITD5)",
                            PBT == "IRRAT30" ~ "Injury-repair associated (IRRAT30)",
                            PBT == "cigt1PBT" ~ "Fibrosis classifier transcripts (ci>1Prob)",
                            PBT == "ctgt1PBT" ~ "Atrophy classifier transcripts (ct>1Prob)",
                            PBT == "txbxCorrCLAGdn_tbb" ~ "Transcripts downregulated in CLAD",
                            PBT == "txbxCorrCLAGup_tbb" ~ "Transcripts upregulated in CLAD",
                            PBT == "SFT" ~ "Surfactant-associated transcripts (SFT)",
                            PBT == "ENDAT" ~ "Endothelial cell-associated transcripts (ENDAT)",
                            PBT == "eDSAST" ~ "Endothelium-expressed DSA-selective transcripts (eDSAST)",
                            PBT == "IGT" ~ "Immunoglobulin transcripts (IGT)",
                            PBT == "BAT" ~ "B cell–associated transcripts (BAT)",
                            PBT == "MCAT" ~ "Mast cell-associated transcripts (MCAT)",
                            PBT == "LT1" ~ "TBB parenchymal transcripts (LT1)",
                            PBT == "LT2" ~ "TBB parenchymal transcripts - no solute carriers (LT2)"
                        ),
                        .before = 1
                    ) %>%
                    dplyr::filter(category %nin% " ") %>%
                    arrange(order2) %>%
                    dplyr::select(-order, -order2, -PBT) %>%
                    dplyr::rename(" " = category, "  " = score)
            }
        )
    )

# flex demon table
df_res_art <- df_res_art %>%
    mutate(
        art_pairwise_table_medianIQR_flex = map2(
            txBxBin, art_pairwise_table_medianIQR_formatted,
            function(txBxBin, art_pairwise_table_medianIQR_formatted) {
                tx <- txBxBin
                N <- df_n %>%
                    dplyr::filter(txBxBin == tx) %>%
                    pull(n)
                title_pairwise <- paste(title_art_pairwise,
                    "(",
                    txBxBin,
                    " posttransplant",
                    paste("; N = ", N %>% sum(), ")", sep = ""),
                    sep = ""
                )

                art_pairwise_table_medianIQR_formatted %>%
                    flextable::flextable() %>%
                    flextable::add_header_row(top = T, values = rep(title_pairwise, ncol_keys(.))) %>%
                    flextable::add_body_row(values = c(" ", " ", paste("N = ", N, sep = "")) %>% as.list()) %>%
                    flextable::add_footer_row(values = footnoteText[[2]], colwidths = ncol_keys(.)) %>%
                    flextable::merge_v(j = 1) %>%
                    flextable::merge_at(i = 2, j = 1:2, part = "header") %>%
                    flextable::merge_at(i = 1, part = "header") %>%
                    flextable::merge_at(i = 1, j = 1:2, part = "body") %>%
                    flextable::border_remove() %>%
                    flextable::border(part = "header", border.left = fp_border(), border.top = fp_border(), border.right = fp_border()) %>%
                    flextable::border(part = "header", border.left = fp_border(), border.top = fp_border(), border.right = fp_border()) %>%
                    flextable::border(i = 1, part = "header", border.bottom = fp_border()) %>%
                    flextable::border(part = "body", border = fp_border()) %>%
                    flextable::border(part = "footer", border.left = fp_border(), border.right = fp_border()) %>%
                    flextable::border(i = 1, part = "footer", border.bottom = fp_border()) %>%
                    flextable::align(align = "center") %>%
                    flextable::align(i = 1:2, align = "center", part = "header") %>%
                    flextable::align(i = 1, j = 1:2, align = "right") %>%
                    flextable::font(fontname = "Arial", part = "all") %>%
                    flextable::fontsize(size = 10, part = "all") %>%
                    flextable::fontsize(size = 10, part = "footer") %>%
                    flextable::fontsize(i = 1, size = 12, part = "header") %>%
                    flextable::bold(part = "header") %>%
                    flextable::bold(i = 1, part = "body") %>%
                    flextable::bold(j = 1, part = "body") %>%
                    flextable::bold(i = ~ FDR < 0.05, j = 2:6, , part = "body") %>%
                    flextable::bg(i = ~ FDR < 0.05, j = 2:6, bg = "grey80", part = "body") %>%
                    flextable::bg(i = 1, bg = "white", part = "body") %>%
                    flextable::padding(padding = 0, part = "all") %>%
                    flextable::width(width = cellWidths, unit = "cm") %>%
                    flextable::width(., width = dim(.)$widths * 26.25 / (flextable_dim(.)$widths), unit = "cm")
                # flextable::height_all(., height = dim(.)$heights * 16.45 / (flextable_dim(.)$heights), unit = "cm")
            }
        )
    )

df_res_art$art_pairwise_table_medianIQR_flex[[1]]


# PRINT THE FLEXTABLES ####
df_res_art$art_pairwise_table_medianIQR_flex[[2]] %>% print(preview = "pptx")


# SAVE THE DATA FOR PLOTTING ####
pgd_pbt_anova_txbx <- df01
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP2.4 BLAD TBB/data/"
save(pgd_pbt_anova_txbx, file = paste(saveDir, "pgd_pbt_anova.RData", sep = ""))


# WRANGLE THE DATA FOR PLOTTING ####
df_test <- df01 %>%
    dplyr::select(txBxBin, PBT, art_con_tidy, medians) %>%
    unnest(art_con_tidy) %>%
    group_by(txBxBin) %>%
    nest() %>%
    mutate(
        FDR = map(
            data,
            ~ p.adjust(.x$p.value, method = "fdr")
        )
    ) %>%
    unnest(everything()) %>%
    dplyr::select(PBT, contrast, p.value, FDR, medians)

df_raw <- df01 %>%
    dplyr::select(txBxBin, data) %>%
    unnest(c(data)) %>%
    pivot_longer(
        cols = c(-PGDstatus, -txBx, -txBxBin),
        names_to = "PBT"
    ) %>%
    dplyr::select(txBxBin, PBT, PGDstatus, value)

df_plot00 <- left_join(df_raw, df_test, by = c("txBxBin", "PBT")) %>%
    mutate(
        score = case_when(
            PBT == "TCMRtPBT" ~ "TCMRProb",
            PBT == "tgt1PBT" ~ "t>1Prob",
            PBT == "igt1PBT" ~ "i>1Prob",
            PBT == "ABMRpmPBT" ~ "ABMRProb",
            PBT == "ggt0PBT" ~ "g>0Prob",
            PBT == "ptcgt0PBT" ~ "ptc>0Prob",
            PBT == "cigt1PBT" ~ "ci>1Prob",
            PBT == "ctgt1PBT" ~ "ct>1Prob",
            TRUE ~ PBT
        ),
        group = paste(PBT, PGDstatus, sep = ":"),
        category =
            case_when(
                PBT %in% Rejectionrelated ~ "Rejection-related",
                PBT %in% TCMRrelated ~ "TCMR-related",
                PBT %in% ABMRrelated ~ "ABMR-related",
                PBT %in% SFT ~ "Surfactant-related",
                PBT %in% Endothelium ~ "Endothelium-related",
                PBT %in% Parenchyma ~ "Parenchyma-related",
                PBT %in% Macrophage ~ "Macrophage-related",
                PBT %in% Injurylate ~ "Late injury-related",
                PBT %in% Injuryrecent ~ "Recent injury-related",
                TRUE ~ " "
            )
    ) %>%
    group_by(txBxBin, category) %>%
    nest() %>%
    tibble()

df_plot01 <- df_plot00 %>%
    mutate(
        ymax = map_dbl(
            data,
            function(data) {
                data %>%
                    group_by(PBT) %>%
                    summarize(ymax = quantile(value, Probs = 0.75) + (IQR(value) * 2)) %>%
                    pull(ymax) %>%
                    max()
            }
        ),
        ymin = map_dbl(
            data,
            function(data) {
                data %>%
                    group_by(PBT) %>%
                    summarize(ymin = quantile(value, Probs = 0.25) - (IQR(value) * 1.85)) %>%
                    pull(ymin) %>%
                    min()
            }
        ),
        ymin =
            case_when(
                category == "ABMR-related" ~ 0.0,
                category == "TCMR-related" ~ 0.4,
                category == "Parenchyma-related" ~ 0.5,
                category == "Late injury-related" ~ 0,
                category == "Recent injury-related" ~ 0.0,
                category == "Surfactant-related" ~ 0.5,
                TRUE ~ ymin
            )
    )


# MAKE VIOLIN PLOTS ####
df_plot02 <- df_plot01 %>%
    mutate(
        gg_violin = pmap(
            list(category, data, ymin, ymax),
            function(category, data, ymin, ymax) {
                df_FDR <- data %>%
                    group_by(score) %>%
                    summarize(
                        FDR = FDR %>% unique(),
                        p.value = p.value %>% unique()
                    ) %>%
                    mutate(
                        sig = case_when(
                            FDR < 0.0001 ~ "***",
                            FDR < 0.001 ~ "**",
                            FDR < 0.05 ~ "*",
                            p.value < 0.05 ~ ".",
                            TRUE ~ ""
                        ),
                        sig_size = ifelse(sig == ".", 5, 2)
                    )
                data %>%
                    mutate(category = category) %>%
                    ggplot(aes(x = score, y = value, fill = PGDstatus)) +
                    # geom_boxplot(outlier.alpha = 0) +
                    geom_violin(scale = "width", adjust = 1.5, alpha = 0.75, draw_quantiles = 0.499) +
                    annotate("text", x = df_FDR$score, y = ymax * 0.90, label = df_FDR$sig, vjust = 1, size = df_FDR$sig_size) +
                    scale_fill_manual(values = c("white", "grey80")) +
                    coord_cartesian(ylim = c(NA, NA)) +
                    # coord_trans(y = log10zero) +
                    labs(fill = NULL, y = NULL) +
                    theme_bw() +
                    theme(
                        legend.position = "none",
                        panel.grid = element_blank(),
                        axis.text.y = element_text(size = 6),
                        axis.text.x = element_text(
                            angle = 90, hjust = 1,
                            vjust = 0.5, size = 8,
                            colour = "black"
                        ),
                        axis.title.x = element_blank(),
                    ) +
                    facet_grid(~category)
            }
        ),
        gg_violin_legend = pmap(
            list(category, data, ymin, ymax),
            function(category, data, ymin, ymax) {
                df_FDR <- data %>%
                    group_by(score) %>%
                    summarize(
                        FDR = FDR %>% unique(),
                        p.value = p.value %>% unique()
                    ) %>%
                    mutate(
                        sig = case_when(
                            FDR < 0.0001 ~ "***",
                            FDR < 0.001 ~ "**",
                            FDR < 0.05 ~ "*",
                            p.value < 0.05 ~ ".",
                            TRUE ~ ""
                        )
                    )
                data %>%
                    mutate(
                        category = category,
                        PGDstatus = PGDstatus %>% factor(labels = c("No PGD", "PGD"))
                    ) %>%
                    ggplot(aes(x = score, y = value, fill = PGDstatus)) +
                    geom_violin(scale = "width", adjust = 1.5, alpha = 0.75, draw_quantiles = 0.499) +
                    annotate("text", x = df_FDR$score, y = ymax * 0.95, label = df_FDR$sig) +
                    scale_fill_manual(values = c("white", "grey80")) +
                    coord_cartesian(ylim = c(ymin, ymax)) +
                    labs(fill = NULL, y = NULL) +
                    theme_bw() +
                    theme(
                        legend.position = "top",
                        panel.grid = element_blank(),
                        axis.text.y = element_text(size = 6),
                        axis.text.x = element_text(
                            angle = 90, hjust = 1,
                            vjust = 0.5, size = 8,
                            colour = "black"
                        ),
                        axis.title.x = element_blank(),
                    ) +
                    facet_grid(~category)
                # ggpubr::get_legend()
            }
        )
    )


# MAKE BOXPLOTS PLOTS ####
df_plot03 <- df_plot02 %>%
    mutate(
        gg_boxplot = pmap(
            list(category, data, ymin, ymax),
            function(category, data, ymin, ymax) {
                df_FDR <- data %>%
                    group_by(score) %>%
                    summarize(
                        FDR = FDR %>% unique(),
                        p.value = p.value %>% unique()
                    ) %>%
                    mutate(
                        sig = case_when(
                            FDR < 0.0001 ~ "***",
                            FDR < 0.001 ~ "**",
                            FDR < 0.05 ~ "*",
                            p.value < 0.05 ~ ".",
                            TRUE ~ ""
                        ),
                        sig_size = ifelse(sig == ".", 5, 4)
                    )
                df_ylims <- data %>%
                    group_by(PBT) %>%
                    nest() %>%
                    tibble() %>%
                    mutate(
                        boxplot_stats = map(data, ~ boxplot.stats(.x$value)$stats),
                        ylims = map(boxplot_stats, ~ tibble(ymin = .x[1], ymax = .x[5]))
                    ) %>%
                    dplyr::select(ylims) %>%
                    unnest(everything())
                ylims <- c(
                    df_ylims %>% pull(ymin) %>% min(),
                    df_ylims %>% pull(ymax) %>% max() * 1.5
                )
                data %>%
                    mutate(category = category) %>%
                    ggplot(aes(x = score, y = value, fill = PGDstatus)) +
                    geom_boxplot(
                        position = position_dodge(0.75),
                        width = 0.75,
                        col = "grey15",
                        outlier.alpha = 0
                    ) +
                    annotate("text", x = df_FDR$score, y = ylims[2], label = df_FDR$sig, vjust = 1, size = 5) +
                    scale_fill_manual(values = c("white", "grey80")) +
                    scale_y_continuous(breaks = c(0.1, 0.3, 0.5, 1, 2, 3)) +
                    # coord_cartesian(ylim = ylims) +
                    coord_trans(y = log10zero, ylim = ylims) +
                    labs(fill = NULL, y = NULL) +
                    theme_bw() +
                    theme(
                        legend.position = "none",
                        panel.grid = element_blank(),
                        axis.text.y = element_text(size = 6),
                        axis.text.x = element_text(
                            angle = 90, hjust = 1,
                            vjust = 0.5, size = 8,
                            colour = "black"
                        ),
                        axis.title.x = element_blank(),
                    ) +
                    facet_grid(~category)
            }
        ),
        gg_boxplot_legend = pmap(
            list(category, data, ymin, ymax),
            function(category, data, ymin, ymax) {
                df_FDR <- data %>%
                    group_by(score) %>%
                    summarize(
                        FDR = FDR %>% unique(),
                        p.value = p.value %>% unique()
                    ) %>%
                    mutate(
                        sig = case_when(
                            FDR < 0.0001 ~ "***",
                            FDR < 0.001 ~ "**",
                            FDR < 0.05 ~ "*",
                            p.value < 0.05 ~ ".",
                            TRUE ~ ""
                        )
                    )
                data %>%
                    mutate(
                        category = category,
                        PGDstatus = PGDstatus %>% factor(labels = c("Post No PGD", "Post PGD"))
                    ) %>%
                    ggplot(aes(x = score, y = value, fill = PGDstatus)) +
                    geom_violin(scale = "width", adjust = 1.5, alpha = 0.75, draw_quantiles = 0.499) +
                    annotate("text", x = df_FDR$score, y = ymax * 0.95, label = df_FDR$sig) +
                    scale_fill_manual(values = c("white", "grey80")) +
                    coord_cartesian(ylim = c(ymin, ymax)) +
                    labs(fill = NULL, y = NULL) +
                    theme_bw() +
                    theme(
                        legend.position = "top",
                        panel.grid = element_blank(),
                        axis.text.y = element_text(size = 6),
                        axis.text.x = element_text(
                            angle = 90, hjust = 1,
                            vjust = 0.5, size = 8,
                            colour = "black"
                        ),
                        axis.title.x = element_blank(),
                    ) +
                    facet_grid(~category)
                # ggpubr::get_legend()
            }
        )
    )


# SIMPLE JOINT PLOTS ####
plot_joint_violin <- ggarrange(
    plotlist = df_plot03 %>%
        dplyr::filter(txBxBin == "<1-year") %>%
        pull(gg_violin),
    common.legend = TRUE
) %>%
    annotate_figure(left = "PBT score (normalized by No PGD)")

plot_joint_boxplot <- ggarrange(
    plotlist = df_plot03 %>%
        dplyr::filter(txBxBin == "<1-year") %>%
        pull(gg_boxplot),
    common.legend = TRUE
) %>%
    annotate_figure(left = "PBT score (normalized by No PGD)")



# CUSTOM JOINT PLOTS ####
df_plot04 <- df_plot03 %>%
    group_by(txBxBin) %>%
    nest() %>%
    mutate(
        plot_row1 = pmap(
            list(data),
            function(data) {
                ggarrange(
                    plotlist = c(
                        data %>%
                            dplyr::filter(category == "Rejection-related") %>%
                            pull(gg_boxplot),
                        data %>%
                            dplyr::filter(category == "ABMR-related") %>%
                            pull(gg_boxplot),
                        data %>%
                            dplyr::filter(category == "TCMR-related") %>%
                            pull(gg_boxplot)
                    ), nrow = 1,
                    align = "hv"
                )
            }
        ),
        plot_row2 = pmap(
            list(data),
            function(data) {
                ggarrange(
                    plotlist = c(
                        data %>%
                            dplyr::filter(category == "Late injury-related") %>%
                            pull(gg_boxplot),
                        data %>%
                            dplyr::filter(category == "Recent injury-related") %>%
                            pull(gg_boxplot)
                    ),
                    align = "hv"
                )
            }
        ),
        plot_row3 = pmap(
            list(data),
            function(data) {
                ggarrange(
                    plotlist = c(
                        data %>%
                            dplyr::filter(category == "Macrophage-related") %>%
                            pull(gg_boxplot),
                        data %>%
                            dplyr::filter(category == "Parenchyma-related") %>%
                            pull(gg_boxplot),
                        data %>%
                            dplyr::filter(category == "Endothelium-related") %>%
                            pull(gg_boxplot),
                        data %>%
                            dplyr::filter(category == "Surfactant-related") %>%
                            pull(gg_boxplot)
                    ), nrow = 1,
                    align = "hv"
                )
            }
        ),
        plot_legend = map(
            data,
            ~ .x$gg_boxplot_legend[[1]] %>%
                get_legend()
        ),
        plot_panel = pmap(
            list(plot_row1, plot_row2, plot_row3, plot_legend),
            function(plot_row1, plot_row2, plot_row3, plot_legend) {
                plot_joint_boxplot_custom <- ggarrange(
                    plot_legend,
                    plot_row1,
                    plot_row2,
                    plot_row3,
                    nrow = 4,
                    heights = c(0.1, 1, 1, 1, 0.5)
                ) %>%
                    annotate_figure(left = "PBT score (normalized to No PGD)")
            }
        )
    )




df_plot04$plot_panel[[1]]


df_plot03$txBxBin %>% unique()


# SAVE THE PLOTS ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP2.4 BLAD TBB/output/"

ggsave(
    plot = df_plot04 %>% dplyr::filter(txBxBin == "all") %>% pull(plot_panel) %>% pluck(1),
    filename = paste(saveDir, "PGDvNoPGD PBT ", " [boxplot].png", sep = ""),
    width = 20,
    height = 18,
    device = "png",
    units = "cm",
    dpi = 300,
    bg = "white"
)

ggsave(
    plot = df_plot04 %>% dplyr::filter(txBxBin == "<1-year") %>% pull(plot_panel) %>% pluck(1),
    filename = paste(saveDir, "PGDvNoPGD PBT lt1yr", " [boxplot].png", sep = ""),
    width = 20,
    height = 18,
    device = "png",
    units = "cm",
    dpi = 300,
    bg = "white"
)

ggsave(
    plot = df_plot04 %>% dplyr::filter(txBxBin == ">=1-year") %>% pull(plot_panel) %>% pluck(1),
    filename = paste(saveDir, "PGDvNoPGD PBT gt1yr", " [boxplot].png", sep = ""),
    width = 20,
    height = 18,
    device = "png",
    units = "cm",
    dpi = 300,
    bg = "white"
)
