# HOUSEKEEPING ####
# CRAN packages
library(tidyverse) # install.packages("tidyverse")
# Custom operators, functions, and datasets
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
source("C:/R/CD38-effect-of-treatment/code/functions/get_slope_function.r")
# Custom operators, functions, and datasets
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
se <- function(x) sd(x) / sqrt(length((x)))
log10zero <- scales::trans_new(
    name = "log10zero",
    transform = function(x) log10(x + 0.01),
    inverse = function(x) 10^x - 0.01
)
corPvalueStudent <- function(cor, nSamples) {
    T <- sqrt(nSamples - 2) * cor / sqrt(1 - cor^2)
    2 * pt(abs(T), nSamples - 2, lower.tail = FALSE)
}
# load reference data
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/data_felzartamab_k1208.RData")


# DEFINE SEED ####
seed <- 314159


# DEFINE CATEGORIES FOR FEATURES ####
# Rejectionrelated <- c("GRIT3", "Rej-RAT", "RejAA_NR")
vars_cfDNA <- c("cfDNA")
vars_abmr <- c("ABMRpm", "ggt0", "ptcgt0", "DSAST", "NKB") # "cggt0", "RejAA_EABMR", "RejAA_FABMR", "RejAA_LABMR")
vars_tcmr <- c("TCMRt", "tgt1", "igt1", "QCAT", "TCB") # , "TCMR-RAT", )
vars_injury <- c("IRRAT30", "IRITD3", "IRITD5")
# Endothelium <- c("ENDAT")
# Parenchyma <- c("KT1", "KT2")
# Macrophage <- c("AMAT1", "QCMAT")
# Injuryrecent <- c("FICOL", "IRRAT30", "IRITD3", "IRITD5")
# Injurylate <- c("IGT", "MCAT", "BAT", "cigt1", "ctgt1")


# DEFINE VARIABLES TO ASSESS ####
# vars <- c("cfDNA", "ABMRpm", "ggt0", "ptcgt0", "DSAST", "NKB", "TCMRt", "tgt1", "igt1", "QCAT", "TCB")
vars <- c(vars_cfDNA, vars_abmr, vars_tcmr, vars_injury)


# DEFINE THE DATA ####
data <- data_felzartamab_k1208 %>%
    dplyr::select(Center, Patient, Felzartamab, Group, Followup, Felzartamab_Group, Felzartamab_Followup, all_of(vars)) %>%
    dplyr::filter(Group != "FU1b", Patient %nin% c(15, 18)) %>%
    left_join(., summarise(., sample_pairs = n(), .by = Felzartamab_Followup), by = "Felzartamab_Followup") %>%
    relocate(sample_pairs, .after = Felzartamab_Followup) %>%
    arrange(Felzartamab, Patient, Group)


# WRANGLE THE PHENOTYPE DATA ####
data_00 <- data %>%
    expand_grid(category = c("cfDNA", "ABMR", "TCMR", "injury")) %>%
    nest(.by = category) %>%
    mutate(
        features = map(
            category,
            function(category) {
                if (category == "cfDNA") {
                    vars_cfDNA
                } else if (category == "ABMR") {
                    vars_abmr
                } else if (category == "TCMR") {
                    vars_tcmr
                } else if (category == "injury") {
                    vars_injury
                }
            }
        ),
        data = pmap(
            list(features, data),
            function(features, data) {
                data %>%
                    dplyr::select(
                        Center, Patient, Felzartamab, Group, Followup, Felzartamab_Group, Felzartamab_Followup, sample_pairs,
                        all_of(features)
                    )
            }
        )
    )


# WRANGLE THE DATA FOR UNIVARIATE TESTS ####
data_01 <- data_00 %>%
    mutate(
        data_univariate = pmap(
            list(data, features),
            function(data, features) {
                data %>%
                    pivot_longer(
                        cols = all_of(features),
                        names_to = "variable",
                        values_to = "value"
                    ) %>%
                    nest(.by = variable) %>%
                    mutate(
                        variable = variable %>%
                            factor(
                                levels = vars
                            ),
                        annotation = case_when(
                            variable %in% vars_cfDNA ~ "cfDNA",
                            variable %in% vars_tcmr ~ "TCMR-related",
                            variable %in% vars_abmr ~ "ABMR-related",
                            variable %in% vars_injury ~ "injury-related",
                            TRUE ~ " "
                        ) %>%
                            factor(
                                levels = c(
                                    "cfDNA",
                                    "ABMR-related",
                                    "TCMR-related",
                                    "injury-related"
                                )
                            ),
                        score = case_when(
                            variable == "cfDNA" ~ "Donor-derived cell-free DNA (dd-cfDNA, cp/mL)",
                            variable == "TCMRt" ~ "TCMR classifier (TCMRProb)",
                            variable == "TCMRt" ~ "TCMR classifier (TCMRProb)",
                            variable == "TCB" ~ "T-cell burden (TCB)",
                            variable == "TCMR-RAT" ~ "TCMR-associated RATs (TCMR-RAT)",
                            variable == "tgt1" ~ "Tubulitis classifier (t>1Prob)",
                            variable == "igt1" ~ "Interstitial infiltrate classifier (i>1Prob)",
                            # variable == "QCAT" ~ "Cytotoxic T cell-associated transcripts (QCAT)",
                            variable == "QCAT" ~ "Cytotoxic T cell-associated (QCAT)",
                            variable == "Rej-RAT" ~ "Rejection-associated transcripts (Rej-RAT)",
                            variable == "GRIT3" ~ "Interferon gamma-inducible transcripts (GRIT3)",
                            variable == "ABMR-RAT" ~ "ABMR-associated RATs (ABMR-RAT)",
                            variable == "ABMRpm" ~ "ABMR classifier (ABMRProb)",
                            variable == "DSAST" ~ "DSA-selective transcripts (DSAST)",
                            variable == "NKB" ~ "NK cell burden (NKB)",
                            variable == "ggt0" ~ "Glomerulitis classifier (g>0Prob)",
                            variable == "cggt0" ~ "Glomerular double contours classifier (cg>0Prob)",
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
                            variable == "KT2" ~ "Kindey parenchymal transcripts - no solute carriers (KT2)",
                            variable == "RejAA_NR" ~ "Archetypal No Rejection score (NR)",
                            variable == "RejAA_EABMR" ~ "Archetypal Early ABMR score (EABMR)",
                            variable == "RejAA_FABMR" ~ "Archetypal Full ABMR score (FABMR)",
                            variable == "RejAA_LABMR" ~ "Archetypal Late ABMR score (LABMR)",
                        ), .before = 1
                    ) %>%
                    arrange(annotation, variable)
            }
        )
    ) %>%
    dplyr::select(category, data_univariate) %>%
    unnest(data_univariate) %>%
    mutate(n_cat = score %>% unique() %>% length(), .by = category, .after = variable)


# UNIVARIATE MEDIANS ####
data_02 <- data_01 %>%
    mutate(
        delta = map(
            data,
            function(data) {
                data %>%
                    reframe(
                        delta = combn(value, 2, diff) %>% as.numeric(),
                        .by = c(Patient, Felzartamab)
                    ) %>%
                    mutate(
                        Followup_pairwise = rep(c("Baseline - Week24", "Baseline - Week52", "Week24 - Week52"), 20) %>%
                            factor(levels = c("Baseline - Week24", "Baseline - Week52", "Week24 - Week52")),
                        .after = Felzartamab
                    ) %>%
                    arrange(Patient, Felzartamab)
            }
        )
    )
data_02$data[[1]]
data_02$delta[[1]]



# REFRAME TO ADD DELTA CFDNA TO INDIVIDUAL COLUMN ####
data_03 <- data_02 %>%
    mutate(
        delta_cfdna = data_02 %>%
            dplyr::filter(category == "cfDNA") %>%
            pull(delta)
    ) %>%
    dplyr::filter(category != "cfDNA")



# UNIVARIATE NONPARAMETRIC TESTS ####
data_04 <- data_03 %>%
    mutate(
        cor = pmap(
            list(variable, delta, delta_cfdna),
            function(variable, delta, delta_cfdna) {
                df <- delta %>%
                    left_join(
                        delta_cfdna,
                        by = c("Patient", "Felzartamab", "Followup_pairwise"),
                        suffix = c("_score", "_cfdna")
                    ) %>%
                    mutate(variable) %>% 
                    nest(.by = c("variable", "Felzartamab", "Followup_pairwise"))
                df %>%
                    mutate(
                        cor = pmap_dbl(
                            list(data),
                            function(data) {
                                cor(x = data$delta_cfdna, y = data$delta_score, use = "p", method = "spearman")
                            }
                        ),
                        p = pmap_dbl(
                            list(data, cor),
                            function(data, cor) {
                                nSamples <- data %>%
                                    nrow()
                                corPvalueStudent(cor, nSamples)
                            }
                        )
                    ) 
            }
        )
    )
data_04$cor



# SAVE DATA NEST ####
felzartamab_cfdna_cor_k1208 <- data_04
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/"
save(felzartamab_cfdna_cor_k1208, file = paste(saveDir, "felzartamab_cfdna_cor_k1208", sep = ""))



# END ####
