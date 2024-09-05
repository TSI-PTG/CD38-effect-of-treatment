# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(vegan) # install.packages("vegan") #for permanova
library(flextable) # install.packages("flextable")
library(officer) # install.packages("officer")
# Bioconductor libraries
library(Biobase) # BiocManager::install("Biobase")
# Custom operators and functions
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
# Suppress pesky dplyr reframe info
options(dplyr.reframe.inform = FALSE)
# load reference set
load("natmed/data/cd38_3Sept24.RData")


# DEFINE SEED ####
seed <- 42


# DEFINE CATEGORIES FOR FEATURES ####
vars_cfDNA <- c("cfDNA_cpml")
vars_abmr <- c("ABMRpm", "ggt0", "ptcgt0", "DSAST", "NKB")
vars_tcmr <- c("TCMRt", "tgt1", "igt1", "QCAT", "TCB")
vars_injury <- c("IRRAT30", "IRITD3", "IRITD5", "cigt1", "ctgt1")
vars_parenchyma <- c("KT1", "KT2")
vars_macrophage <- c("AMAT1", "QCMAT")


# DEFINE VARIABLES TO ASSESS ####
vars <- c(vars_cfDNA, vars_abmr, vars_tcmr, vars_macrophage, vars_injury, vars_parenchyma)



# DEFINE THE data ####
data <- set %>%
    pData() %>%
    dplyr::select(Center, Patient, Felzartamab, Group, Followup, Felzartamab_Group, Felzartamab_Followup, all_of(vars)) %>%
    dplyr::filter(Patient %nin% c(15, 18)) %>%
    left_join(., summarise(., sample_pairs = n(), .by = Felzartamab_Followup), by = "Felzartamab_Followup") %>%
    relocate(sample_pairs, .after = Felzartamab_Followup) %>%
    arrange(Felzartamab, Patient, Group)


# WRANGLE THE PHENOTYPE DATA ####
df00 <- data %>%
    expand_grid(category = c("cfDNA", "ABMR", "TCMR", "macrophage", "injury", "parenchyma")) %>%
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
                } else if (category == "macrophage") {
                    vars_macrophage
                } else if (category == "injury") {
                    vars_injury
                } else if (category == "parenchyma") {
                    vars_parenchyma
                }
            }
        ),
        data = pmap(
            list(features, data),
            function(features, data) {
                data %>%
                    dplyr::select(
                        Center, Patient, Felzartamab, Group, Followup, Felzartamab_Group, Felzartamab_Followup,
                        all_of(features)
                    )
            }
        )
    )


# PERMANOVA ####
permanova <- df00 %>%
    mutate(
        permanova = pmap(
            list(data, features),
            function(data, features) {
                set.seed(seed)
                features <- data %>% dplyr::select(dplyr::all_of(features))
                Followup <- data$Followup
                Felzartamab <- data$Felzartamab
                res <- vegan::adonis2(
                    features ~ Followup * Felzartamab,
                    data = data,
                    method = "euclidean",
                    permutations = 1e+06
                )
                attr(res, "F.perm") <- NULL
                return(res)
            }
        )
    )
names(permanova$permanova) <- permanova$category


# EXPORT RESULTS ####
saveDir <- "natmed/results/"
save(permanova, file = paste(saveDir, "permanova.RData", sep = ""))
