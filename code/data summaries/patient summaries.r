# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(flextable) # install.packages("flextable") #for table outputs
library(officer) # install.packages("officer")
# Bioconductor libraries
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
# source plot function
source("C:/R/CD38-effect-of-treatment/code/functions/plot.gg_violin_interaction.r")
source("C:/R/CD38-effect-of-treatment/code/functions/plot.gg_patient_pairs_interaction.r")
# load reference set
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/data_felzartamab_k1208.RData")



# DEFINE SET ####
data_K1208 <- data_felzartamab_k1208


# PATIENT SUMMARIES ####
patient_summary <- data_K1208 %>%
    dplyr::select(CEL, Center, Patient, Felzartamab, Group, Followup) %>%
    dplyr::filter(CEL %nin% c("FBN003_NBN010_B2_(PrimeView).CEL", "FVI022_FVI022_B2_(PrimeView).CEL")) %>%
    dplyr::select(-Group) %>%
    pivot_wider(names_from = Followup, values_from = CEL) %>%
    arrange(Patient, Felzartamab)

patient_summary %>%
    mutate(
        across(
            c("Day0", "Week12", "Week24", "Week52"),
            ~ ifelse(. %>% is.na() | . == "", NA, "X")
        )
    ) %>%
    flextable::flextable() %>%
    flextable::add_header_row(values = rep("Felzartamab study population", ncol_keys(.))) %>%
    flextable::merge_h(part = "header") %>%
    flextable::merge_v(j = 1:2) %>%
    flextable::fontsize(size = 8, part = "all") %>%
    flextable::align(align = "center", part = "all") %>%
    flextable::bg(bg = "white", part = "all") %>%
    flextable::colformat_double(j = 2:3, digits = 2) %>%
    flextable::border_remove() %>%
    flextable::bold(part = "header") %>%
    flextable::padding(padding = 0, part = "all") %>%
    flextable::border(border = fp_border(), part = "all") %>%
    flextable::autofit()
# %>%
# print(preview = "pptx")

# table with cel id's
patient_summary %>%
    flextable::flextable() %>%
    flextable::add_header_row(values = rep("Felzartamab study population", ncol_keys(.))) %>%
    flextable::merge_h(part = "header") %>%
    flextable::merge_v(j = 1:2) %>%
    flextable::fontsize(size = 8, part = "all") %>%
    flextable::align(align = "center", part = "all") %>%
    flextable::bg(bg = "white", part = "all") %>%
    flextable::colformat_double(j = 2:3, digits = 2) %>%
    flextable::border_remove() %>%
    flextable::bold(part = "header") %>%
    flextable::padding(padding = 0, part = "all") %>%
    flextable::border(border = fp_border(), part = "all") %>%
    flextable::autofit()

data_K1208 %>%
    dplyr::filter(Patient %nin% c(15, 18)) %>%
    dplyr::select(Felzartamab, Followup) %>%
    table() %>%
    as_tibble() %>%
    pivot_wider(names_from = Followup, values_from = n)
