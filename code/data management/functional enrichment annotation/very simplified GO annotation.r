# IMMUNE PROCESSES ####
immune_response <- paste(
    c(
        "immune", "immunity", "cytokine", "leukocyte", "cell activation",
        "interaction", "virus", "symbiont", "defense response",
        "cell killing", "cell surface receptor", "interspecies"
    ),
    collapse = "|"
)
infection_response <- paste(c("virus", "symbiont", "defense response"), collapse = "|")


# STRESS RESPONSE ####
stress_response <- paste(c("cellular response to stress", "stress"), collapse = "|")


# FOREIGN STIMULUS ####
exogenous_stimlulus <- paste(c("external stimulus", "interaction"), collapse = "|")


# ENDOGENOUS STIMULUS ####
endogenous_stimulus <- paste(c("endogenous stimulus", "cellular response to endogenous stimulus"), collapse = "|")


# INJURY PROCESSES ####
inflammation <- paste(c("inflam"), collapse = "|")
injury <- paste(c("injury"), collapse = "|")


# CELLULAR PROCESSES ####
cell_cycle <- paste(c("cycle"), collapse = "|")
cell_signalling <- paste(c("communication", "signal", "signalling"), collapse = "|")
cell_mobilization <- paste(c("locomotion", "loco", "migration", "motility"), collapse = "|")
cellular_development <- paste(c(
    "chromosome", "organelle fission", "organization", "segregation", "division",
    "development", "neurogenesis", "generation", "morphogenesis", "differentiation", "component",
    "biological_process", "hemopoiesis", "localization", "endocytosis"
), collapse = "|")
cellular_regulation <- paste(c("cellular process", "regulation of cellular", "positive regulation of response to stimulus"), collapse = "|")


# METABOLIC PROCESSES ####
protein_metabolism <- paste(
    c("protein metabolic", "protein metabolism"),
    collapse = "|"
)

nitrogen_metabolism <- paste(
    c("nitrogen", "response to organonitrogen compound"),
    collapse = "|"
)

xenobiotic_metabolism <- paste(
    c("organic substance", "macromolecule metabolic", "response to organic cyclic compound"),
    collapse = "|"
)

metabolic_response <- paste(c("metabolism", "metabolic", "catabolic", "cellular metabolic"), collapse = "|")
general_metabolic_response <- paste(protein_metabolism, nitrogen_metabolism, xenobiotic_metabolism, metabolic_response,  sep = "|")


# PROTEIN SYNTHESIS ####
protein_synthesis <- paste(
    c("positive regulation of macromolecule biosynthetic process"),
    collapse = "|"
)


# MERGE ANNOTATIONS ####
cellular_development_and_metabolism <- (general_metabolic_response, cellular_development, cellular_regulation)



# DEFINE THE GO ANNOTATION LEVELS ####
GO_annotation_levels_truncated <- c(
    "immune response", "response to infection","stress response",
    "response to exogenous stimulus", "response to endogenous stimulus",
    "inflammation", "injury response",
    "cell cycling", "cell signalling", "cell mobilization", "cell development", "cellular regulation",
    "metabolic response",
    "protein_synthesis"
)
