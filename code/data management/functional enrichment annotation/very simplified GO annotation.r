# IMMUNE PROCESSES ####
immune_response <- paste(
    c(
        "immune", "immunity", "cytokine", "leukocyte", "cell activation",
        "interaction", "virus", "symbiont", "defense response",
        "cell killing", "cell surface receptor", "interspecies", "other organism",
        "immunoglobulin production"
    ),
    collapse = "|"
)
infection_response <- paste(c("virus", "symbiont", "defense response"), collapse = "|")


# STRESS RESPONSE ####
stress_response <- paste(c("cellular response to stress", "stress"), collapse = "|")


# FOREIGN STIMULUS ####
exogenous_stimulus <- paste(c("external stimulus", "interaction", "external biotic stimulus", "biotic stimulus"), collapse = "|")


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
cellular_regulation <- paste(c(
    "cellular process", "regulation of cellular", "positive regulation of response to stimulus",
    "positive regulation of biological process"
), collapse = "|")


rna_transcription <- paste(c(
    "RNA biosynthetic process", "templated transcription", "RNA splicing",
    "transcription by RNA", "nucleobase-containing"
), collapse = "|")



# METABOLIC PROCESSES ####
protein_metabolism <- paste(
    c("protein metabolic", "protein metabolism", "regulation of protein modification process"),
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

organic_metabolism <- paste(
    c(
        "compound biosynthetic", "macromolecule metabolic",
        "macromolecule biosynthetic"
    ),
    collapse = "|"
)

metabolic_response <- paste(c("metabolism", "metabolic", "catabolic", "cellular metabolic", "biosynthetic process"), collapse = "|")



# PROTEIN SYNTHESIS ####
protein_synthesis <- paste(
    c("positive regulation of macromolecule biosynthetic process"),
    collapse = "|"
)


# HOMEOSTASIS ####
homeostasis <- paste(c("homeostasis"),collapse = "|")


# MERGE ANNOTATIONS ####
general_metabolic_response <- paste(protein_metabolism, nitrogen_metabolism,
    xenobiotic_metabolism, metabolic_response, organic_metabolism,
    sep = "|"
)
cellular_development_and_metabolism <- paste(general_metabolic_response, cell_cycle, cell_mobilization, cellular_development, cellular_regulation, protein_synthesis, sep = "|")
response_to_stimulus <- paste(endogenous_stimulus, exogenous_stimulus, stress_response, sep = "|")
cell_signalling_and_RNA_transcription <- paste(cell_signalling, rna_transcription, sep = "|")




# DEFINE THE GO ANNOTATION LEVELS ####
GO_annotation_levels_truncated <- c(
    "immune response", "response to infection",
    "inflammation", "injury response",
    "response to exogenous/endogenous stimulus",
    "cell signalling and RNA transcription",
    "cell development,\nmobilization,\nand metabolism", 
    "homeostasis"
)
