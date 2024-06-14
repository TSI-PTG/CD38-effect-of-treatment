# PATHWAY KEYS ####
immune_response <- paste(
    c(
        "immune", "immunity", "cytokine", "leukocyte", "cell activation",
        "response to", "interaction", "virus", "symbiont", "defense response",
        "cell killing"
    ),
    collapse = "|"
)
cell_cycle <- paste(c("cycle"), collapse = "|")
inflammation <- paste(c("inflam"), collapse = "|")
injury <- paste(c("injury"), collapse = "|")
external_stimulus <- paste(c("response to", "interaction"), collapse = "|")
reg_cellular_processes <- paste(c("regulation of"), collapse = "|")
cellular_development <- paste(c(
    "chromosome", "organelle fission", "organization", "segregation", "division",
    "development", "neurogenesis", "generation", "morphogenesis", "differentiation", "component",
    "cellular process", "biological_process"
), collapse = "|")
cellular_communication <- paste(c("communication", "signal", "signalling"), collapse = "|")
infection_response <- paste(c("virus", "symbiont", "defense response"), collapse = "|")
metabolic_response <- paste(c("metabolism", "metabolic", "catabolic"), collapse = "|")