# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(openxlsx) # install.packages("openxlsx")
# Custom operators, functions, and datasets
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
# load gene lists
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/genes_NK_GEP.RData")
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/ABMR_NK_genes.RData")
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/ABMR_endothelial_genes.RData")
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/ABMR_activity_genes.RData")
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/ABMR_IFNG_genes.RData")


# WARNGLE THE GENESETS ####
ABMRactivity_genesets <- tibble(
    name = c(
        "ABMR activity genes",
        "IFNG-inducible ABMR activity genes",
        "NK cell-expressed ABMR activity genes",
        "ABMR-associated endothelial genes"
    ),
    abbreviation = c("AAG", "IIAAG", "NKAAG", "AEG"),
    geneset = list(
        AAG = genes_ABMR_activity,
        IIAAG = genes_ABMR_IFNG,
        NKAAG = genes_ABMR_NK,
        AEG = genes_ABMR_endothelial
    )
)


# EXPORT THE DATA AS AN EXCEL SHEET ####
saveDir1 <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/output/"
openxlsx::write.xlsx(ABMRactivity_genesets$geneset,
    asTable = TRUE,
    file = paste(saveDir1, "ABMR activity genesets",
        # Sys.Date(),
        # format(Sys.time(), "_%I%M%p"),
        ".xlsx",
        sep = ""
    )
)
