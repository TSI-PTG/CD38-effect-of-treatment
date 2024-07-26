# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
# Custom operators, functions, and datasets
"%nin%" <- function(a, b) match(a, b, nomatch = 0) == 0
# load existing PBT list
load("Z:/DATA/Datalocks/Other data/PBTlist219_14Dec2020_JR.RData")
# load new PBT list
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/PBTlist219_ABMRactivity_26Jul2024_PTG.RData")
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/PBTlist219_AKIinduced_26Jul2024_PTG.RData")


# COMBINED PBT LISTS ####
PBTlist219_new <- c(PBTlist219, PBTlist219_ABMRactivity, PBTlist219_AKIinduced)


# SAVE AN UPDATED PBTkust219 ####
PBTlist219 <- PBTlist219_new
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/"
save(PBTlist219, file = paste(saveDir, "PBTlist219_26Jul2024_PTG.RData", sep = ""))
saveDir2 <- "Z:/DATA/Datalocks/Other data/"
save(PBTlist219, file = paste(saveDir2, "PBTlist219_26Jul2024_PTG.RData", sep = ""))
