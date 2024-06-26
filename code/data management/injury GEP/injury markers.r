# HOUSEKEEPING ####
# CRAN libraries
library(tidyverse) # install.packages("tidyverse")
library(readxl) # install.packages("readxl")
library(circlize) # install.packages("circlize")
# load AFFYMAP
load("Z:/DATA/Datalocks/Other data/affymap219_21Oct2019_1306_JR.RData") # for labeling genes
# load hinze single cell results normal vs AKI
single_cell_path <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/Hinze single cell state injury markers.xlsx"
sheetnames <- excel_sheets(single_cell_path)
mylist <- lapply(sheetnames, read_excel, path = single_cell_path, col_types = c(rep("text", 4), rep("numeric", 5)))
names(mylist) <- sheetnames
hinze <- mylist %>% tibble(celltype = names(.), data = .)
# load mean expression in K1208
load("Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/mean_expression_K1208_MMDx.RData")


# MEAN EXPRESSION BY PROBE ####
mean_exprs_by_probe <- means_K1208 %>%
    mutate(
        mean_exprs = means_K1208 %>%
            dplyr::select(-AffyID:-PBT) %>%
            rowMeans()
    ) %>%
    dplyr::slice_max(mean_exprs, by = Symb)

affy_genes <- mean_exprs_by_probe %>%
    dplyr::filter(Symb != "") %>%
    distinct(Symb, .keep_all = TRUE) %>%
    dplyr::select(AffyID, Symb, Gene, PBT)


# WRANGLE GENE EXPRESSION PROFILES ####
genes_injury_markers <- hinze %>%
    mutate(
        data = map(
            data,
            function(data) {
                data %>%
                    left_join(affy_genes, by = "Symb") %>%
                    mutate(contrast = "AKI vs control")  %>% 
                    relocate(AffyID, Symb, Gene, PBT, contrast, .before = log2FC)
            }
        )
    )
genes_injury_markers$data[[1]]


# EXPORT THE GENE SET ####
saveDir <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/data/"
save(genes_injury_markers, file = paste(saveDir, "Hinze_injury_markers.RData", sep = ""))


# EXPORT THE DATA AS AN EXCEL SHEET ####
saveDir1 <- "Z:/MISC/Phil/AA All papers in progress/A GC papers/AP1.0A CD38 molecular effects Matthias PFH/output/"
openxlsx::write.xlsx(genes_injury_markers$data,
    asTable = TRUE,
    file = paste(saveDir1, "Hinze_injury_markers_26Jun24",
        # Sys.Date(),
        # format(Sys.time(), "_%I%M%p"),
        ".xlsx",
        sep = ""
    )
)
