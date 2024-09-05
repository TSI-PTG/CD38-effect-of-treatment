# THIS SCRIPT PROCESSES THE RAW GEO DOWNLOAD #


# HOUSEKEEPING ####
library(tidyverse) # install.packages("tidyverse")
library(R.utils) # install.packages("R.utils")
library(BiocManager) # install.packages("BiocManager")
library(affy) # BiocManager::install("affy")
library(Biobase) # BiocManager::install("Biobase")



library(GEOquery) # BiocManager::install("GEOquery")
library(httr) # install.packages("httr")
library(curl)

gse <- getGEO("GSE275824", GSEMatrix = TRUE)


geo_accession <- "GSE275824"
token <- 'qpcneasuxdyffcj'



# Set up the URL for the GEO series
url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE275824&format=file"

# url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE275824&format=file&token=qpcneasuxdyffcj"

response <- GET(url, add_headers(Authorization = paste("Bearer", token)), config = list(verbose = TRUE))


handle <- new_handle()
handle_setheaders(handle, Authorization = paste("Bearer", token))
handle_setopt(handle, verbose = TRUE) # Enable verbose output


curl_download(url, destfile = paste(dir_download, "test_download", sep = ""))


curl_download(url, "downloaded_data.gz", handle = handle)


# DEFINE LOCAL DIRECTORIES ####
# note: at present the GEO access is private and restricted to reviewers
# pleae adjust to you local directories
dir_download <- "C:/GSE275824/"
dir_save <- "C:/GSE275824/CEL/"


# EXTRACT THE INITIAL TARBALL OF CEL FILES ####
untar(
    tarfile = paste(dir_download, "GSE275824_RAW.tar", sep = ""),
    exdir = dir_save
)


# EXTACT THE INDIVIDUALLY ZIPPED CEL FILES ####
# Get the list of all .gz files in the directory
gz_files <- list.files(dir_save, pattern = "*.gz", full.names = TRUE)
# Unzip each .gz file
lapply(gz_files, gunzip, remove = TRUE)


# EXTRACT THE PHENOTYPE DATA ####
gunzip(
    filename = paste(dir_download, "GSE275824_complete_data.txt.gz", sep = ""),
    remove = FALSE,
    overwrite = TRUE
)


# READ THE CEL FILES AND CREATED AN AffyBatch object ####
cel_files <- list.files(dir_save, pattern = "*.CEL", full.names = TRUE)
affy_data <- ReadAffy(filenames = cel_files)
affy_data
set <- ExpressionSet(assayData = affy_data  %>% exprs)
set




# READ THE PHENOTYPE DATA ####
pheno_data <- read.table(
    paste(dir_download, "GSE275824_complete_data.txt", sep = ""),
    header = TRUE,
    sep = "\t",
    quote = ""
) %>%
    as_tibble

















# NORMALIZE THE EXPRESSION SET ####
eset <- exprs(affy_data) %>%
    as_tibble(rownames = "AffyID")

eset  %>% str



eset <- ExpressionSet(
    assayData = exprs(affy_data)
    # phenoData = AnnotatedDataFrame(pheno_data)
)

