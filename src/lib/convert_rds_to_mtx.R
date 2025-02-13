library(HDF5Array)
library(Matrix)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
input_path <- args[1]

# Extract base name (e.g., "GSE154659_C57_Raw_counts") from input file name
base_name <- tools::file_path_sans_ext(basename(input_path))
base_name <- gsub("_Raw_counts$", "", base_name)

# Define output paths based on the base name
base_dir <- dirname(input_path)
output_mtx <- file.path(base_dir, paste0(base_name, "_Raw_counts.mtx"))
output_genes <- file.path(base_dir, paste0(base_name, "_genes.csv"))
output_barcodes <- file.path(base_dir, paste0(base_name, "_barcodes.csv"))

# Load the data from the uncompressed .RDS file
data <- readRDS(input_path)
print(class(data))

# Convert to dataframe
df <- as.data.frame(as.matrix(data))

# Extract row & column names
col_names <- colnames(df)
row_names <- rownames(df)

if (!inherits(data, "sparseMatrix")) {
    data <- as(data, "sparseMatrix")
}

# Write matrix and metadata files
writeMM(data, file = output_mtx)
write.csv(row_names, file = output_genes, row.names = FALSE, col.names = FALSE)
write.csv(col_names, file = output_barcodes, row.names = FALSE, col.names = FALSE)
