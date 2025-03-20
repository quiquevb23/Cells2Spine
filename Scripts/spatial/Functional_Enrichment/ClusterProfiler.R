#!/usr/bin/env Rscript

# Load necessary libraries
library(jsonlite)  # For reading JSON config
library(clusterProfiler)
library(org.Rn.eg.db)  # Genome-wide annotation for rat
library(dplyr)
library(readr)

# Load config file
config <- fromJSON("config.json")
input_dir <- config$input_dir

# Function to perform enrichment analysis
run_enrichment <- function(gene_list, category, direction, output_dir) {
  if (length(gene_list) == 0) {
    return(NULL)
  }
  
  enrich_result <- enrichGO(
    gene = gene_list,
    OrgDb = org.Rn.eg.db,
    keyType = "SYMBOL",
    ont = category,
    pvalueCutoff = 0.05
  )
  
  if (!is.null(enrich_result)) {
    output_file <- file.path(output_dir, paste0("go_enrichment_", direction, "_", category, ".csv"))
    write.csv(as.data.frame(enrich_result), output_file, row.names = FALSE)
  }
}

# Loop through directories in input_dir
dirs <- list.dirs(input_dir, recursive = FALSE)
dirs <- dirs[!grepl("Spatial", basename(dirs))]  # Exclude directories starting with "Spatial"

for (dir_path in dirs) {
  # Get list of CSV files
  csv_files <- list.files(dir_path, pattern = "pseudobulk\\.csv$", full.names = TRUE)

  for (file_path in csv_files) {
    # Read the DEGs CSV file
    deg_df <- read_csv(file_path)

    # Ensure output directory exists inside the input directory
    output_enr_dir <- file.path(dir_path, "ClusterProfilerEnrichment")
    if (!dir.exists(output_enr_dir)) {
      dir.create(output_enr_dir, recursive = TRUE)
    }

    # Filter significant DEGs (FDR < 0.05)
    filtered_degs <- deg_df %>% filter(FDR < 0.05)

    # Separate upregulated and downregulated genes
    upregulated_genes <- filtered_degs %>% filter(logFC > 0) %>% pull(genes)
    downregulated_genes <- filtered_degs %>% filter(logFC < 0) %>% pull(genes)

    # Perform enrichment analysis for each ontology and save results in the same directory
    for (category in c("BP", "CC", "MF")) {
      run_enrichment(upregulated_genes, category, "upregulated", output_enr_dir)
      run_enrichment(downregulated_genes, category, "downregulated", output_enr_dir)
    }
  }
}

cat("Functional enrichment analysis completed. Results saved in input directory.\n")
