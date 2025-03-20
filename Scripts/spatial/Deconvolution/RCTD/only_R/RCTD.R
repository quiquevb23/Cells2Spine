#R script to run RCTD for each of spatial samples

library(renv)
renv::restore(prompt = FALSE)

library(jsonlite)
library(spacexr)
library(Matrix)

# Load config values
config <- fromJSON("config.json")
tmp_dir <- config$tmp_directory # Read tmp_dir 
input_dir <- config$input_directory #Read input_dir to save results RDS

# Function to load reference data
load_reference_data <- function(label, tmp_dir) {
    label_dir <- file.path(tmp_dir, label)

    counts_file <- file.path(label_dir, paste0("counts.mtx"))
    cell_types_file <- file.path(label_dir, paste0("cell_types.csv"))
    genes_file <- file.path(label_dir, paste0("genes.csv"))
    barcodes_file <-  file.path(label_dir, paste0("cells.csv"))
    
    # Load row (genes) and column (cell) names, skipping the first line
    genes <- read.csv(genes_file, header=FALSE, stringsAsFactors=FALSE)[,1]
    barcodes <- read.csv(barcodes_file, header=FALSE,	stringsAsFactors=FALSE)[,1]

    # Load counts matrix 
    counts <- readMM(counts_file)
    
    counts <- as(counts, "CsparseMatrix")  # Convert to sparse format
    
    # Ensure correct dimensions before assigning names
    if (length(genes) != nrow(counts)) {
        stop("Error: Number of genes does not match row count in counts matrix!")
    }
    if (length(barcodes) != ncol(counts)) {
        stop("Error: Number of barcodes does not match column count in counts matrix!")
    }
    
    # Handle duplicate gene names
    if (any(duplicated(genes))) {
        warning("Duplicate gene names detected! Making names unique.")
        genes <- make.unique(genes)  # Appends ".1", ".2", etc. to duplicates
    }

    # Handle duplicate gene names
    if (any(duplicated(barcodes))) {
        warning("Duplicate gene names detected! Making names unique.")
        barcodes <- make.unique(barcodes)  # Appends ".1", ".2", etc. to duplicates
    }

    # Assign row and column names
    rownames(counts) <- genes
    colnames(counts) <- barcodes

    # Load cell types metadata (No header → Manually define column names)
    cell_types_df <- read.csv(cell_types_file, header=FALSE, stringsAsFactors=FALSE)
    colnames(cell_types_df) <- c("barcodes", "cell_type")  # Manually assign names

    # Ensure unique barcode names in metadata
    if (any(duplicated(cell_types_df$barcodes))) {
        warning("Duplicate barcodes in metadata detected! Making them unique.")
        cell_types_df$barcodes <- make.unique(cell_types_df$barcodes)
    }

    # Ensure barcodes match the counts matrix
    cell_types_df <- cell_types_df[cell_types_df$barcodes %in% colnames(counts), ]
    # Convert to factor
    cell_types <- as.factor(cell_types_df$cell_type)
    names(cell_types) <- cell_types_df$barcodes

    #cell_types <- cell_types_df$cell_type
    #names(cell_types) <- cell_types_df$barcodes  # Match barcodes to colnames(counts)
    #cell_types <- as.factor(cell_types)
    # Print number of barcodes per cell type

    cat("Cell Type Counts:\n")
    cell_type_counts <- table(cell_types)
    print(cell_type_counts)

    # Remove cell types with fewer than 25 barcodes
    valid_cell_types <- names(cell_type_counts[cell_type_counts >= 25])
    cell_types_df <- cell_types_df[cell_types_df$cell_type %in% valid_cell_types, ]
    
    # Update counts matrix and cell types
    counts <- counts[, cell_types_df$barcodes, drop=FALSE]
    cell_types <- factor(cell_types_df$cell_type)
    names(cell_types) <- cell_types_df$barcodes
    
    nUMI <- colSums(counts) # Sum over rows (genes) to get per-cell nUMI

    return(Reference(counts, cell_types, nUMI))
}

# Load reference data for both labels
reference_healthy <- load_reference_data("healthy", tmp_dir)
reference_injured <- load_reference_data("injured", tmp_dir)

## Examine reference object (optional)
print(dim(reference_healthy@counts)) #observe Digital Gene Expression matrix
table(reference_healthy@cell_types) #number of occurences for each cell type

# Function to load spatial data
load_spatial_data <- function(sample_name, tmp_dir) {
  counts_file <- file.path(tmp_dir, paste0(sample_name, "_counts.csv"))
  coords_file <- file.path(tmp_dir, paste0(sample_name, "_coordinates.csv"))
  
  counts <- read.csv(counts_file, row.names = 1, check.names = FALSE)
  
  # Load coordinates (barcodes are in the first column)
  coords <- read.csv(coords_file, row.names = 1, check.names = FALSE)
  
  nUMI <- colSums(counts)
  
  return(SpatialRNA(coords, counts, nUMI))
}

# Load spatial data for each sample
sample_names <- config$samples  # List of sample names from config
pucks <- lapply(sample_names, function(sample) load_spatial_data(sample, tmp_dir))

# Run RCTD for each sample with the appropriate reference
results <- mapply(function(puck, sample_name) {
  reference <- if (sample_name %in% c("Spatial_1", "Spatial_2")) reference_healthy else reference_injured
  myRCTD <- create.RCTD(puck, reference, max_cores=4)
  myRCTD <- run.RCTD(myRCTD, doublet_mode='full')
  return(myRCTD)
}, pucks, sample_names, SIMPLIFY = FALSE)

# Save results
for (i in seq_along(sample_names)) {
  sample_output_dir <- file.path(input_dir, "indiv_samples", sample_names[i], "outs", "matrices")

  # Ensure the directory exists
  dir.create(sample_output_dir, recursive = TRUE, showWarnings = FALSE)
 
  saveRDS(results[[i]], file.path(sample_output_dir, paste0(sample_names[i], "_RCTD_results.rds")))
 
  # Convert sparse matrix to a regular data frame before saving as CSV
  # weights_df <- as.data.frame(as.matrix(results[[i]]@results$weights))
  # write.csv(weights_df, file.path(sample_output_dir, paste0(sample_names[i], "_RCTD_weights.csv")))
}

print("RCTD analysis completed successfully.")
