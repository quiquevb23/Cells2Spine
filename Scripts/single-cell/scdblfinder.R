# Load the renv package and activate the environment
if (!requireNamespace("renv", quietly = TRUE)) {
    install.packages("renv")
}
library(renv)
renv::restore()

print("Starting R script for scDblFinder analysis.")

# Load necessary libraries

library(SingleCellExperiment)
library(scDblFinder)
library(readr)

# Define file paths

base_dir <- "/storage/gge/Quique/Cells2SpineData/Pilot/single-cell/matrices/"

# Loop through all directories in base_dir (each directory is a sample)
sample_dirs <- list.dirs(base_dir, recursive = FALSE)
print(paste("Found", length(sample_dirs), "samples to process."))

# Iterate over each sample directory
for (sample_dir in sample_dirs) {
    print(paste("Processing sample directory:", sample_dir))
    
    # Find all CSV files in the sample directory
    csv_files <- list.files(sample_dir, pattern = "_data_mat.csv$", full.names = TRUE)
    
    # Iterate over each CSV file in the sample directory
    for (csv_file in csv_files) {
        print(paste("Processing CSV file:", csv_file))
        
        # Load the CSV file into a matrix
        data_mat <- as.matrix(read.csv(csv_file, row.names = 1))  # Adjust row.names if needed
        print(head(data_mat))
        
        print("Creating SingleCellExperiment object.")
        
        # Run scDblFinder on the data
        sce <- SingleCellExperiment(list(counts = data_mat)) 
 
        print("Running scDblFinder analysis.")

        results <- scDblFinder(sce)
        doublet_score <- results$scDblFinder.score
        doublet_class <- results$scDblFinder.class
        
        # Define output directory and save the result
        output_dir <- dirname(csv_file)
        save(doublet_score, doublet_class, file = file.path(output_dir, "scDblFinder_results.RData"))
        
        print(paste("Saved scDblFinder results for:", csv_file))
    }
}
