from gprofiler import GProfiler
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Load the CSV file into a Pandas DataFrame
deg_csv_path = "/home/quiquevb/Cells2Spine/Cells2Spine/Outputs/spatial/DEA_dorsal/DEA_dorsal_pseudobulk.csv"
deg_df = pd.read_csv(deg_csv_path)

# Filter significant DEGs (optional)
# For example, keeping only DEGs with FDR < 0.05
filtered_degs = deg_df[deg_df['FDR'] < 0.05]  # Adjust threshold as needed

# Separate upregulated and downregulated genes based on logFC
upregulated_genes = filtered_degs[filtered_degs['logFC'] > 0]['genes'].tolist()  # Replace 'logFC' with your log fold change column name
downregulated_genes = filtered_degs[filtered_degs['logFC'] < 0]['genes'].tolist()  # Same here

# Initialize g:Profiler
gp = GProfiler(return_dataframe=True)

# Perform enrichment analysis for upregulated genes (for rat)
upreg_enrichment_results = gp.profile(
    organism='rnorvegicus',  # Rat organism code
    query=upregulated_genes,
    ordered=True
)

# Perform enrichment analysis for downregulated genes (for rat)
downreg_enrichment_results = gp.profile(
    organism='rnorvegicus',  # Rat organism code
    query=downregulated_genes,
    ordered=True
)

# Display top enrichment results for upregulated and downregulated genes
print("Enrichment Results for Upregulated Genes:")
print(upreg_enrichment_results.head())

print("Enrichment Results for Downregulated Genes:")
print(downreg_enrichment_results.head())

# Save the enrichment results to CSV files
upreg_enrichment_results.to_csv("/home/quiquevb/Cells2Spine/Cells2Spine/Outputs/spatial/DEA_dorsal/upregulated_enrichment_results.csv", index=False)
downreg_enrichment_results.to_csv("/home/quiquevb/Cells2Spine/Cells2Spine/Outputs/spatial/DEA_dorsal/downregulated_enrichment_results.csv", index=False)


# Load upregulated and downregulated enrichment results (assuming they were saved previously)
upreg_enrichment_path = "/home/quiquevb/Cells2Spine/Cells2Spine/Outputs/spatial/DEA_dorsal/upregulated_enrichment_results.csv"
downreg_enrichment_path = "/home/quiquevb/Cells2Spine/Cells2Spine/Outputs/spatial/DEA_dorsal/downregulated_enrichment_results.csv"

upreg_enrichment = pd.read_csv(upreg_enrichment_path)
downreg_enrichment = pd.read_csv(downreg_enrichment_path)

# Top 10 enriched terms (adjust this number if needed)
#top_n = 10

# Select top N enriched terms for both upregulated and downregulated genes based on the p-value
upreg_enrichment['-log10(p-value)'] = -np.log10(upreg_enrichment['p_value'])
downreg_enrichment['-log10(p-value)'] = -np.log10(downreg_enrichment['p_value'])

# Select top terms based on p-value
#upreg_top_terms = upreg_enrichment[['name', '-log10(p-value)']].sort_values(by='-log10(p-value)', ascending=False).head(top_n)
#downreg_top_terms = downreg_enrichment[['name', '-log10(p-value)']].sort_values(by='-log10(p-value)', ascending=False).head(top_n)

# Add a column to distinguish upregulated and downregulated genes
#upreg_top_terms['regulation'] = 'Upregulated'
#downreg_top_terms['regulation'] = 'Downregulated'

upreg_enrichment['regulation'] = 'Upregulated'
downreg_enrichment['regulation'] = 'Downregulated'

# Set up the plot
plt.figure(figsize=(12, 6))

# Combine both upregulated and downregulated terms into a single DataFrame
combined_enrichment = pd.concat([upreg_enrichment, downreg_enrichment])

# Add a column to distinguish upregulated and downregulated genes
combined_enrichment['direction'] = combined_enrichment['regulation'].apply(lambda x: 'Upregulated' if x == 'Upregulated' else 'Downregulated')

# Log-transform the p-values (negative log10)
combined_enrichment['neg_log10_p_value'] = -np.log10(combined_enrichment['p_value'])

# Plot for upregulated genes (blue)
upregulated_terms = combined_enrichment[combined_enrichment['direction'] == 'Upregulated']
plt.scatter(range(len(upregulated_terms)), upregulated_terms['neg_log10_p_value'], 
            color='blue', edgecolors='black', alpha=0.6, label='Upregulated')

# Plot for downregulated genes (red)
downregulated_terms = combined_enrichment[combined_enrichment['direction'] == 'Downregulated']
plt.scatter(range(len(downregulated_terms)), downregulated_terms['neg_log10_p_value'], 
            color='red', edgecolors='black', alpha=0.6, label='Downregulated')

# Add a threshold line for p-value < 0.05
plt.axhline(y=-np.log10(0.05), color='green', linestyle='--', label='p-value = 0.05')

# Customize the plot
plt.title('Manhattan Plot of Enrichment Results for Upregulated and Downregulated Genes')
plt.xlabel('Functional Terms')
plt.ylabel('-log10(p-value)')
plt.xticks(rotation=90)  # Rotate x-axis labels for readability
plt.legend()

# Save the plot
plt.tight_layout()
plt.savefig('/home/quiquevb/Cells2Spine/Cells2Spine/Outputs/spatial/DEA_dorsal/manhattan_plot_up_down_full.png')
plt.show()