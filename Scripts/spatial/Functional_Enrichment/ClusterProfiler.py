'''
Script to do Functional Enrichment analysis (ORA) on list of DEGs
'''

import pandas as pd
from rpy2 import robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
import os
import argparse


def parse_args():
    parser = argparse.ArgumentParser(description="Process directories for single-cell data.")
    
    # Define the arguments for base_dir and output_base_dir
    parser.add_argument('--base_dir', type=str, required=True, 
                        help="Base directory for input data.")
    parser.add_argument('--output_base_dir', type=str, required=True, 
                        help="Base directory for output data.")
    
    # Parse the arguments
    return parser.parse_args()

args = parse_args()

base_dir = args.base_dir
output_base_dir = args.output_base_dir

pandas2ri.activate()

# Ensure that renv is activated in the R environment
ro.r('''
    library(renv)
    renv::restore()
''')

# Load R packages and functions
ro.r('''
    library(clusterProfiler)
    library(org.Rn.eg.db) #Genome wide annotation for rat
''')

# Set the output directory path
output_dir = os.path.join(output_base_dir, "DEA_areas")

# Loop through directories in the output directory
for dir_name in os.listdir(output_dir):
    dir_path = os.path.join(output_dir, dir_name)
    
    # Check if it is a directory and does not start with "Spatial"
    if os.path.isdir(dir_path) and not dir_name.startswith("Spatial"):
        for file_name in os.listdir(dir_path):
            # Check if the file is a CSV
            if file_name.endswith("pseudobulk.csv"):
                file_path = os.path.join(dir_path, file_name)
                deg_df = pd.read_csv(file_path)

                output_enr_dir = os.path.join(dir_path, "ClusterProfilerEnrichment")
                if not os.path.exists(output_enr_dir):
                    os.makedirs(output_enr_dir)

                # Filter significant DEGs (optional)
                # For example, keeping only DEGs with FDR < 0.05
                filtered_degs = deg_df[deg_df['FDR'] < 0.05]  # Adjust threshold as needed

                # Separate upregulated and downregulated genes based on logFC
                upregulated_genes = filtered_degs[filtered_degs['logFC'] > 0]['genes'].tolist()  # Replace 'logFC' with your log fold change column name
                downregulated_genes = filtered_degs[filtered_degs['logFC'] < 0]['genes'].tolist()  # Same here

                # Convert the gene lists to R vectors
                upregulated_genes_r = ro.StrVector(upregulated_genes)
                downregulated_genes_r = ro.StrVector(downregulated_genes)

                with localconverter(ro.default_converter + pandas2ri.converter):
                    # Perform GO enrichment analysis for all ontologies (BP, CC, MF)
                    go_enrich_up_bp = ro.r['enrichGO'](gene = upregulated_genes_r, OrgDb = 'org.Rn.eg.db', keyType = 'SYMBOL', ont = 'BP', pvalueCutoff = 0.05)
                    go_enrich_up_cc = ro.r['enrichGO'](gene = upregulated_genes_r, OrgDb = 'org.Rn.eg.db', keyType = 'SYMBOL', ont = 'CC', pvalueCutoff = 0.05)
                    go_enrich_up_mf = ro.r['enrichGO'](gene = upregulated_genes_r, OrgDb = 'org.Rn.eg.db', keyType = 'SYMBOL', ont = 'MF', pvalueCutoff = 0.05)

                    # Perform GO enrichment analysis for downregulated genes (same ontologies)
                    go_enrich_down_bp = ro.r['enrichGO'](gene = downregulated_genes_r, OrgDb = 'org.Rn.eg.db', keyType = 'SYMBOL', ont = 'BP', pvalueCutoff = 0.05)
                    go_enrich_down_cc = ro.r['enrichGO'](gene = downregulated_genes_r, OrgDb = 'org.Rn.eg.db', keyType = 'SYMBOL', ont = 'CC', pvalueCutoff = 0.05)
                    go_enrich_down_mf = ro.r['enrichGO'](gene = downregulated_genes_r, OrgDb = 'org.Rn.eg.db', keyType = 'SYMBOL', ont = 'MF', pvalueCutoff = 0.05)

                # Extract the result slot
                go_enrich_up_bp_result = ro.r['slot'](go_enrich_up_bp, "result")
                go_enrich_up_cc_result = ro.r['slot'](go_enrich_up_cc, "result")
                go_enrich_up_mf_result = ro.r['slot'](go_enrich_up_mf, "result")
                go_enrich_down_bp_result = ro.r['slot'](go_enrich_down_bp, "result")
                go_enrich_down_cc_result = ro.r['slot'](go_enrich_down_cc, "result")
                go_enrich_down_mf_result = ro.r['slot'](go_enrich_down_mf, "result")

                # Convert the R data frames back to Pandas DataFrames after the R function call
                with localconverter(ro.default_converter + pandas2ri.converter):
                    go_enrich_up_bp_df = pandas2ri.rpy2py(go_enrich_up_bp_result)
                    go_enrich_up_cc_df = pandas2ri.rpy2py(go_enrich_up_cc_result)
                    go_enrich_up_mf_df = pandas2ri.rpy2py(go_enrich_up_mf_result)

                    go_enrich_down_bp_df = pandas2ri.rpy2py(go_enrich_down_bp_result)
                    go_enrich_down_cc_df = pandas2ri.rpy2py(go_enrich_down_cc_result)
                    go_enrich_down_mf_df = pandas2ri.rpy2py(go_enrich_down_mf_result)

                # Print or save the results
                print("GO Enrichment for Upregulated Genes (BP):")
                print(go_enrich_up_bp_df.head)

                print("GO Enrichment for Upregulated Genes (CC):")
                print(go_enrich_up_cc_df.head)

                print("GO Enrichment for Upregulated Genes (MF):")
                print(go_enrich_up_mf_df.head)

                # Optionally, print downregulated gene results
                print("GO Enrichment for Downregulated Genes (BP):")
                print(go_enrich_down_bp_df.head)

                print("GO Enrichment for Downregulated Genes (CC):")
                print(go_enrich_down_cc_df.head)

                print("GO Enrichment for Downregulated Genes (MF):")
                print(go_enrich_down_mf_df.head)

                # Save the output DataFrames to CSV files in the output directory
                go_enrich_up_bp_df.to_csv(f"{output_enr_dir}/go_enrichment_upregulated_BP.csv", index=False)
                go_enrich_up_cc_df.to_csv(f"{output_enr_dir}/go_enrichment_upregulated_CC.csv", index=False)
                go_enrich_up_mf_df.to_csv(f"{output_enr_dir}/go_enrichment_upregulated_MF.csv", index=False)

                go_enrich_down_bp_df.to_csv(f"{output_enr_dir}/go_enrichment_downregulated_BP.csv", index=False)
                go_enrich_down_cc_df.to_csv(f"{output_enr_dir}/go_enrichment_downregulated_CC.csv", index=False)
                go_enrich_down_mf_df.to_csv(f"{output_enr_dir}/go_enrichment_downregulated_MF.csv", index=False)
