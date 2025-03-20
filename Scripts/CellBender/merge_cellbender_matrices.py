import scanpy as sc
import anndata as ad
import os

base_path = "/storage/gge/Quique/Cells2SpineData/Pilot/single-cell/full_seq"

# Define base directories for each sample
sample_dirs = {
    "Single_Cell_1": ["HFNHJDSXC/outs/outs/cellbender_filtered_filtered.h5", "HWYVKDSXC/outs/outs/cellbender_filtered_filtered.h5"],
    "Single_Cell__3": ["H2K7NDSXF/outs/outs/cellbender_filtered_filtered.h5", "HFNHJDSXC/outs/outs/cellbender_filtered_filtered.h5"]
}


for sample, runs in sample_dirs.items():
    print(f"\nProcessing {sample}...")

    adata_list = []
    
    for run in runs:
        run_path = os.path.join(base_path, sample, run)
        if os.path.exists(run_path):
            print(f"  Loading: {run_path}")
            adata = sc.read_10x_h5(run_path)
            adata.obs["batch"] = run  # Keep track of which run each cell came from
            adata_list.append(adata)
        else:
            print(f"  Warning: {run_path} not found, skipping.")

    if adata_list:
        # Merge runs
        adata_combined = ad.concat(adata_list)
        output_path = "/storage/gge/Quique/Cells2SpineData/Pilot/single-cell/full_seq_cellbender_aggr/indiv_samples"
        sample_name = "Sample_1" if sample == "Single_Cell_1" else "Sample_3"
        # Save merged data as an AnnData file
        merged_output_path = os.path.join(output_path, sample_name, f"{sample_name}_cellbender_filtered.h5ad")
        adata_combined.write(merged_output_path)
        print(f"✅ Merged dataset saved at: {merged_output_path}")
    else:
        print(f"❌ No valid files found for {sample}, skipping merging.")

print("\n🎉 Merging completed!")


