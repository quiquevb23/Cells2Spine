# Cells2Spine
Scripts to analyze single-cell and spatial transcriptomics data from mice spinal cords upon injury

Steps for single-cell data analysis:

- Count matrix generation
- CellBender preprocessing of raw matrices to eliminate empty droplets
- Doublet identification with scDblFinder
- Integration of datasets from both experimental conditions
- QC metrics and optional additional filtering, normalization, dimensionality reduction and clustering
- Manual cluster annotation at different resolutions to identify cell types and subtypes
- Cell subtype differential abundance testing among conditions
- Perturbation analysis with Augur to find out most responsive cell types to injury
- Differential Gene Expression in pairs for the same subtype in both conditions to identify marker genes caused by or responsive to the injury
- GSEA and PA between pairs of the same subtype to identify pathways and biological functions up/down regulated due to the injury

Steps for spatial data analysis:

- Count matrix generation
- Preprocessing
- Registration to common coordinate reference
- Integration of all datasets (or only of those from same condition)
- QC metrics and optional additional filtering, normalization, dimensionality reduction and clustering:
	- feature selection: HVGs vs SVGs
	- clustering: non-spatial vs spatial (GraphST)
- Spot deconvolution from single-cell data
- Differential expression between clusters/spatial domains for one condition to identify markers for regions
- ... (continue)
