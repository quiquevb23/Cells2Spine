The way to run the scripts for single-cell analysis of the samples is:

	- create_h5 if not created from count matrix from CellRanger (when we have barcodes, mtx and features)
	- qc_metrics to obtain qc_metrics
	- filter_scdblfinder to filter based on MAD and calculate scores of Doublets
	- normalize_clustering the name says it all
		- here we can run harmony integrate to integrate 2 or more samples or
		- continue with annotation independently

*The scripts filtering alone and scDblFinder.R can be run to separate the steps, but scDblFinder not working


