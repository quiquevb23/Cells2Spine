The way to run the scripts for single-cell analysis of the samples is:

	- create_h5 if not created from count matrix from CellRanger (when we have barcodes, mtx and features)
	- qc_metrics to obtain qc_metrics
	- filter_scdblfinder to filter based on MAD and calculate scores of Doublets
	- normalize_HVG.py
        - clustering.py:
		performs DIM reduction (PCA) and saves them so that we can either:
			- run harmony integrate to integrate 2 or more samples (harmony_integrate.py)
                        - continue with clustering individual samples
	- annotation for indiv samples and annotation_joined for integrated samples

*The scripts filtering alone and scDblFinder.R can be run to separate the steps, but scDblFinder not working

Additionally we have XGBoost, Celltypist, integration_annotation to perform reference based annotation

