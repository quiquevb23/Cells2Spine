- First run QC metrics to check

- Run SCTtransform for normalizing spots individually more accurately, better than sc.normalize

	- Then continue with clustering with GraphST.py to apply clusters on individual samples

	- Or run shared_domains.py to run GraphST on integrated embeddings (with Harmony) 

Alternatively run PASTE pairwise alignment to align samples

Then run GraphST with vertical alignment to find common domains
