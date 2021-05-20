# VeTra
A tool for infering scRNA_seq trajectory from RNA velocity

# Citation
https://www.biorxiv.org/content/10.1101/2020.09.01.277095v1

# Dependency
  	openmpi
	JPype
	rpy2
	princurve
 	sklearn
  	scipy
  	numpy
  	statsmodels
  

The script works on python.
R 3.6(or higher version) and "princurve" package in R are required.
TF(Transfer entroy) inference needs openmpi installation on Linux OS. If you only want to infer trajectories of datasets, it is unnessary to install "openmpi" and "JPype".

## 1. Run VeTenet using expression data in a csv file and RNA velocity files in txt format
#### Initialize an example
	import VeTra as vt
	ex1 = vt.VeTra("embedding.txt", "delta_embedding.txt")
	
#### Execute trajectory inference 
	
	ex1.vetra(deltaThreshold=12, WCCsizeCutoff=5, clusternumber=3)


#### Run TENET for all trajectories 
	
	ex1.run_tenet_tf(expression="chroman_exp_filtered.csv", thread= 15, history_len = 1, species = 'mouse', bulk_run=True)
	
#### Make GRNs for  all trajectories 
	
	ex1.makeGRN_tf(method= "links", threshold = 1000, bulk_run=True)


###### (1) expression_file - a csv file with N cells in the rows and M genes in the columns (same format with wishbone pseudotime package).

		GENE_1	GENE_2	GENE_3	...	GENE_M

	CELL_1	

	CELL_2

	CELL_3

	.
	.
	.

	CELL_N


###### (2) embedding_file - a txt file with N cells in rows and 2D embedding coordinates in columns. The 2D embedding can be results of PCA, tSNE, or UMAP.

	-1.662885950152473424e-02 -1.019793607352199594e-01
	.
	.
	.
	2.932140719083742297e-02 2.628113801391315230e-01

###### (3) delta_embedding_file - a txt file with N cells in rows and 2D velocity vectorss in columns. Users can run "velocyto"(http://velocyto.org/) or "scVelo" (https://github.com/theislab/scvelo) to get delta_embedding.

	-6.344244118975589153e+00 1.268898329669120084e+00
	.
	.
	.
	-3.511296625112249714e+00 -3.950214438779450082e-02

###### (4) history_length - the length of history. In the benchmark data TENET provides best result when the length of history set to 1.

#### Output

	TE_result_matrix.txt - TEij, M genes x M genes matrix representing the causal relationship from GENEi to GENEj.

	TE	GENE_1	GENE_2	GENE_3	...	GENE_M
	GENE_1	0	0.05	0.02	...	0.004
	GENE_2	0.01	0	0.04	...	0.12
	GENE_3	0.003	0.003	0	...	0.001
	.
	.
	.
	GENE_M	0.34	0.012	0.032	...	0


