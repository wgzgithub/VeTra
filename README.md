# VeTra
A tool for infering scRNA_seq trajectory from RNA velocity

# Citation
https://doi.org/10.1093/bioinformatics/btab364

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
#### Initialize an example(please make your script in the same folder of VeTra to import VeTra successfully)
	import VeTra as vt
	ex1 = vt.VeTra("embedding.txt", "delta_embedding.txt")
	
#### Execute trajectory inference 
	
	ex1.vetra(deltaThreshold=12, WCCsizeCutoff=5, clusternumber=3,cosine_thres=0.7, expand=2)
	#The default of cosine_thres is 0.7, which is a threshold to pick up similar neighbors of each cell. High cosine_thres usually finds less trajectories. 
	#The default of expand is 2. If you find the trajectories inferred dont cover all possbile cells, you will change expand larger.
	#deltaThreshold determines how many neighbors selected for each cell. It depends on how many cells in total. For example, deltaThreshold = 20~50 is good 	for the number of cells equalling to 3000

#### Run TENET for all trajectories 
	
	ex1.run_tenet_tf(expression="chroman_exp_filtered.csv", thread= 15, history_len = 1, species = 'mouse', bulk_run=True)
	
###### (1) history_length - the length of history. In the benchmark data TENET provides best result when the length of history set to 1.

	
#### Make GRNs for  all trajectories 
	
	ex1.makeGRN_tf(method= "links", threshold = 1000, bulk_run=True)


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


