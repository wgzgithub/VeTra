# VeTra
A tool for infering scRNA_seq trajectory from RNA velocity
# Dependency
Python package:   rpy2
<br>R package:   princurve  

The script works on python.
R 3.6(or higher version) and "princurve" package in R are required.

# Usage
```
python VeTra [embedding] [delta_embedding] [deltaThreshold] [WCCsizeCutoff] [Clusternumber] [outfolder]
```
Command for example provided:
```
python VeTra embedding.txt delta_embedding.txt 15 60 3 ./results
```
# Input
1.embedding - a txt file with N cells in rows and 2D embedding coordinates in columns. The 2D embedding can be results of PCA, tSNE, or UMAP.
#### example  
```
-6.344244118975589153e+00 1.268898329669120084e+00
-3.511296625112249714e+00 -3.950214438779450082e-02
.
.
-2.453315327236579968e+00 -7.753647990937119205e-01
```

2.delta_embedding - a txt file with N cells in rows and 2D velocity vectorss in columns. Users can run "velocyto"(http://velocyto.org/) or "scVelo" (https://github.com/theislab/scvelo) to get delta_embedding.  
#### example  
```
-1.662885950152473424e-02 -1.019793607352199594e-01
2.932140719083742297e-02 2.628113801391315230e-01
.
.
-2.781920407232985060e-02 -1.070944132601359122e-02
```
3.deltaThreshold - The number of closest neighbor cells

4.WCCsizeCutoff - The minimum size of weakly connective components.  

5.Clusternumber - The number of clusters for hireachical clustering. It means how many trajectories can be found.  

6.outfolder - A absolute or relative path for output. 


# Output

1. 2D embedding figures colored by grouping of cells in the same stream of trajectory. Ex) cell_group1.pdf

2. 2D embedding figures colored by pseudo-time ordering. Ex) cell_group1with_pseudotime.pdf

3. Files for the lists of selected cells for each group of cells. Ex) cell_group1_selected_cell.txt

4. Files for the pseudo-time ordering information for each group of cells. Ex) cell_group1_pseudotime.txt

