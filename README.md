# VelocityTrajectory
A tool for infering scRNA_seq trajectory from RNA velocity
# Dependency
Python package:   rpy2
<br>R package:   princurve  

The script wokrs on python but it depends on a R package "princurve".   
Users should install R 3.6(or higher version) and install "princurve" for your system.

# Usage
```
python run.py [embedding] [delta_embedding] [deltaThreshold] [SPnumber] [graphClusterCutoff] [WCCsizeCutoff] [trajectoryClusterSize] [outfolder]
```
# Input
1.embedding - a txt file with N cells in  rows and 2 embedding coordinates in columns.  
#### example  
```
-6.344244118975589153e+00 1.268898329669120084e+00
-3.511296625112249714e+00 -3.950214438779450082e-02
.
.
-2.453315327236579968e+00 -7.753647990937119205e-01
```

2.delta_embedding - a txt file with N cells in rows and 2 velocity coodinates in columns. Users can run "velocyto"(http://velocyto.org/) to get delta_embedding .  
#### example  
```
-1.662885950152473424e-02 -1.019793607352199594e-01
2.932140719083742297e-02 2.628113801391315230e-01
.
.
-2.781920407232985060e-02 -1.070944132601359122e-02
```
3.deltaThreshold - a threshold for constructing graph. The default is 5.  

4.SPnumber - a threshold for constrcuting the long shortest paths. The default is 1000.  

5.graphClusterCutoff -  a threshold for the inconsistency coefficients (or inconsistent values) of nodes in the tree(agglomerative clusters from linkages)  

6.WCCsizeCutoff - 

7.trajectoryClusterSize -  

8.outfolder - A absolute or relative path for output. 


# Output

1.There are 6 csv format files including shortest path and clusters for embeeding graph.

2.cluster_x_cell_select.txt - a txt file of cell selection data for each cluster.("x" is order number of cluster)

3.cluster_x_pseudotime - a txt file of pseudotime data with N time points in the same order as the N cells for each cluster.("x" is order number of cluster)

