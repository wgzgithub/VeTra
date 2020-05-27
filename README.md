# VelocityTrajectory
A tool for infering scRNA trajectory from RNA velocity
# Dependency
Python package:   rpy2
<br>R package:   princurve  

The script wokrs on python but it depends on a R package "princurve".   
You should install R_3.6(or higher version) and install "princurve" for your system.

# Usage
python runVelocityTrajectory.py embedding,delta_embedding,   deltaThreshold,SPnumber,graphClusterCutoff,WCCsizeCutoff,thres_cellshape,outfolder
# Input
1.embedding_file - a txt file with N cells in  rows and 2 embedding coordinates in columns.  
#### example  
```
-6.344244118975589153e+00 1.268898329669120084e+00
-3.511296625112249714e+00 -3.950214438779450082e-02
.
.
-2.453315327236579968e+00 -7.753647990937119205e-01
```

2.delta_embedding - a txt file with N cells in rows and 2 velocity coodinates in columns. You can run "velocyto"(http://velocyto.org/) to get delta_embedding .  
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

5.graphClusterCutoff - a threshold for clustering close shortest paths. The default is 0.01.  

6.WCCsizeCutoff - 
