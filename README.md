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
1.embedding_file - a txt file with N rows and 2 columns. It can be obtained from "velocyto"(The attribute "embedding").  
example  
'''python
-6.344244118975589153e+00 1.268898329669120084e+00
-3.511296625112249714e+00 -3.950214438779450082e-02
'''

2.delta_embedding - a txt file with N rows and 2columns. It can be obtained from "velocyto"(The attribute "delta_embedding").  

3.deltaThreshold - a threshold for constructing graph. The default is 5.  

4.SPnumber - a threshold for constrcuting the long shortest paths. The default is 1000.  

5.graphClusterCutoff - a threshold for clustering close shortest paths. The default is 0.01.  

6.WCCsizeCutoff - 
