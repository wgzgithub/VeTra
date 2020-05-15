# VelocityTrajectory
A tool for infering scRNA trajectory from RNA velocity
# Dependency
Python package: rpy2
R package: princurve
# Usage
python runVelocityTrajectory.py embedding.txt delta_embedding.txt 5 1000 0.01 100
# Input
1.embedding_file - a txt file with N rows and 2 columns. It can be obtained from "velocyto"(The attribute "embedding").
2.delta_embedding - a txt file with N rows and 2columns. It can be obtained from "velocyto"(The attribute "delta_embedding")
3.deltaThreshold - a threshold for constructing graph. The default is 5.
4.SPnumber - a threshold for constrcuting the long shortest paths. The default is 1000.
5.graphClusterCutoff - a threshold for clustering close shortest paths. The default is 0.01
6.WCCsizeCutoff - 
