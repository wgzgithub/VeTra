# VelocityTrajectory
A tool for infering scRNA_seq trajectory from RNA velocity
# Dependency
Python package:   rpy2
<br>R package:   princurve  

The script wokrs on python but it depends on a R package "princurve".   
Users should install R 3.6(or higher version) and install "princurve" for your system.

# Usage
```
python run.py [embedding] [delta_embedding] [deltaThreshold] [WCCsizeCutoff] [Clusternumber] [outfolder]
```
Command for example provided:
```
python run.py embedding.txt delta_embedding.txt 15 60 3 ./results
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

2.delta_embedding - a txt file with N cells in rows and 2 velocity coodinates in columns. Users can run "velocyto"(http://velocyto.org/) to get delta_embedding.  
#### example  
```
-1.662885950152473424e-02 -1.019793607352199594e-01
2.932140719083742297e-02 2.628113801391315230e-01
.
.
-2.781920407232985060e-02 -1.070944132601359122e-02
```
3.deltaThreshold - The number of near cells for a each single cell.

4.WCCsizeCutoff - The minimum size of weakly connective components.  

5.Clusternumber - Maxmium clusters for hireachy cluster. It equals to how many pathways can be found.  

6.outfolder - A absolute or relative path for output. 


# Output

1. Some figiurs for pathways with pseudo time.

2. Some files record the pathways containing cell indexs.

3. Some files record the pseudo time for each single pathway.

