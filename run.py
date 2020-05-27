#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  6 11:25:56 2020

@author: Kim Junil, Guangzheng Weng
"""

num_job=$2

import numpy as np
import sys
import os
import VelocityTrajectory as v
from matplotlib import pyplot as plt
import matplotlib
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
import rpy2.robjects.numpy2ri
matplotlib.use('Agg')






file_name = ""
embedding,delta_embedding, deltaThreshold,SPnumber,graphClusterCutoff,WCCsizeCutoff,thres_cellshape,outfolder = sys.argv[1:]


embedding = np.loadtxt(embedding)
delta_embedding = np.loadtxt(delta_embedding)
deltaThreshold = int(float(deltaThreshold))
SPnumber = int(float(SPnumber))
graphClusterCutoff = float(graphClusterCutoff)
WCCsizeCutoff = int(float(WCCsizeCutoff))
thres_cellshape = int(float(thres_cellshape))
re = v.main(embedding, delta_embedding, deltaThreshold, SPnumber, graphClusterCutoff, WCCsizeCutoff)
embeddingGraphShortestpath=re[0]
embeddingGraphSPclusters=re[1]
embeddingGraphSPclusters = np.delete(embeddingGraphSPclusters,0,axis=0)
embeddingGraphSPstart=re[2]
embeddingGraphSPend=re[3]
graphSource=re[4]
graphTarget=re[5]
embeddingSource=embedding[graphSource,:]
embeddingTarget=embedding[graphTarget,:]

# =============================================================================
#create a folder to save all return values
# =============================================================================

isExists=os.path.exists(outfolder+'')
if not isExists:
    os.makedirs(outfolder) 
else:
    pass

np.savetxt(outfolder+'/embeddingGraphShortestpath.csv', embeddingGraphShortestpath, delimiter = '\t')
np.savetxt(outfolder+'/embeddingGraphSPclusters.csv', embeddingGraphSPclusters, delimiter = '\t')
np.savetxt(outfolder+'/embeddingGraphSPstart.csv', embeddingGraphSPstart, delimiter = '\t')
np.savetxt(outfolder+'/embeddingGraphSPend.csv', embeddingGraphSPend, delimiter = '\t')
np.savetxt(outfolder+'/graphSource.csv', graphSource, delimiter = '\t')
np.savetxt(outfolder+'/graphTarget.csv', graphTarget, delimiter = '\t')


   



# =============================================================================
#plot results
# =============================================================================

uniqueClusters=np.unique(embeddingGraphSPclusters)
cellSizeClusters=[]
fig = plt.gcf()
fig.set_size_inches(12, 12)
for clusterIndex in range(len(uniqueClusters)):
    index = np.where(embeddingGraphSPclusters==clusterIndex+1)
    tmp = np.sum(embeddingGraphShortestpath[index[0],:], axis=0)
    cellSizeClusters.append(sum(tmp!=0))

cellSizeClusters = np.array(cellSizeClusters)
sortindex = np.argsort(-cellSizeClusters)+1
cellSizeSorted = -np.sort(-cellSizeClusters)

plt.quiver(embeddingSource[:,0],embeddingSource[:,1],embeddingTarget[:,0]-embeddingSource[:,0],
           embeddingTarget[:,1]-embeddingSource[:,1],0,
        scale=8, scale_units='inches', color='g')


for clusterIndex in sortindex:
    index = np.where(embeddingGraphSPclusters==clusterIndex)
    tmp = np.sum(embeddingGraphShortestpath[index[0],:], axis=0)
    cellIndex=np.where(tmp!=0)[0]
    if cellIndex.shape[0]>thres_cellshape:
        plt.scatter(embedding[cellIndex,0],embedding[cellIndex,1], marker=".")
        

plt.savefig(outfolder+'/'+file_name+'clusters'+'.pdf')
for clusterIndex in sortindex:
    index = np.where(embeddingGraphSPclusters==clusterIndex)
    tmp = np.sum(embeddingGraphShortestpath[index[0],:], axis=0)
    cellIndex=np.where(tmp!=0)[0]
    if cellIndex.shape[0]>thres_cellshape:
        plt.figure(clusterIndex+100)
        plt.quiver(embeddingSource[:,0],embeddingSource[:,1],embeddingTarget[:,0]-embeddingSource[:,0],embeddingTarget[:,1]-embeddingSource[:,1],0)
        plt.scatter(embedding[cellIndex,0],embedding[cellIndex,1], marker=".");
        plt.savefig(outfolder+'/'+file_name+'cluster_'+str(clusterIndex)+'.pdf');

for clusterIndex in sortindex:
    index = np.where(embeddingGraphSPclusters==clusterIndex)
    tmp = np.sum(embeddingGraphShortestpath[index[0],:], axis=0)
    cellIndex=np.where(tmp!=0)[0]
    i=cellIndex

    new = np.vstack((embedding[cellIndex,0], embedding[cellIndex,1]))
    new = np.vstack((new, delta_embedding[cellIndex,0]))
    new = np.vstack((new, delta_embedding[cellIndex,1]))
    new = new.T
    print(new.shape)  
    #np.savetxt('/home/gzw/work/test_data/'+str(clusterIndex)+'.txt', ( embedding[i,0].T,embedding[i,1].T,delta_embedding[i,0].T,delta_embedding[i,1].T)  )
    np.savetxt(outfolder+'/'+'cluster_'+str(clusterIndex)+'.txt', new)

# =============================================================================
#run principla curve
#return and save lamdba.txt
# =============================================================================


rpy2.robjects.numpy2ri.activate()
princurve = importr('princurve', on_conflict="warn")

for clusterIndex in sortindex:
    
    PrincipalCurve = np.loadtxt(outfolder+"//"+'cluster_'+str(clusterIndex)+".txt", delimiter=' ')
    PrincipalCurve = PrincipalCurve[:,1:3]
    B=PrincipalCurve
    nr,nc = B.shape
    Br = ro.r.matrix(B, nrow=nr, ncol=nc)
    ro.r.assign("B", Br)
    PrincipalCurveLambda = np.array(princurve.principal_curve(B).rx2('lambda'))
    np.savetxt(outfolder+'/'+'cluster_'+str(clusterIndex)+'.lambda.txt', PrincipalCurveLambda)


fig = plt.gcf()
fig.set_size_inches(12, 12)

plt.quiver(embeddingSource[:,0],embeddingSource[:,1],embeddingTarget[:,0]-embeddingSource[:,0],
           embeddingTarget[:,1]-embeddingSource[:,1],
        scale=7, scale_units='inches', color='chartreuse')

for clusterIndex in sortindex:

    PrincipalCurve = np.loadtxt(outfolder+"//"+'cluster_'+str(clusterIndex)+".txt", delimiter=' ')
    PrincipalCurveLambda = np.loadtxt(outfolder+"//"+'cluster_'+str(clusterIndex)+".lambda.txt")
    index = np.where(embeddingGraphSPclusters==clusterIndex)
    tmp = np.sum(embeddingGraphShortestpath[index[0],:], axis=0)
    cellIndex=np.where(tmp!=0)[0]
    startIndex = np.zeros(len(cellIndex))
    endIndex = np.zeros(len(cellIndex))

    for i in range(len(cellIndex)):
        tmp_start = np.where(embeddingGraphSPstart!=0)
        tmp_end = np.where(embeddingGraphSPend!=0)
        if cellIndex[i] in tmp_start[0]:
            startIndex[i] = 1
        elif cellIndex[i] in tmp_end[0]:
            endIndex[i] = 1
    
    if np.mean(PrincipalCurveLambda[np.where(PrincipalCurveLambda!=0)[0]])>np.mean(PrincipalCurveLambda[np.where(endIndex!=0)[0]]):
        PrincipalCurveLambda=max(PrincipalCurveLambda)-PrincipalCurveLambda
        plt.scatter(embedding[cellIndex,0],embedding[cellIndex,1], c=PrincipalCurveLambda, marker='o')
    
plt.savefig(outfolder+'/'+file_name+'final.pdf')

# =============================================================================
#run TENET
#return and save lamdba.txt
# =============================================================================

tmp_te = np.where(cellSizeSorted>thres_cellshape)
select_gene = []
pseudotime = []

for clusterIndex in sortindex[tmp_te]:
    PrincipalCurve = np.loadtxt(outfolder+"//"+'cluster_'+str(clusterIndex)+".txt", delimiter=' ')
    PrincipalCurveLambda = np.loadtxt(outfolder+"//"+'cluster_'+str(clusterIndex)+".lambda.txt")
    index = np.where(embeddingGraphSPclusters==clusterIndex)
    tmp = np.sum(embeddingGraphShortestpath[index[0],:], axis=0)
    cellIndex=np.where(tmp!=0)[0]
    startIndex = np.zeros(len(cellIndex))
    endIndex = np.zeros(len(cellIndex))

    for i in range(len(cellIndex)):
        tmp_start = np.where(embeddingGraphSPstart!=0)
        tmp_end = np.where(embeddingGraphSPend!=0)
        if cellIndex[i] in tmp_start[0]:
            startIndex[i] = 1
        elif cellIndex[i] in tmp_end[0]:
            endIndex[i] = 1
    
    if np.mean(PrincipalCurveLambda[np.where(PrincipalCurveLambda!=0)[0]])>np.mean(PrincipalCurveLambda[np.where(endIndex!=0)[0]]):
        PrincipalCurveLambda=max(PrincipalCurveLambda)-PrincipalCurveLambda
    
    save_tmp = []
    for i in range(embeddingGraphShortestpath.shape[1]):
        b = sum(cellIndex==i)
        save_tmp.append(b)
    
    save_tmp = np.array(save_tmp)
    
    np.savetxt(outfolder+'//'+file_name+'cluster_'+str(clusterIndex)+'_cell_select.txt',save_tmp,fmt='%d')

      
    save_tmp = []
    for i in range(embeddingGraphShortestpath.shape[1]):
        if sum(cellIndex==i):
            tmp = np.where(cellIndex==i)
            save_tmp.append(PrincipalCurveLambda[tmp])
            
        else:
            save_tmp.append(0)
    save_tmp = np.array(save_tmp)
       
    np.savetxt(outfolder+'/'+file_name+'cluster_'+str(clusterIndex)+'_pseudotime.txt', save_tmp,fmt='%f')



