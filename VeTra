#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  6 11:25:56 2020

@author: Kim Junil, Guangzheng Weng
"""


import numpy as np
import sys
import os
import VelocityTrajectory as v
from matplotlib import pyplot as plt
import matplotlib
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
import rpy2.robjects.numpy2ri
from scipy import spatial
from sklearn.metrics.pairwise import cosine_similarity
matplotlib.use('Agg')
rpy2.robjects.numpy2ri.activate()
princurve = importr('princurve', on_conflict="warn")

def scale_(princ1):
    maxva = princ1.max()
    minva = princ1.min()
    ranges = maxva - minva
    
    new_leg = (princ1-minva)/ranges
    
    return new_leg
    
    
def find_sink_aera(wcc):
    cellIndex = wcc
    new = np.vstack((embedding[cellIndex,0], embedding[cellIndex,1]))
    new = np.vstack((new, merge_embed_delta[cellIndex,0]))
    new = np.vstack((new, merge_embed_delta[cellIndex,1]))
    new = new.T
    
    PrincipalCurve = new[:,:]
    B=PrincipalCurve
    nr,nc = B.shape
    Br = ro.r.matrix(B, nrow=nr, ncol=nc)
    ro.r.assign("B", Br)
    PrincipalCurveLambda = np.array(princurve.principal_curve(B).rx2('lambda'))
    cutsize = int(len(PrincipalCurveLambda)*0.1)
    index = np.argsort(PrincipalCurveLambda)
    inver_index = np.argsort(-PrincipalCurveLambda)
    
    mid = int(len(PrincipalCurveLambda)/2)
    mid_index = index[int(mid-cutsize/2):int(mid+cutsize/2)]
    mid_cell_index = cellIndex[mid_index]
    mid_ave_vec = np.average(delta_embedding[mid_cell_index,:], axis=0)
    
    aera1_index = index[:cutsize]
    aera1_cell_index = cellIndex[aera1_index]
    aera2_index = index[-cutsize:]
    aera2_cell_index = cellIndex[aera2_index]
    
    aera1_ave_emd = np.average(embedding[aera1_cell_index,:], axis=0)
    aera2_ave_emd = np.average(embedding[aera2_cell_index,:], axis=0)
    
    
    aera1_cell_bound = cellIndex[aera1_index[0]]
    aera2_cell_bound = cellIndex[aera2_index[-1]]
    
    
    cosine_j = 1 - spatial.distance.cosine(embedding[aera2_cell_bound,:]-embedding[aera1_cell_bound,:], mid_ave_vec)
    
    if cosine_j > 0.2:
        #print("good")
        #head_cell is an index for principalambda
        head_cell = aera1_index[0]
        tail_cell = aera2_index[-1]
        
        head_cell_bound = aera1_cell_bound
        tail_cell_bound = aera2_cell_bound
        
        head_aera = cellIndex[aera1_index[:5]]
        tail_aera = cellIndex[aera2_index[-5:]]
       
    else:
        #print("bad")
        tail_cell = aera1_index[0]
        head_cell = aera2_index[-1]
        
        
        tail_cell_bound = aera1_cell_bound
        head_cell_bound = aera2_cell_bound
        
        tail_aera = cellIndex[aera1_index[:5]]
        head_aera = cellIndex[aera2_index[-5:]]
    
    if PrincipalCurveLambda[head_cell] - PrincipalCurveLambda[tail_cell] > 0:
        PrincipalCurveLambda[index] =  PrincipalCurveLambda[inver_index]
    
    return PrincipalCurveLambda




file_name = ""
embedding,delta_embedding, deltaThreshold,WCCsizeCutoff,clusternumber,outfolder = sys.argv[1:]


embedding = np.loadtxt(embedding)
delta_embedding = np.loadtxt(delta_embedding)
merge_embed_delta = np.concatenate((embedding , delta_embedding), axis = 1)
deltaThreshold = int(float(deltaThreshold))
clusternumber = int(float(clusternumber))
WCCsizeCutoff = int(float(WCCsizeCutoff))
re = v.main(embedding,delta_embedding, deltaThreshold,WCCsizeCutoff,clusternumber)
embeddingGraphShortestpath=re[0]
embeddingGraphSPclusters=re[1]
graphSource=re[2]
graphTarget=re[3]
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

#np.savetxt(outfolder+'/embeddingGraphShortestpath.csv', embeddingGraphShortestpath)
#np.savetxt(outfolder+'/embeddingGraphSPclusters.csv', embeddingGraphSPclusters, delimiter = '\t', fmt='%f')
#np.savetxt(outfolder+'/graphSource.csv', graphSource, delimiter = '\t', fmt='%f')
#np.savetxt(outfolder+'/graphTarget.csv', graphTarget, delimiter = '\t', fmt='%f')
#


# =============================================================================
#plot results
# =============================================================================

#plot whole pathway in one figure

print(len(embeddingGraphShortestpath))
print(embeddingGraphSPclusters)
co = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf','black', 'peru','lime']
plt.figure(None,(20,20))
plt.quiver(embedding[:,0],embedding[:,1],delta_embedding[:,0],delta_embedding[:,1],color='#61c0bf')
for k in range(len(embeddingGraphShortestpath)):
    tmp = embeddingGraphShortestpath[k]
    plt.scatter(embedding[tmp,0],embedding[tmp,1],  marker='o', color =co[int(embeddingGraphSPclusters[k])])

plt.savefig(outfolder+'/'+file_name+'clusters'+'.pdf')
#plot every single pathway
for k in range(len(embeddingGraphShortestpath)):
    tmp = embeddingGraphShortestpath[k]
    plt.figure(None,(20,20))
    plt.quiver(embedding[:,0],embedding[:,1],delta_embedding[:,0],delta_embedding[:,1],color='#61c0bf')
    plt.scatter(embedding[tmp,0],embedding[tmp,1],  marker='o', color =co[int(embeddingGraphSPclusters[k])])
    plt.savefig(outfolder+'/'+file_name+'cluster_'+str(k+1)+'.pdf')

#run principal curve for each pathway
pseudo_time = []
for k in range(len(embeddingGraphShortestpath)):
    pathway = embeddingGraphShortestpath[k]
    PrincipalCurveLambda = find_sink_aera(pathway)
    scaled_Lambda = scale_(PrincipalCurveLambda)
    np.savetxt(outfolder+'/'+'pathway_'+str(k+1)+'.lambda.txt', scaled_Lambda, fmt='%f')
    np.savetxt(outfolder+'/'+'pathway_'+str(k+1)+'_selected_cell.txt', pathway, fmt='%d')
    pseudo_time.append(scaled_Lambda)
    
    plt.figure(None,(20,20))
    plt.quiver(embedding[:,0],embedding[:,1],delta_embedding[:,0],delta_embedding[:,1],color='#61c0bf')
    plt.scatter(embedding[pathway,0],embedding[pathway,1],  marker='o', c = scaled_Lambda)
    cb = plt.colorbar()
    plt.xticks([])
    plt.yticks([])
    plt.xlabel("UMAP1", fontsize = 45)
    plt.ylabel("UMAP2", fontsize = 45)
    cb.ax.tick_params(labelsize=30)
    plt.savefig(outfolder+'/'+file_name+'pathway'+str(k+1)+'.pdf', dpi=600)
    
#plot pseudo time for whole paths
plt.figure(None,(20,20))
plt.quiver(embedding[:,0],embedding[:,1],delta_embedding[:,0],delta_embedding[:,1],color='#61c0bf')
for k in range(len(pseudo_time)):
    pathway = embeddingGraphShortestpath[k] 
    plt.scatter(embedding[pathway,0],embedding[pathway,1],  marker='o', c = pseudo_time[k])

cb = plt.colorbar()
plt.xticks([])
plt.yticks([])
plt.xlabel("UMAP1", fontsize = 45)
plt.ylabel("UMAP2", fontsize = 45)
cb.ax.tick_params(labelsize=30)
plt.savefig(outfolder+'/'+file_name+'whole_paths'+'.pdf', dpi=600)
