# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 09:51:05 2020

@author: gzw
"""

from tqdm import tqdm
from numpy import *
import numpy as np
from scipy.spatial.distance import cdist
from scipy.spatial.distance import squareform
from scipy.spatial.distance import pdist
import networkx as nx
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import shortest_path
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from collections import Counter
import scipy.cluster.hierarchy as sch
from sklearn.cluster import AgglomerativeClustering

def VelocityTrajectory(embedding,delta_embedding,deltaThreshold,SPnumber,graphClusterCutoff,WCCsizeCutoff):
       
    
    #Construct network based on velocity vectors
    merge_embed_delta = np.concatenate((embedding , delta_embedding), axis = 1)
    embeddingVectorDistances = squareform(pdist(merge_embed_delta,'seuclidean'))
    embedding_toward = embedding+delta_embedding
    embedding_distance = cdist(embedding_toward,embedding,metric='euclidean')
    
    graphSource=[]
    graphTarget=[]
    
    for i in range(embedding.shape[0]):
        sortIndex = np.argsort(embedding_distance[i,:])
        for j in range(deltaThreshold):
            if sortIndex[j]!= i and ((embedding[sortIndex[j],:]- embedding[i,:])*delta_embedding[i,:]>0).sum() ==2:
                graphSource.append(i)
                graphTarget.append(sortIndex[j])
                
    embeddingGraph = nx.DiGraph()
    edges = []
    for i in range(len(graphSource)):
        each_edge = [graphSource[i], graphTarget[i]]
        edges.append(each_edge)
    
    embeddingGraph.add_edges_from(edges)
    
    
        
    bins = [i for i in nx.weakly_connected_components(embeddingGraph)]
    numbersBins = [len(c) for c in nx.weakly_connected_components(embeddingGraph)]
    numbersBinsSorted = np.array(sorted(numbersBins,reverse = True))
    numbersBins= np.array(numbersBins)
    #sortIndex: get index which WCC > threshold
    sortIndex = np.argsort(-numbersBins)
    sortIndex=np.where(numbersBins>WCCsizeCutoff)
    sortIndex = sortIndex[0]
    sortIndex.tolist()
    #Initialize matrixs
    embeddingGraphShortestpathTotal=np.zeros((0,embedding.shape[0]))
    embeddingGraphSPclustersTotal=np.zeros(1)
    embeddingGraphSPstart=np.zeros(embedding.shape[0]).T
    embeddingGraphSPend=np.zeros(embedding.shape[0]).T
    
    
    WCCindex=sortIndex
    for i in tqdm(WCCindex):
        subgra_nodes = list(bins[i])
        subgra_nodes = sorted(subgra_nodes)
        embeddingGraphWCC=embeddingGraph.subgraph(subgra_nodes)
        nodeIndex = sorted(embeddingGraphWCC.nodes)
        
        #Identify long shortest paths
        embeddingGraphDistances=[i for i in nx.all_pairs_shortest_path_length(embeddingGraphWCC)]
        mat = nx.to_scipy_sparse_matrix(embeddingGraphWCC, nodelist=nodeIndex)
        embeddingGraphDistances= shortest_path(csgraph=mat)
        embeddingGraphDistances[np.isinf(embeddingGraphDistances)]=-1
        embeddingGraphShortestpath=np.zeros((0, embedding.shape[0]))
        new = []
        
        for distanceIndex in range(int(embeddingGraphDistances.max())):
            if len(embeddingGraphShortestpath)>SPnumber:
                break
            else:
                [filterIndex1,filterIndex2]=np.where(embeddingGraphDistances==embeddingGraphDistances.max()-distanceIndex)
                newindex =np.argsort(filterIndex2)
                filterIndex1=filterIndex1[newindex]
                filterIndex2=filterIndex2[newindex]
                
                embeddingGraphShortestpathTemp=np.zeros((len(filterIndex1), embedding.shape[0]))
            
                for i in range(len(filterIndex1)):
                    SPtemp=nx.shortest_path(embeddingGraphWCC,subgra_nodes[filterIndex1[i]],subgra_nodes[filterIndex2[i]])
                    embeddingGraphShortestpathTemp[i,SPtemp]=1
                    embeddingGraphSPstart[SPtemp[0]]=1
                    embeddingGraphSPend[SPtemp[-1]]=1
                
                
                embeddingGraphShortestpath=np.concatenate((embeddingGraphShortestpath, embeddingGraphShortestpathTemp), axis=0)
                new.append(embeddingGraphShortestpathTemp)
    
                if len(filterIndex1)>1:
                    #Combine very close shortest paths
                    embeddingGraphSPlinkages=linkage(embeddingGraphShortestpath, method='single', metric='Hamming')
                    embeddingGraphSPclusters=fcluster(embeddingGraphSPlinkages,t=graphClusterCutoff)
                    embeddingGraphShortestpath2=np.zeros((embeddingGraphSPclusters.max(),embedding.shape[0]))
                    
                    for clusterIndex in range(embeddingGraphSPclusters.max()):
                        #Notice: clusterIndex from 0 to max. but min cluster number is 1.
                        mat1 = np.sum(embeddingGraphShortestpath[embeddingGraphSPclusters==clusterIndex+1,:], axis=0)
                        cc = embeddingGraphShortestpath[embeddingGraphSPclusters==clusterIndex+1,:]
                        
                        
                        embeddingGraphShortestpath2[clusterIndex, np.where(mat1>0)] = 1
                    embeddingGraphShortestpath=embeddingGraphShortestpath2;
        
        #Cluster paths based on vectors
        embeddingGraphSPdistances=[]
        count = 0
        for i in range(len(embeddingGraphShortestpath)):
            for j in range(len(embeddingGraphShortestpath)):
                if i<j:
                    
                    temp_index1 = np.where(embeddingGraphShortestpath[i,:]!=0)[0]
                    temp_index2 = np.where(embeddingGraphShortestpath[j,:]!=0)[0]
                    temp = np.zeros((len(temp_index1), len(temp_index2)))
                    count+=1
                    for n in range(len(temp_index1)):
                        for m in range(len(temp_index2)):
                            temp[n,m] = embeddingVectorDistances[temp_index1[n], temp_index2[m]]
                    
                    c_min = np.min(temp,axis=0).tolist()
                    r_min = np.min(temp.T,axis=0).tolist()
                    tmp_distance = mean(c_min+r_min)
                    
                    embeddingGraphSPdistances.append(tmp_distance)
                    
        embeddingGraphSPlinkages=linkage(embeddingGraphSPdistances,'complete')
        
        embeddingGraphSPdistances=squareform(embeddingGraphSPdistances)
        
        embeddingGraphSPclusters=fcluster(embeddingGraphSPlinkages,t=1, criterion='maxclust')
        WithinDistance=np.mean(embeddingGraphSPdistances)
        
        for NumberOfClusters in range(2,11):
            embeddingGraphSPclusters=fcluster(embeddingGraphSPlinkages,t=NumberOfClusters,criterion='maxclust')
            WithinDistanceNext=[]
            BetweenDistanceNext=[]
            
            for i in range(max(embeddingGraphSPclusters)):
                
                for j in range(max(embeddingGraphSPclusters)):
                    
                    
                    index1 = np.where(embeddingGraphSPclusters==i+1)[0]
                    index2 = np.where(embeddingGraphSPclusters==j+1)[0]
                    
                    temp_mat = np.zeros((len(index1), len(index2)))
                    for n in range(len(index1)):
                        for m in range(len(index2)):
                            temp_mat[n,m] = embeddingGraphSPdistances[index1[n], index2[m]]
                    
                    
                    distanceTemp=np.mean(temp_mat)
                    
                    if i==j:
                        WithinDistanceNext.append(distanceTemp)
                    elif i<j:
                        BetweenDistanceNext.append(distanceTemp)
                 
            if NumberOfClusters==2:
                if WithinDistance<max(WithinDistanceNext):
                    NumberOfClusters=NumberOfClusters-1
                    embeddingGraphSPclusters=fcluster(embeddingGraphSPlinkages, t=NumberOfClusters,criterion='maxclust')
                    break
                else:
                    WithinDistance=WithinDistanceNext
                    BetweenDistance=BetweenDistanceNext
            
            else:
                judge1 = np.median(BetweenDistanceNext)/np.median(WithinDistanceNext)
                judge2 = np.median(BetweenDistance)/np.median(WithinDistance)
                if judge1 < judge2:
                    NumberOfClusters=NumberOfClusters-1
                    embeddingGraphSPclusters=fcluster(embeddingGraphSPlinkages,t=NumberOfClusters,criterion='maxclust')
                    break
                else:
                    WithinDistance=WithinDistanceNext
                    BetweenDistance=BetweenDistanceNext
        
        #print(embeddingGraphSPclusters.shape)
        if embeddingGraphSPclustersTotal.any():
            #print("ok")
            embeddingGraphSPclusters=embeddingGraphSPclusters+np.max(embeddingGraphSPclustersTotal)
        
        
        embeddingGraphShortestpathTotal=np.concatenate((embeddingGraphShortestpathTotal,embeddingGraphShortestpath),axis=0)
        embeddingGraphSPclustersTotal=np.concatenate((embeddingGraphSPclustersTotal,embeddingGraphSPclusters),axis=0)
    
        
    embeddingGraphShortestpath=embeddingGraphShortestpathTotal
    embeddingGraphSPclusters=embeddingGraphSPclustersTotal

    
    return embeddingGraphShortestpath,embeddingGraphSPclusters,embeddingGraphSPstart,embeddingGraphSPend,graphSource,graphTarget