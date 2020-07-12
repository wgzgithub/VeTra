# -*- coding: utf-8 -*-
"""
Created on Fri Jun 26 10:33:27 2020

@author: Guangzheng Weng&Junil Kim
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
import getopt
import sys
from scipy import spatial
from sklearn.metrics.pairwise import cosine_similarity
import matplotlib
# from rpy2.robjects.packages import importr
# import rpy2.robjects as ro
# import rpy2.robjects.numpy2ri
# rpy2.robjects.numpy2ri.activate()
# princurve = importr('princurve', on_conflict="warn")
matplotlib.use('Agg')

def scale_(princ1):
    maxva = princ1.max()
    minva = princ1.min()
    ranges = maxva - minva
    
    new_leg = (princ1-minva)/ranges
    
    return new_leg
    


def main(embedding,delta_embedding, deltaThreshold,WCCsizeCutoff,clusternumber):

    delta_embedding = np.loadtxt("delta_embedding.txt")
    merge_embed_delta = np.concatenate((embedding , delta_embedding), axis = 1)
    embeddingVectorDistances = squareform(pdist(merge_embed_delta,'euclidean'))
    embedding_toward = embedding+delta_embedding
    embedding_distance = cdist(embedding_toward,embedding,metric='euclidean')
    embed_points_distance = squareform(pdist(embedding,'euclidean'))

    graphSource=[]
    graphTarget=[]


    for i in range(embedding.shape[0]):
        sortIndex = np.argsort(embedding_distance[i,:])
        two_points_index = np.argsort(embed_points_distance[i,:])
        cosine_list = []
        for j in range(deltaThreshold):
            if sortIndex[j] != i:
                pair_vector = delta_embedding[i,:]
                tmp = sortIndex[j]
                arrow_vector = delta_embedding[tmp,:]
                cosine = 1 - spatial.distance.cosine(pair_vector, arrow_vector)
                #print(cosine)
                if cosine < 0:
                    cosine = 0
                    cosine_list.append((cosine,j))
                else:
                    cosine_list.append((cosine,j))
        
        
        
        arrow_root_cos_list =[]
        for k in cosine_list:
            if k[0] > 0.5:
                pair_vector = embedding[sortIndex[k[1]],:]- embedding[i,:]
                arrow_vector = delta_embedding[i,:]
                cosine_n = 1 - spatial.distance.cosine(pair_vector, arrow_vector)
                if cosine_n < 0:
                    cosine_n = 0
                    arrow_root_cos_list.append((cosine_n, k[1]))
                else:
                    arrow_root_cos_list.append((cosine_n, k[1]))
        
                    
        if len(arrow_root_cos_list) > 0:
            arrow_root_cosines = [b[0] for b in arrow_root_cos_list]
            arrow_root_cosines = np.array(arrow_root_cosines)
            arrow_root_index = [b[1] for b in arrow_root_cos_list]
            arrow_root_index = np.array(arrow_root_index)
            cosine_index = np.argsort(-arrow_root_cosines)[0]
            
            original_j = sortIndex[arrow_root_index[cosine_index]]
            max_cosine = arrow_root_cosines[cosine_index]
            if max_cosine > 0.5:
                
                    
                graphSource.append(i)
                graphTarget.append(original_j)
           
    embeddingGraph = nx.DiGraph()
    edges = []

    for i in range(len(graphSource)):

        each_edge = (graphSource[i], graphTarget[i])
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
    embeddingGraphSPdistances=[]

    for i in sortIndex:
        space_path = np.zeros((1, embedding.shape[0]))
        subgra_nodes = list(bins[i])
        subgra_nodes = sorted(subgra_nodes)
        embeddingGraphWCC=embeddingGraph.subgraph(subgra_nodes)
        nodeIndex = sorted(embeddingGraphWCC.nodes)
        nodeIndex = np.array(nodeIndex)
        space_path[0, nodeIndex] = 1
        embeddingGraphShortestpathTotal=np.concatenate((embeddingGraphShortestpathTotal,space_path),axis=0)

    if embeddingGraphShortestpathTotal.shape[0] > 1:
        for i in range(len(embeddingGraphShortestpathTotal)):
            for j in range(len(embeddingGraphShortestpathTotal)):
                if i<j:
                    
                    p1_index =  np.where(embeddingGraphShortestpathTotal[i,:]>0)[0]
                    p2_index = np.where(embeddingGraphShortestpathTotal[j,:]>0)[0]
                    
                    p1 = merge_embed_delta[p1_index,:]
                    p2 = merge_embed_delta[p2_index,:]
                    
                    if len(p1) < len(p2):
                        
                        two_path_dis = cdist(p1,p2,metric='euclidean')
                        min_dis = np.min(two_path_dis, axis=1)
                        max_dis = min_dis.max()
                        embeddingGraphSPdistances.append(max_dis)
                    else:
                        two_path_dis = cdist(p2,p1,metric='euclidean')
                        min_dis = np.min(two_path_dis, axis=1)
                        max_dis = min_dis.max()
                        embeddingGraphSPdistances.append(max_dis)
                    #distance, path = fastdtw(p1, p2, dist=euclidean)
                    #embeddingGraphSPdistances.append(distance)
                    

                    
        embeddingGraphSPlinkages=linkage(embeddingGraphSPdistances,'complete')
        embeddingGraphSPclusters=fcluster(embeddingGraphSPlinkages,t=clusternumber,criterion='maxclust')
    else:
        embeddingGraphSPclusters = [1]
        
    deltaThreshold = deltaThreshold * 3
    new_cluster_list = []

    for i in range(embeddingGraphSPclusters.max()):
        same_cluster_index = np.where(embeddingGraphSPclusters == i+1)[0]
        same_cluster_path = embeddingGraphShortestpathTotal[same_cluster_index,:]
        tmp = np.sum(same_cluster_path, axis = 0)
        same_cluster_cells = np.where(tmp > 0)[0]
        
        append_cells = []
        for j in same_cluster_cells:
            
            
            two_p_index = np.argsort(embed_points_distance[j,:])
            
            for k in range(deltaThreshold):
                vec1 = delta_embedding[j,:]
                vec2 = delta_embedding[two_p_index[k]]
                cosine_ = 1 - spatial.distance.cosine(vec1, vec2)
                if cosine_ > 0.7 and two_p_index[k] not in same_cluster_cells:
                    append_cells.append(two_p_index[k])
                    
        new_cluster = np.concatenate((same_cluster_cells, append_cells))
        
        new_cluster = np.unique(new_cluster)
        new_cluster_list.append(new_cluster)
    
    new_embeddingGraphShortestpathTotal = np.array(new_cluster_list)
    return new_embeddingGraphShortestpathTotal, embeddingGraphSPclusters, graphSource, graphTarget
