#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 24 17:32:08 2020

@author: gzw
"""

import os, glob,shutil
import numpy as np
import sys
from matplotlib import pyplot as plt
import matplotlib
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
import rpy2.robjects.numpy2ri
from scipy import spatial
from sklearn.metrics.pairwise import cosine_similarity
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
import numpy
import statsmodels.sandbox.stats.multicomp
import scipy.stats
import sys
import pandas as pd
import time



#matplotlib.use('Agg')
rpy2.robjects.numpy2ri.activate()
princurve = importr('princurve', on_conflict="warn")



comm = "./TENET expression_data.csv 10 trajectory.txt cell_select.txt 1"

#os.system(comm)


class VeTra:
    
    def __init__(self, embedding,delta_embedding):
        
        '''initialize necessary files'''
        self.embedding = np.loadtxt(embedding)
        self.delta_embedding = np.loadtxt(delta_embedding)
    
    def scale_(self, princ1):
        maxva = princ1.max()
        minva = princ1.min()
        ranges = maxva - minva
        
        new_leg = (princ1-minva)/ranges
        
        return new_leg
    
    
    def find_sink_aera(self, wcc):
        embedding = self.embedding
        merge_embed_delta = self.merge_embed_delta
        delta_embedding = self.delta_embedding
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



    def vetra(self,deltaThreshold,WCCsizeCutoff,clusternumber,cosine_thres=0.7, expand=2):

        #clean trajectories files created by last task
        os.system("rm cell_select*")
        os.system("rm trajectory*")


        
        embedding = self.embedding
        delta_embedding = self.delta_embedding
        merge_embed_delta = np.concatenate((embedding , delta_embedding), axis = 1)
        self.merge_embed_delta = merge_embed_delta
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
                if k[0] > cosine_thres:
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
                if max_cosine > cosine_thres:
                    
                        
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
            
        deltaThreshold = int(deltaThreshold * expand)
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
        
        embeddingGraphShortestpath = new_embeddingGraphShortestpathTotal
        self.embeddingGraphShortestpath = embeddingGraphShortestpath
        self.embeddingGraphSPclusters = embeddingGraphSPclusters
        self.graphSource = graphSource
        self.graphTarget = graphTarget
        
        outfolder = "./"
        file_name = ""
        
        isExists=os.path.exists(outfolder+'')
        if not isExists:
            os.makedirs(outfolder) 
        else:
            pass
        co = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf','black', 'peru','lime']
        
        pseudo_time = []
        
        
        self.number_of_trajectory = len(embeddingGraphShortestpath)
        for k in range(len(embeddingGraphShortestpath)):
            cell_select = np.zeros(len(embedding))
            trajectory  = np.zeros(len(embedding))
            pathway = embeddingGraphShortestpath[k]
            PrincipalCurveLambda = self.find_sink_aera(pathway)
            scaled_Lambda = self.scale_(PrincipalCurveLambda)
            
            cell_select[pathway] =1
            trajectory[pathway] = scaled_Lambda
            
            
            np.savetxt(outfolder+'/'+'trajectory_'+ str(k+1)+'.txt', trajectory, fmt='%f')
            np.savetxt(outfolder+'/'+'cell_select_'+str(k+1)+'.txt', cell_select, fmt='%d')
            
            
            pseudo_time.append(scaled_Lambda)
            
            plt.figure(None,(20,10))
            plt.quiver(embedding[:,0],embedding[:,1],delta_embedding[:,0],delta_embedding[:,1],color='#61c0bf')
            plt.scatter(embedding[pathway,0],embedding[pathway,1],  marker='o', c = scaled_Lambda,s=300)
            plt.axis("off")
            plt.title('trajectory_'+str(k+1),fontsize=25)
            plt.savefig(outfolder+'/'+file_name+'trajectory_'+str(k+1)+'.pdf', dpi=350)


        outfolder_ = "./TI_results"
        
        isExists=os.path.exists(outfolder_+'')
        if not isExists:
            os.makedirs(outfolder_)
            for i in range(len(embeddingGraphShortestpath)):
                shutil.move('trajectory_' + str(i+1) + '.pdf', './TI_results/')
                #os.system("mv trajectory_*pdf ./TI_results/")
                shutil.copy('trajectory_' + str(i+1) + '.txt', './TI_results/')
                shutil.copy('cell_select_' + str(i+1) + '.txt', './TI_results/') 
            
        else:
            shutil.rmtree(outfolder_, ignore_errors=True)
            os.makedirs(outfolder_)
            #os.system("mv trajectory_*pdf ./TI_results/")
            for i in range(len(embeddingGraphShortestpath)):
                shutil.move('trajectory_' + str(i+1) + '.pdf', './TI_results/')
                shutil.copy('trajectory_' + str(i+1) + '.txt', './TI_results/')
                shutil.copy('cell_select_' + str(i+1) + '.txt', './TI_results/') 
            
        
        print("Trajectory inference results saved in "+outfolder_)
        
    def divideIntoNstrand(self, listTemp, n):
        twoList = [ [] for i in range(n)]
        for i,e in enumerate(listTemp):
            twoList[i%n].append(e)
            
        len_each_list = [len(i) for i in twoList]
        return len_each_list
            



    def run_tenet_tf(self, expression, thread, history_len,species, bulk_run=False, trajectory=None, cell_select=None):    
        
        if bulk_run == True:
            #clean result files created by last task
            os.system("rm TE_result_matrix*")
            os.system("rm TE_result_all*")
            
            for i in range(self.number_of_trajectory):

                

                cores_run = thread                
                each_trajectory = 'trajectory_'+str(i+1)+'.txt'
                each_cell_select = 'cell_select_' + str(i+1)+'.txt'
                task_id = str(i+1)
                
                cmd = "./TENET_TF" + " "
                cmd = cmd + expression +" "
                cmd = cmd + str(cores_run) +" "
                cmd = cmd + each_trajectory + " "
                cmd = cmd + each_cell_select +" "
                cmd = cmd + str(history_len)+ " "
                cmd = cmd + species + " "
                print(cmd)
                os.system(cmd)
                os.rename("TE_result_matrix.txt",  "TE_result_matrix_"+task_id+".txt")
                os.rename("TE_result_all.csv",  "TE_result_all_"+task_id+".csv")
                time.sleep(5)



     


                
        else:
            self.tmp_trajectory = np.loadtxt(trajectory)
            self.tmp_cell_select = np.loadtxt(cell_select)
            
            cmd = "./TENET_TF" + " "
            cmd = cmd + expression +" "
            cmd = cmd + str(thread) +" "
            cmd = cmd + trajectory + " "
            cmd = cmd + cell_select +" "
            cmd = cmd + str(history_len)+ " "
            cmd = cmd + species
            print(cmd)
            os.system(cmd)
        

    def run_tenet(self, expression, thread, history_len,bulk_run=False, trajectory=None, cell_select=None):    
        
        if bulk_run == True:

            #remove files created by last task
            os.system("rm TE_result_matrix*")
            os.system("rm TE_result_all*")
            
            for i in range(self.number_of_trajectory):

                

                cores_run = thread                
                each_trajectory = 'trajectory_'+str(i+1)+'.txt'
                each_cell_select = 'cell_select_' + str(i+1)+'.txt'
                task_id = str(i+1)
                
                cmd = "./TENET" + " "
                cmd = cmd + expression +" "
                cmd = cmd + str(cores_run) +" "
                cmd = cmd + each_trajectory + " "
                cmd = cmd + each_cell_select +" "
                cmd = cmd + str(history_len)
                
                print(cmd)
                os.system(cmd)
                os.rename("TE_result_matrix.txt",  "TE_result_matrix_"+task_id+".txt")
                
                time.sleep(5)



     


                
        else:
            self.tmp_trajectory = np.loadtxt(trajectory)
            self.tmp_cell_select = np.loadtxt(cell_select)
            
            cmd = "./TENET_TF" + " "
            cmd = cmd + expression +" "
            cmd = cmd + str(thread) +" "
            cmd = cmd + trajectory + " "
            cmd = cmd + cell_select +" "
            cmd = cmd + str(history_len)
            
            print(cmd)
            os.system(cmd)


    def makeGRN_fdr(self,threshold, file_name,task_id):
        file_name=file_name

        ifile = open(file_name)
        line = ifile.readline()
        temp = line.split()
        gene_name=[]
        for i in range(len(temp)-1):
            gene_name.append(temp[i+1])

        cutOff=0
        sourceIndex=0
        TEnetwork=[]
        source=[]
        TE=[]
        target=[]
        for line in ifile:
            temp = line.split()
            for targetIndex in range(len(temp)-1):
                if float(temp[targetIndex+1])>cutOff:            
                    source.append(gene_name[sourceIndex])
                    TE.append(float(temp[targetIndex+1]))
                    target.append(gene_name[targetIndex])
            sourceIndex=sourceIndex+1
        ifile.close()

        TEzscore=(TE-numpy.mean(TE))/numpy.std(TE)
        TEpvalue=1-scipy.stats.norm.cdf(TEzscore)
        TEfdr=statsmodels.sandbox.stats.multicomp.multipletests(TEpvalue,alpha=0.05,method='fdr_bh')

        fdrCutoff=threshold
        ofile = open(file_name.replace(".txt",".fdr")+str(fdrCutoff)+".sif","w")
        for i in range(len(source)):
            if TEfdr[1][i]<fdrCutoff:
                ofile.write(source[i]+"\t"+str(TE[i])+"\t"+target[i]+"\n")
        ofile.close()

        file_name_save = file_name.replace(".txt",".fdr")+str(fdrCutoff)+".sif"
        #self.countOutdegree(file_name_save, task_id)

        return file_name_save, task_id

    def makeGRN_links(self, threshold, file_name, task_id):

        file_name=file_name

        ifile = open(file_name)
        line = ifile.readline()
        temp = line.split()
        gene_name=[]
        for i in range(len(temp)-1):
            gene_name.append(temp[i+1])

        cutOff=0
        sourceIndex=0
        TEnetwork=[]
        source=[]
        TE=[]
        target=[]
        for line in ifile:
            temp = line.split()
            for targetIndex in range(len(temp)-1):
                if float(temp[targetIndex+1])>cutOff:            
                    source.append(gene_name[sourceIndex])
                    TE.append(float(temp[targetIndex+1]))
                    target.append(gene_name[targetIndex])
            sourceIndex=sourceIndex+1
        ifile.close()

        TE=numpy.array(TE)

        TEsortIndex=numpy.argsort(TE)

        NumberOfLinks=threshold
        ofile = open(file_name.replace(".txt",".NumberOfLinks")+str(NumberOfLinks)+".sif","w")
        for i in range(int(NumberOfLinks)):
            ofile.write(source[TEsortIndex[-i-1]]+"\t"+str(TE[-i-1])+"\t"+target[TEsortIndex[-i-1]]+"\n")
        ofile.close()

        file_name_save = file_name.replace(".txt",".NumberOfLinks")+str(NumberOfLinks)+".sif"
        #self.countOutdegree(file_name_save, task_id)

        return file_name_save, task_id


    def makeGRN_TF_fdr(self,threshold,file_name, task_id):
        
        '''Make GRN by fdr'''
        
        ifile = open("gene_names")
        gene_names=[]
        for line in ifile:
            gene_names.append(line.replace("\n",""))
        ifile.close()

        file_name = file_name
        ifile = open(file_name)
        cutOff=0
        source=[]
        TE=[]
        target=[]
        for line in ifile:
            temp = line.replace("\n","").split(",")
            if float(temp[1])>cutOff:
                source.append(int(temp[0]))
                TE.append(float(temp[2]))
                target.append(int(temp[1]))
        ifile.close()

        TEzscore=(TE-numpy.mean(TE))/numpy.std(TE)
        TEpvalue=1-scipy.stats.norm.cdf(TEzscore)
        TEfdr=statsmodels.sandbox.stats.multicomp.multipletests(TEpvalue,alpha=0.05,method='fdr_bh')

        fdrCutoff=threshold
        ofile = open(file_name.replace(".csv",".fdr")+str(fdrCutoff)+".sif","w")
        for i in range(len(source)):
            if TEfdr[1][i]<fdrCutoff:
                ofile.write(gene_names[source[i]-1]+"\t"+str(TE[i])+"\t"+gene_names[target[i]-1]+"\n")
        ofile.close()
        file_name_save = file_name.replace(".csv",".fdr")+str(fdrCutoff)+".sif"
       
        #self.countOutdegree(file_name_save, task_id)

        return file_name_save, task_id
        
        
    def makeGRN_TF_links(self,threshold,file_name, task_id):
        '''make GRN by number of links'''
        
        ifile = open("gene_names")
        gene_names=[]
        for line in ifile:
            gene_names.append(line.replace("\n",""))
        ifile.close()

        
        ifile = open(file_name)
        cutOff=0
        source=[]
        TE=[]
        target=[]
        for line in ifile:
            temp = line.replace("\n","").split(",")
            if float(temp[1])>cutOff:
                source.append(int(temp[0]))
                TE.append(float(temp[2]))
                target.append(int(temp[1]))
        ifile.close()

        TE=numpy.array(TE)
        TEsortIndex=numpy.argsort(TE)

        NumberOfLinks=threshold
        ofile = open(file_name.replace(".csv",".NumberOfLinks")+str(NumberOfLinks)+".sif","w")
        for i in range(int(NumberOfLinks)):
            ofile.write(gene_names[source[TEsortIndex[-i-1]]-1]+"\t"+str(TE[-i-1])+"\t"+gene_names[target[TEsortIndex[-i-1]]-1]+"\n")
        ofile.close()
                
        
        file_name_save = file_name.replace(".csv",".NumberOfLinks")+str(NumberOfLinks)+".sif"
        #self.countOutdegree(file_name_save, task_id)


        return file_name_save, task_id
        
    
    def makeGRN_tf(self, method, threshold, TE_file=None, cell_select_file=None, trajectory_file=None, bulk_run = False):
        '''Users choose method to make GRN by giving threshold'''
        
        
        if bulk_run == True:
            print(self.number_of_trajectory)
            for i in range(self.number_of_trajectory):
                
                file_name = "TE_result_all_"+str(i+1)+".csv"
                if method == "fdr":
                    
                    print(file_name)
                    file_name_save, task_id = self.makeGRN_TF_fdr(threshold,file_name, str(i+1))
                    self.countOutdegree(file_name_save, task_id)
                
                elif method == "links":
                    
                    file_name_save, task_id= self.makeGRN_TF_links(threshold, file_name, str(i+1))
                    self.countOutdegree(file_name_save, task_id)
            #self.clean_save_results()
                
                
                
        else:
            file_name = TE_file
            if method == "fdr":
                
                file_name_, task_id = self.makeGRN_TF_fdr(threshold,file_name, "")
                self.countOutdegree_single_run(file_name_, cell_select_file, trajectory_file)
                
                
            
            elif method == "links":
                
                file_name_, task_id = self.makeGRN_TF_links(threshold,file_name, "")
                self.countOutdegree_single_run(file_name_, cell_select_file, trajectory_file)
                
                

        self.clean_save_results()


    def makeGRN(self, method, threshold, TE_file=None, cell_select_file=None, trajectory_file=None, bulk_run = False):
        '''Users choose method to make GRN by giving threshold'''
        
        
        if bulk_run == True:
            print(self.number_of_trajectory)
            for i in range(self.number_of_trajectory):
                file_name = "TE_result_matrix_" + str(i+1) + ".txt"
                if method == "fdr":
                    
                    print(file_name)
                    file_name_save, task_id = self.makeGRN_fdr(threshold,file_name, str(i+1))
                    self.countOutdegree(file_name_save, task_id)
                
                elif method == "links":
                    
                    file_name_save, task_id = self.makeGRN_links(threshold,file_name,str(i+1))
                    self.countOutdegree(file_name_save, task_id)

        
                
                
                
        else:
            file_name = TE_file
            if method == "fdr":

                file_name_, task_id = self.makeGRN_fdr(threshold,file_name, "")
                self.countOutdegree_single_run(file_name_, cell_select_file, trajectory_file)
                
            
            elif method == "links":
                
                file_name_, task_id  = self.makeGRN_links(threshold,file_name, "")
                self.countOutdegree_single_run(file_name_, cell_select_file, trajectory_file)
               
                
        self.clean_save_results()
            
    def countOutdegree(self, file_name, task_id):
        
        
        ifile = open(file_name)
        TFlist=[];TFlistOutdegree=[]
        for line in ifile:
            temp = line.split()
            if temp[0] not in TFlist:
                TFlist.append(temp[0])
                TFlistOutdegree.append(1)
            else:
                TFlistOutdegree[TFlist.index(temp[0])]=TFlistOutdegree[TFlist.index(temp[0])]+1
        TFlistOutdegreeIndex=numpy.argsort(TFlistOutdegree)
        
        out_folder = ""
        ofile = open(out_folder + file_name+".outdegree.txt","w")
        for i in range(len(TFlist)):
            ofile.write(TFlist[TFlistOutdegreeIndex[-i-1]]+"\t"+str(TFlistOutdegree[TFlistOutdegreeIndex[-i-1]])+"\n")
        ofile.close()
        
        
        filename = file_name+".outdegree.txt"
        traject_id = task_id
    
        tmp_cell_select = np.loadtxt("cell_select_" + traject_id + ".txt")
        tmp_trajectory = np.loadtxt("trajectory_" + traject_id + ".txt")

        TF = np.loadtxt(filename,dtype=object)
        pathway = np.where(tmp_cell_select==1)[0]
        trajectory = tmp_trajectory[pathway]

        fig = plt.figure(None,(20,10))
        left, bottom, width, height = 0.02, 0.15, 0.7, 0.7
        ax1 = fig.add_axes([left, bottom, width, height])
        ax1.quiver(self.embedding[:,0],self.embedding[:,1],self.delta_embedding[:,0],self.delta_embedding[:,1],color='black')
        ax1.scatter(self.embedding[pathway,0],self.embedding[pathway,1],  marker='o',c=trajectory, s=100)
        ax1.axis("off")


        TF_value = TF[:,1].astype(int)
        TF_name = TF[:,0]

        display_cut = int(0.35*len(TF_name))
        if display_cut > 31:
            display_cut = 31
        TF_value = TF_value[:display_cut]
        TF_name = TF_name[:display_cut]

        TF_value = np.flip(TF_value)
        TF_name = np.flip(TF_name)

        name_value = []
        for i in range(len(TF_name)):
            tmp = str(TF_name[i]) + " "+ "(" + str(TF_value[i]) + ")"
            name_value.append(tmp)
        name_value = np.array(name_value)

        left, bottom, width, height = 0.8, 0.2, 0.15, 0.7
        ax2 = fig.add_axes([left, bottom, width, height])
        b = ax2.barh(range(len(TF_value)), TF_value,  color = "#6699CC")

        ax2.set_xticks(())
        ax2.set_yticks(range(len(name_value)))
        ax2.set_yticklabels(name_value, fontsize=18)

        

        ax2.set_title("Number of targets\n(GRN by TENET)", fontsize=18)
        figsave_name = filename.replace("txt", "pdf")
        fig.savefig(figsave_name)


    def countOutdegree_single_run(self, file_name, cell_select_file, trajectory_file):
        
        
        ifile = open(file_name)
        TFlist=[];TFlistOutdegree=[]
        for line in ifile:
            temp = line.split()
            if temp[0] not in TFlist:
                TFlist.append(temp[0])
                TFlistOutdegree.append(1)
            else:
                TFlistOutdegree[TFlist.index(temp[0])]=TFlistOutdegree[TFlist.index(temp[0])]+1
        TFlistOutdegreeIndex=numpy.argsort(TFlistOutdegree)
        
        out_folder = ""
        ofile = open(out_folder + file_name+".outdegree.txt","w")
        for i in range(len(TFlist)):
            ofile.write(TFlist[TFlistOutdegreeIndex[-i-1]]+"\t"+str(TFlistOutdegree[TFlistOutdegreeIndex[-i-1]])+"\n")
        ofile.close()
        
        
        filename = file_name+".outdegree.txt"
    
        tmp_cell_select = np.loadtxt(cell_select_file)
        tmp_trajectory = np.loadtxt(trajectory_file)

        TF = np.loadtxt(filename,dtype=object)
        pathway = np.where(tmp_cell_select==1)[0]
        trajectory = tmp_trajectory[pathway]

        fig = plt.figure(None,(20,10))
        left, bottom, width, height = 0.02, 0.15, 0.7, 0.7
        ax1 = fig.add_axes([left, bottom, width, height])
        ax1.quiver(self.embedding[:,0],self.embedding[:,1],self.delta_embedding[:,0],self.delta_embedding[:,1],color='black')
        ax1.scatter(self.embedding[pathway,0],self.embedding[pathway,1],  marker='o',c=trajectory, s=100)
        ax1.axis("off")


        TF_value = TF[:,1].astype(int)
        TF_name = TF[:,0]

        display_cut = int(0.35*len(TF_name))
        if display_cut > 31:
            display_cut = 31
        TF_value = TF_value[:display_cut]
        TF_name = TF_name[:display_cut]

        TF_value = np.flip(TF_value)
        TF_name = np.flip(TF_name)

        name_value = []
        for i in range(len(TF_name)):
            tmp = str(TF_name[i]) + " "+ "(" + str(TF_value[i]) + ")"
            name_value.append(tmp)
        name_value = np.array(name_value)

        left, bottom, width, height = 0.8, 0.2, 0.15, 0.7
        ax2 = fig.add_axes([left, bottom, width, height])
        b = ax2.barh(range(len(TF_value)), TF_value,  color = "#6699CC")

        ax2.set_xticks(())
        ax2.set_yticks(range(len(name_value)))
        ax2.set_yticklabels(name_value, fontsize=18)

        

        ax2.set_title("Number of targets\n(GRN by TENET)", fontsize=18)
        figsave_name = filename.replace("txt", "pdf")
        fig.savefig(figsave_name)
        


    def clean_save_results(self):
        outfolder = "./TE_results"
        
        isExists=os.path.exists(outfolder+'')
        if not isExists:
            os.makedirs(outfolder) 
            os.system("mv TE_result_matrix*sif* ./TE_results/")
            os.system("mv trajectory_*pdf ./TE_results/")
            #os.system("mv cell_select* ./results/")
            os.system("mv TE_result_all*sif* ./TE_results/")
            print("saved in ./ressults")
        else:
            shutil.rmtree(outfolder, ignore_errors=True)
            os.makedirs(outfolder)

            #shutil.move(“oldpos”,”newpos”)
            
            os.system("mv TE_result_matrix*sif* ./TE_results/")
            os.system("mv trajectory_*pdf ./TE_results/")
            #os.system("mv cell_select* ./results/")
            os.system("mv TE_result_all*sif* ./TE_results/")
            print("saved in ./TE_results")
        
        
        
