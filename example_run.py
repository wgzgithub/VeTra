#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 24 21:56:47 2020

@author: gzw
"""

import VeTra as vt



#initialize an example 
ex1 = vt.VeTra("embedding.txt", "delta_embedding.txt")
#run VeTra to get trajectories with pseudotime
#All files all saved in folder "vetra/results"
ex1.vetra(12, 5, 3)
#run TENET for all trajectories
#ex1.run_tenet_tf("chroman_exp_filtered.csv", 15, 1, 'mouse', bulk_run=True)
# run downstream analysis
#usuer can choose the methods and species here
#ex1.makeGRN_tf("links",1000, bulk_run=True)




