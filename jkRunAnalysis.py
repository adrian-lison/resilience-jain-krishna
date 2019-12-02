# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 09:22:04 2019

@author: Adrian
"""

from jkNetwork import jkNetwork as jkN
from jkAnalysis import jkAnalysis as jkA

import numpy as np
import pandas as pd

import os
#
#try:
#    task_id_o = int(os.environ['SLURM_ARRAY_TASK_ID']) 
#except:
#    task_id_o = 0#(100*1+28)*5
#    print(f"Could not load task_id from environment (set task_id={task_id_o}).")
#    
results = list()

for task_id in range(0,50000): #range(task_id_o*20,(1+task_id_o)*20):
    m_min=0.2
    m_stepsize=0.05
    grown_per_m=100 #how many grown networks exist for each parameter level
    runs_per_grown=100 #how many runs are made per grown network
    grown_task_id=int(task_id/runs_per_grown) # task_id of grown network which will be the parent
    m=np.round(m_min+int(grown_task_id/grown_per_m)*m_stepsize,2)
    run=grown_task_id % grown_per_m #run number of grown network
    
    if((task_id % 100) == 0): print(task_id)
    if(not os.path.isfile(f"jkNetworks/jkN G{str(m).replace('.','_')}-{run}-{task_id - runs_per_grown*grown_task_id}.jk")):
        #print(f"Network G{str(m).replace('.','_')}-{run}-{task_id - runs_per_grown*grown_task_id} not found.")
        continue
    
    analysisNetwork = jkN.loadCompletely(network_id=f"G{str(m).replace('.','_')}-{run}-{task_id - runs_per_grown*grown_task_id}")
         
    ana=analysisNetwork.analysis
    for i,threshold in enumerate(ana.thresholds):
        results.append({"parent":f"G{str(m).replace('.','_')}-{run}",
                  "run":task_id - runs_per_grown*grown_task_id,
                  "threshold":ana.thresholds[i],
                  "tTresh":ana.tTresh[i],
                  "tRecov":ana.tRecov[i],
                  "tMin":ana.tMin[i],
                  "lifetime":ana.lifetime[i],
                  "recov_time":ana.recov_time[i],
                  "tPrecol":ana.tPrecol[i],
                  "tReorg":ana.tReorg[i],
                  "s1Itg":ana.s1Itg[i]
                 })

df=pd.DataFrame(results,columns=results[0].keys())
df.to_csv("Analysis/runAna.csv")
