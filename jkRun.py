# -*- coding: utf-8 -*-
"""
Created on Sun Mar 24 14:05:22 2019

@author: Adrian
"""

from jkNetwork import jkNetwork as jkN
from jkAnalysis import jkAnalysis as jkA

import numpy as np

import os

try:
    task_id_o = int(os.environ['SLURM_ARRAY_TASK_ID']) 
except:
    task_id_o = 20 #(100*5+100)*100+96
    print("Could not load task_id from environment (set task_id={task_id_o}).")

for task_id in range(task_id_o*20,(1+task_id_o)*20):
    m_min=0.2
    m_stepsize=0.05
    grown_per_m=100 #how many grown networks exist for each parameter level
    runs_per_grown=100 #how many runs are made per grown network
    grown_task_id=int(task_id/runs_per_grown) # task_id of grown network which will be the parent
    m=np.round(m_min+int(grown_task_id/grown_per_m)*m_stepsize,2)
    run=grown_task_id % grown_per_m #run number of grown network

    
    if(os.path.isfile(f"jkNetworks/jkN G{str(m).replace('.','_')}-{run}-{task_id - runs_per_grown*grown_task_id}.jk")):
        continue
    try:
        parentNetwork = jkN.loadCompletely(network_id=f"G{str(m).replace('.','_')}-{run}")
        runNetwork = jkN(parent=parentNetwork,Id=f"{task_id - runs_per_grown*grown_task_id}",seed=task_id)
        runNetwork.run_till_catastrophe_and_recovery(thres_p=0.1)
        runNetwork.runAnalysis(thresholds=[0.75,0.5,0.25,0.1])
        #import pdb; pdb.set_trace()
        #jkN.plotNetwork(parentNetwork.C)
        #runNetwork.plotLogs(plotAna=[0,1])
        runNetwork.saveCompletely()
    except Exception as exep:
        print(exep)
