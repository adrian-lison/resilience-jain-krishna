# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 11:11:25 2019

@author: Adrian
"""

from jkNetwork import jkNetwork as jkN
from jkTopology import jkTopology as jkT

import numpy as np
import pandas as pd
        
m_min=0.2
m_stepsize=0.05
runs_per_m=100

results=list()

for task_id in range(0,500):
    if(task_id % 50 ==0): print(task_id)
    m=np.round(m_min+(int(task_id/runs_per_m)*m_stepsize),2)
    run=task_id % runs_per_m   
    try:
        jkGrown = jkN.loadCompletely(network_id=f"G{str(m).replace('.','_')}-{run}")
    except Exception as excep:
        print(excep)
        continue
    
    jkTop=jkT(jkGrown)
    results.append({"network":jkGrown.id,
              "nw_m":jkTop.m,
              "nw_nodes":jkTop.s,
              "lambda1":jkTop.lambda1,
              "nw_edges":jkTop.nw_edges,              
              "nw_comple":jkTop.nw_comple,
              "nw_between_central_mean":jkTop.nw_between_central_mean,
              "nw_between_central_max":jkTop.nw_between_central_max,
              "nw_between_central_skew":jkTop.nw_between_central_skew,
              "core_nodes":len(jkTop.core_nodes),
              "core_edges":jkTop.core_edges,
              "core_comple":jkTop.core_comple,
              "core_between_central_mean":jkTop.core_between_central_mean,
              "core_between_central_max":jkTop.core_between_central_max,
              "core_between_central_skew":jkTop.core_between_central_skew,
              "core_subgraphs":len(jkTop.core_bs_nodes),
              "av_core_path_mult":jkTop.av_core_path_mult,
              "av_core_path_mult_complete":jkTop.av_core_path_mult_complete,
              "av_core_cycles":jkTop.av_core_cyclelen
                 })
    
df=pd.DataFrame(results,columns=results[0].keys())
df.to_csv(f"Analysis/growAna.csv")
