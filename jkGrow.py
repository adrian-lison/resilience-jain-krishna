# -*- coding: utf-8 -*-
"""
Created on Sat Mar 23 11:39:22 2019

@author: Adrian
run 0-10999
"""

from jkNetwork import jkNetwork as jkN

import numpy as np

import os

try:
    task_id = int(os.environ['SLURM_ARRAY_TASK_ID']) 
except:
    print("Could not load task_id from environment (set task_id=0).")
    task_id = 0

m_min=0.2
m_stepsize=0.05
runs_per_m=100
m=np.round(m_min+(int(task_id/runs_per_m)*m_stepsize),2)
run=task_id % runs_per_m   

jkN1 = jkN(s=100,m=m,Id=f"G{str(m).replace('.','_')}-{run}",seed=task_id)
jkN1.run_growth_phase()
jkN1.saveCompletely()

