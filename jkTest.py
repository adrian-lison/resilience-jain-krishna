# -*- coding: utf-8 -*-
"""
Created on Thu Sep 20 14:59:00 2018

@author: Adrian
"""

# %% Imports and Functions
from jkNetwork import jkNetwork as jkN
from jkAnalysis import jkAnalysis as jkA

import time as time
#import copy
import numpy as np

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats

def createGrownNetwork(s):
    jkN1 = jkN(s=s,Id="Grown")
    jkN1.run_growth_phase()
    return jkN1

def saveEverything(network,optionalName=None):
    if(optionalName is None): filename = time.strftime("%y_%m_%d %H_%M_%S")
    else: filename = optionalName
    
    network.saveDynamicGraph(optionalName=filename)
    network.saveCurrentNetwork(optionalName=filename)
    network.saveCurrentLog(optionalName=filename)
    return filename

def loadLog_and_Network(s,file):
    jkN1 = jkN(s=s)
    jkN1.loadLog(file)
    jkN1.loadNetwork(file)
    return jkN1

def testDevelopment(initialNetwork,runs,steps_per_run):
    """Runs test iterations starting from an initial network
    
    Parameters:
    initialNetwork (jkNetwork): The initial network to start with
    runs (numeric): Number of simulations to be run
    steps_per_iterations (numeric): Number of time steps for the simulations

    Returns:
        Numpy Array: Array with two columns, one for time to collapse and one
            for integral, each row stands for one simulation
    """
    sample_paths=list()
    for i in range(runs):
        #Initialize with parent network
        testJain=jkN(parent=initialNetwork,Id=i)
        #Iterate fixed number of steps
        testJain.iterate(steps_per_run)
        #Analyse
        testJain.runAnalysis()
        sample_paths.append(testJain)
    return sample_paths


# %% Code

if(1==0):
    '''Testing'''
    # Create Network (and save it)
    createdJain = createGrownNetwork(100)
    createdJain.plotLogs(plotlinks=True, plots1=True,plotLambda1=True)
    
    savedfile = saveEverything(createdJain)
    
    # Load Network
    loadedJain = loadLog_and_Network(100,savedfile)
    loadedJain.plotLogs(plotlinks=True, plots1=True,plotLambda1=True)
    
    # Test Runs
    results=testDevelopment(loadedJain,10,2000)
    
    
    # Test series
    nNetworks=50
    nSpecies=50
    nRuns=100
    nTimeSteps=2000
    dfT2Collapse = pd.DataFrame(index=range(nRuns))
    dfIntegral = pd.DataFrame(index=range(nRuns))
    for i in range(nNetworks):
        jainN = createGrownNetwork(nSpecies)
        savedfile = saveEverything(jainN,f"Test2 Network {i}")
        results=testDevelopment(jainN,nRuns,nTimeSteps)
        dfT2Collapse[i] = results[:,0]
        dfIntegral[i] = results[:,1]
        time.sleep(1)
        
        '''ToDo: progress indication'''
        
    # save 
    dfIntegral.to_csv("Analysis/Test 2 Integral.csv")
    dfT2Collapse.to_csv("Analysis/Test 2 T2Collapse.csv")
    
    
    '''Analysis'''       
    dfIntegral = pd.read_csv("Analysis/Test 1 Integral.csv",index_col=0)
    dfT2Collapse = pd.read_csv("Analysis/Test 1 T2Collapse.csv",index_col=0)
    
    #Analysis of T2Colllapse
    dfT2Collapse = dfT2Collapse.replace(np.infty,2000)
    flattend = dfT2Collapse.values.flatten()
    distribution=sns.distplot(flattend,kde=False)
    
    meanT2Collapse = dfT2Collapse.values.mean(0)
    
    # Analysis of specific fully grown network: collapse times
    CollapseTimes=dfT2Collapse["5"]
    g=sns.distplot(CollapseTimes,kde=False,fit=stats.expon,bins=30,label="Distribution of mean time to collapse")
    
    percs = np.linspace(0,100,21)
    qn_a = np.percentile(np.random.exponential(scale=1,size=1000), percs)
    qn_b = np.percentile(CollapseTimes, percs)
    sns.scatterplot(x=qn_a,y=qn_b)
    x=np.linspace(min(qn_a),max(qn_a),100)
    plt.plot(x,x*(qn_b[5]-qn_b[15])/(qn_a[5]-qn_a[15]), linewidth=2)
    
    sns.scatterplot(x=dfT2Collapse.mean(),y=dfT2Collapse.std())
    plt.xlabel('Mean T2C')
    plt.ylabel('Std T2C')
    np.corrcoef(dfT2Collapse.mean(),dfT2Collapse.std())
    
    np.std(dfT2Collapse.mean())/flattend.mean()
    
    sns.distplot(dfT2Collapse.std())
    
    np.mean(np.std(dfT2Collapse))/np.std(flattend)
    np.mean(dfT2Collapse)/np.mean(flattend)
    
    # Analysis of Integral
    dfIntegral_descr=dfIntegral.describe()
    sns.distplot(dfIntegral_descr.mean())
    
    # Correlation
    sns.scatterplot(x=dfT2Collapse.mean(),y=dfIntegral_descr.mean())
    plt.xlabel('Mean T2C')
    plt.ylabel('Mean Integral')
    np.corrcoef(x=dfT2Collapse.mean(),y=dfIntegral_descr.mean())
    
    # Load Initial Networks
    nNetworks=50
    nSpecies=50
    lNetworks=list()
    for i in range(nNetworks):
        lNetworks.append(loadLog_and_Network(nSpecies,f"Test1 Network {i}"))
    
    lLinks=[nw.linkLog[-1] for nw in lNetworks]
    
    
    sns.scatterplot(x=lLinks,y=dfT2Collapse.mean())
    np.std(lLinks)/np.mean(lLinks)
    np.std(dfT2Collapse.mean())/np.mean(dfT2Collapse.mean())
    
    for s in range(50,100):
        mat=np.ones((s,s),dtype=bool)
        np.fill_diagonal(mat,False)
        print(f"{s}: {np.round(jkN.getPFE(mat)[1],2)}")
    
    #testing
    
    Jain3=jkN(parent=createGrownNetwork(50),Id=1)
    Jain3.run_till_catastrophe_and_recovery(thres_p=0.1)
    Jain3.runAnalysis([40,20,5])  
    Jain3.plotLogs(plotS1=True,plotLinks=True,plotLambda1=True,plotAna=[0,1])   
     
    Jain3.saveCompletely()
    
    Jain4 = jkN.loadCompletely("jkN Grown-1")
    
    
    jkN1 = jkN(s=100,m=0.25,Id="Grown")
    jkN1.run_growth_phase()
    jkN1.plotLogs(plotS1=True,plotLinks=True,plotLambda1=True)  
    
    jkN1.iterate(400)
    jkN.plotNetwork(jkN1.C)
    
    jkN2=jkN(parent=jkN1,Id=1,seed=42)
    jkN2.iterate(500)  
    
    jkN3=jkN(parent=jkN1,Id=1,seed=42) 
    jkN3.iterate(500)
    
    np.all(jkN2.C==jkN3.C)
    
    #kN2.countlogs==jkN3.countlogs
    np.all(jkN2.C==jkN3.C)
    np.all(jkN2.restoreGraph(20)==jkN3.restoreGraph(20))
    
    np.where(jkN2.restoreGraph(-1)!=jkN2.C)

    jkN2.graphHist==jkN3.graphHist
    jkN2.linkLog==jkN3.linkLog
    
    
    
    
    ##
    from jkNetwork import jkNetwork as jkN
    jkN.plotNetwork(jkN1.C)
    
    
    
    jkN2.iterate(1) 
    jkN3.iterate(1)
    
    jkN2.plotLogs(plotS1=True,plotLinks=True,plotLambda1=True)
    jkN3.plotLogs(plotS1=True,plotLinks=True,plotLambda1=True) 
    
    jkN2=jkN(parent=jkN1,Id=1)
    jkN2.iterate(steps=10000)
    jkN2.runAnalysis(thresholds=[0.5,0.1])
    jkN2.analysis=jkA(jkN2,thresholds=[0.5,0.1])
    j.plotLogs(plotS1=True,plotLinks=True,plotLambda1=True,plotAna=[0,1])   
     
    jkN2.saveCompletely()
    jkN2.saveCurrentNetwork("testA")
    
    jkN2.analysis.save()
    
    jkN2.saveCurrentLog("testB",withGraphHist=False)

    '''Other testing stuff
    #savedfile= '18_11_10 13_24_47'
      
    len(loadedJain.s1Log)
    
    startJain=copy.deepcopy(loadedJain)
    startJain.run_till_catastrophe()
    len(startJain.s1Log)  
    
    for i in range(500):
        loadedJain.iterate()
    loadedJain.plotLogs(plotlinks=True, plots1=True,plotLambda1=True)
    
    s1Log=loadedJain.s1Log'''
