# -*- coding: utf-8 -*-
"""
Created on Sat Nov 10 13:47:38 2018

@author: Sascha
"""

import numpy as np
import scipy as sp
import scipy.linalg as la
import scipy.integrate as itg
import pickle
import time
import networkx as nx

class jkAnalysis:
    
    def __init__(self,network,thresholds=[]):
        if(len(thresholds)==0):
            thresholds=[0.5*network.s,0.25*network.s,0.1*network.s]
        else:
            thresholds=list(np.array(thresholds)*network.s)
        Nthresholds = len(thresholds)
        self.thresholds=thresholds
        self.networkId=network.id
        self.starttime=0 #preliminary
        
        self.nStart = network.restoreGraph(self.starttime)
        s1=np.array(network.s1Log)
        
        self.org_phase = np.where(s1>=network.s-1)[0] #all times at which the network is above the threshold
        
        self.above_threshold=list()
        self.tTresh=jkAnalysis.naArray(Nthresholds) #time at which the treshold is undershoot for the first time
        self.nTresh=list()
        self.tRecov=jkAnalysis.naArray(Nthresholds) #time at which the ttresholds is reached again
        self.nRecov=list()
        self.tMin=jkAnalysis.naArray(Nthresholds)
        self.nMin=list()
        self.lifetime=jkAnalysis.naArray(Nthresholds)
        self.recov_time=jkAnalysis.naArray(Nthresholds)
        self.tPrecol=jkAnalysis.naArray(Nthresholds)
        self.nPrecol=list()
        self.tReorg=jkAnalysis.naArray(Nthresholds)
        self.nReorg=list()
        self.s1Itg=jkAnalysis.naArray(Nthresholds)    

        for i,threshold in enumerate(thresholds):
            if not np.any(s1<=threshold):
                self.above_threshold.append(None)
                self.s1Itg[i]=self.s1Integral(len(s1),s1)
                continue
            
            #all points in time at which network is above threshold
            self.above_threshold.append(np.where(s1>threshold)[0]) 
            
            #point in time when threshold is undershot the first time
            self.tTresh[i]=np.where(s1<=threshold)[0][0]
            self.nTresh.append(network.restoreGraph(self.tTresh[i]))
            
            #time until collapse/undershoot
            self.lifetime[i]=self.tTresh[i]-self.starttime
            
            #integral over time
            self.s1Itg[i]=self.s1Integral(self.tTresh[i],s1)
            
            #point in time when network is fully organized the last time before collapse
            self.tPrecol[i]=self.findTprecol(self.tTresh[i])
            self.nPrecol.append(network.restoreGraph(self.tPrecol[i]))
            
            #point in time when network is recovered (above threshold again)
            self.tRecov[i]=self.findRecov(self.tTresh[i],self.above_threshold[i])
            if(np.isnan(self.tRecov[i])): continue
            self.nRecov.append(network.restoreGraph(self.tRecov[i]))
            
            #time needed until revocered (above threshold again)
            self.recov_time[i]=self.tRecov[i]-self.tTresh[i]
            
            #point in time when network is fully reorganized again
            self.tReorg[i]=self.findReorg(self.tTresh[i])
            if(np.isnan(self.tReorg[i])): continue
            self.nReorg.append(network.restoreGraph(self.tReorg[i]))
                       
            #point in time when population minimum is reached
            self.tMin[i]=self.findMin(self.tPrecol[i],self.tRecov[i],s1)
            self.nMin.append(network.restoreGraph(self.tMin[i]))
            
    
    @staticmethod        
    def naArray(dim):
        a = np.empty(dim)
        a[:] = np.nan
        return a

    def s1Integral(self,Ttresh,s1):
        #import pdb; pdb.set_trace()
        s1=np.array(s1[self.starttime:self.starttime+int(Ttresh)]) 
        return itg.simps(s1)  
        
    def findTprecol(self,Ttresh):     
        """Finds point in time when network is fully organized the last time before collapse"""
        #remove everything that occurs after Ttresh
        org_before = np.delete(self.org_phase,np.where(self.org_phase>=Ttresh))    
        return np.amax(org_before)
    
    def findRecov(self,Ttresh, above_threshold):
        """Finds point in time when treshold is reached again after undershoot"""
        above_after = np.delete(above_threshold,np.where(above_threshold<Ttresh)) #remove everything that occurs bevfore Ttresh
        if(len(above_after)==0): return np.nan
        return np.amin(above_after)
    
    def findReorg(self,Ttresh):
        """Finds point in time when network is fully reorganized again"""
        org_after=np.delete(self.org_phase,np.where(self.org_phase<=Ttresh))
        if(len(org_after)==0): return np.nan
        return np.amin(org_after)
    
    def findMin(self,Tprecol,Trecov,s1):
        """Finds point in time when population minimum is reached"""
        
        collapse=s1[Tprecol.astype(int):Trecov.astype(int)]
        return np.argmin(collapse)+Tprecol.astype(int)
    
    def save(self,optionalName=None):        
        if(optionalName is not None): filename = optionalName + ".ana"
        if(self.networkId is not None): filename = f"jkA {self.networkId}" + ".ana"      
        else: filename = "jkA " + time.strftime("%y_%m_%d %H_%M_%S") + ".ana" 
        ana_file = open('Analysis/'+filename, 'wb')
        pickle.dump(self,ana_file)
        ana_file.close()