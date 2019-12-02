# -*- coding: utf-8 -*-
"""
Created on Thu Sep 20 12:07:52 2018

@author: Adrian
"""
import pylab as pl
#import networkx as nx
import numpy as np
import scipy.linalg as la
import scipy.integrate as itg
import scipy.stats as st
import time as time

import pickle

from gexfDynamic import DynamicGraph
from jkAnalysis import jkAnalysis as jkA
from jkTopology import jkTopology as jkT
import networkx as nx

import plotly.graph_objs as go
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
  

class jkNetwork:
	'''Class to simulate networks according to the Jain-Krishna model.
	
	If parent network is given, the topology is copied so that the
	simulation can be continued from the last state.
	
	Attributes
	----------
	
	id : str
		Unique identifier of network.
	s : int
		Number of nodes/species.
	m : float
		Expected percentage of links formed during mutation.
	p : float
		Probability of formation for a link.
	timestep :
		Current timestep of simulation.
	C : ndarray(boolean,boolean)
		Adjacency matrix.
	x : ndarray(float)
		Current vector of relative population sizes for each species.
	x_old : ndarray(float)
		Previous step vector of relative population sizes for each species.
	PFE : ndarray(float)
		Current Perron-Frobeniues Eigenvector of the graph.
	lambda1 : float
		Largest eigenvalue of adjacency matrix of graph.

	analysis : jkAnalysis
		Object which computes measures of resilience for the simulation run.
	
	countlogs : int
		Current number of logged time steps.
	linkLog : list of int
		Log of number of links in the graph at each step.
	maxLog : list of float
		Log of relative size of biggest species at each step.
	s1Log : list of int
		Log of number of species with size bigger than zero at each step.
	complLog : list of float
		Log of complexity measure of the graph at each step.
	lambdaLog : list of float
		Log of largest eigenvalue of adjacency matrix of the graph at each step.
	graphHist : dict
		Log of evolution of links in the graph over time.
	logPopulation : boolean
		Log population over time. 
	popHist : dict
		Log of population vector at each each time step.
	
	
	Parameters
	----------
	s : int, optional
		Number of nodes/species.
	m : float, optional
		Expected percentage of links formed during mutation.
	logPop : boolean, optional
		Log the population time at each time step.
	parent : jkNetwork, optional
		A network from which the topology should be inherited.
	id : str, optional
		Unique identifier for this network. If parent network exists, id is derived from parent using concatenation.
	seed : float, optional
		Seed for initialization of random generator.
	pCollapsePars : (boolean,float,int)
		Parameters for computation of p(collapse in next time step) - (log or not?,threshold,samples)
	'''


    def __init__(self,s=50,m=0.25,logPop=False,parent=None,id=None,seed=None,pCollapsePars=(False,0,0)):
		
        if(seed is not None): np.random.seed(seed)
        #print(np.random.get_state())
        
        #Number of different species
        if parent:
            self.s=parent.s # number of nodes/species; same as parent network
            self.m=parent.m # rate of links of species after mutation; same as parent network
            self.id=f"{parent.id}-{id}" # if parent network exists, derive id from parent
        else:
            self.s=s # number of nodes/species
            self.m=m # rate of links of species after mutation
            self.id=id # identifier of network
            
        self.p=self.m/(self.s-1) # formation probability for each new link
        
        self.logPopulation=logPop # should the population be recoded over time?
        self.initialize_Log()
        self.initialize(parent)
    
        self.t_wait = 5000 # time to wait before mutation
        self.x_old = None
		
		self.pCollapsePars = pCollapsePars
        
        
    """___________________________________________________________"""
    """Functions for Jain-Krishna Model"""
    
    def initialize(self,parent):
        """Initialize network"""
        
        self.timestep = -1
        self.countlogs = 0
        
        #Set C to adjacency matrix of erdos-renyi graph
        if parent:
            self.C=parent.C.copy()
            self.x=parent.x.copy()
            self.PFE=parent.PFE.copy()
            self.lambda1=parent.lambda1.copy()
        else:
            self.C=np.random.choice(a=[False,True],size=(self.s,self.s),p=[1-self.p,self.p])
            np.fill_diagonal(self.C,False)
        
            #Set x to initial population values (all size of 1/s)
            self.x = np.array([1/self.s]*self.s)
            
            self.PFE, self.lambda1 = jkNetwork.getPFE(self.C)
            self.updatePopulationsPFE()  
        
        # Log nodes
        if(self.logPopulation): self.logPopDev(self.x,np.zeros(self.s))
            
        # Log the original edges
        for n in range(self.s):
            rowDiff = self.C[n,:]
            columnDiff = np.zeros(1,dtype=bool)
            #import pdb; pdb.set_trace()
            self.logGraphDev(n, rowDiff, columnDiff)
        
        self.log_iteration()
            
    def iterate(self,steps=1):
        for step in range(steps):
            self.timestep += 1
        
            # Let a random species from the set of least fit species mutate        
            setOfLeastFitNodes = self.getSetOfLeastFit()
            while(len(setOfLeastFitNodes)==0):
                setOfLeastFitNodes=self.getSetOfLeastFit()
             
            self.mutatingSpecies = np.random.choice(setOfLeastFitNodes) 
            self.mutateSpecies(self.mutatingSpecies)
            
            #Update PFE
            self.PFE, self.lambda1 = jkNetwork.getPFE(self.C) 
            # Let the species grow according to the growth rule
            self.updatePopulationsPFE()
            
            self.log_iteration()
			
			#log collapse probability
			if(pCollapsePars[0]):
				self.pCollapseLog.append(self.getCollapse_prob(pCollapsePars[1],pCollapsePars[2]))

	def iterate_species(self,species):
		self.timestep += 1
		 
		self.mutatingSpecies = species 
		self.mutateSpecies(self.mutatingSpecies)
		
		#Update PFE
		self.PFE, self.lambda1 = jkNetwork.getPFE(self.C) 
		# Let the species grow according to the growth rule
		self.updatePopulationsPFE()
		
		self.log_iteration()     
        
    def step_back(self):
        self.timestep -= 1
        
        # Take back species growth
        x_new = self.x.copy()
        self.x = self.x_old
        if(self.logPopulation): self.unlogPopDev(x_new, self.x_old)
    
        # Take back last mutation of species            
        self.unmutateSpecies(self.mutatingSpecies)
        
        self.unlog_iteration()
		
	def sample_next_pop(self,species):
		self.iterate_species(species)
		pop = self.x.copy()
		self.step_back()
		return(pop)
		
	def getCollapse_prob(self,threshold,samples):
		'''Approximate the probability of a collapse in the next time step'''
		setOfLeastFitNodes = self.getSetOfLeastFit()
		while(len(setOfLeastFitNodes)==0):
			setOfLeastFitNodes=self.getSetOfLeastFit()
    	
		# Compute the overall probability of collapse as the average of the collapse probabilities of all species in the
		# set of least fit nodes. To do so, sample from the population distribution of the next time step 
		# and estimate the probability of collapse for each species separately.
		return(mean(
			[sum(
				[len(np.where(sample_next_pop(species)>0.001)[0])<=threshold*self.s for i in range(samples)]
				)/samples 
				for species in setOfLeastFitNodes]))
		

    def fast_td(self,x, t0): #fast time development
        return np.dot(self.C,x)-x*np.sum(np.dot(self.C,x))
    
    def updatePopulations(self):
        '''Changes the relative population of every species according to the growth rule (catalytic growth)'''
      
        self.x_old = self.x #save populations before update
         
        sol=itg.odeint(self.fast_td,self.x,np.array([0,self.t_wait])) 
        self.x=sol[1]
        
        #populations = list()
        #y = np.dot(np.linalg.matrix_power(np.identity(s)+C,t_wait),y)
        # Calculate new relative size
        #y = np.divide(y,np.sum(y)) 
        
        if(self.logPopulation):        
            self.logPopDev(self.x,self.x_old)
    
    def updatePopulationsPFE(self):
        '''Changes the relative population of every species according to the growth rule (catalytic growth)'''
        
        self.x_old = self.x #save populations before update
        
        self.x=self.PFE/np.sum(self.PFE)
                
        if(self.logPopulation):        
            self.logPopDev(self.x,self.x_old)              
        return(True)
        
    def getPFE(C):
        '''Restore the relative population of every species according to the growth rule (catalytic growth)'''
   
        #get eigen value of C with largest real part, per theorem, this is real and >=0
        evalues, evectors = la.eig(C) #get array of eigen values and vectors
        #print("eigen values", evalues)
        lamb1 = np.amax(evalues.real)
        #print("lamb1: "+str(lamb1))
        if(np.around(lamb1.imag,5)!=0 or lamb1.real<0):
            raise Exception("Non-zero imaginary part or negative real part of lambda1.")
        
        #get all eigenvectors, for the case lambda1 is degenerate
        #list of indices of eigen vectors t lambda1
        min_ind = [i for i,ev in enumerate(evalues) if (np.around(ev,5)==np.around(lamb1,5))]
        #print("numb. of PFEs: ",len(min_ind))
        
        #find eigenvector with no negative entries
        for j in min_ind:
            if((np.all(np.around(x_i,5)>=0 and np.around(x_i.imag,5)==0) for x_i in evectors[:,j])):
                x_lamb=evectors[:,j]
                #print("y_lambd: "+ str(y_lamb))
                return(x_lamb.real/np.sum(x_lamb.real),evalues[j].real)
        print(evectors)
        raise Exception("No eigenvector without negative entries.") #no such eigen vector found
           
    def getSetOfLeastFit(self):
        '''Returns the set of the least fit species, that is: the ones with the smallest population'''
        return np.where(self.x==np.amin(self.x))[0]
        
    def getSetOfLeastFit2(self):
        '''Returns the set of the least fit species, randomly chosen from all species with probabilities proportional to population'''   
        rand = np.random.random_sample(size=len(self.x))
        return np.where(self.x < rand)[0]
    
    def mutateSpecies(self,n,log=True):
        '''Mutates a species: removes all links of the species and randomly creates new links'''
        
        self.mutatingColumn = self.C[:,n].copy()
        self.mutatingRow = self.C[n,:].copy()
        
        # Remove all existing relationships of species n
        # Randomly create new relationships of species n
        self.C[:,n]=np.random.choice([True,False],p=[self.p,1-self.p],size=self.s)
        self.C[n,:]=np.random.choice([True,False],p=[self.p,1-self.p],size=self.s)  
        self.C[n,n] = False
    
        if(log):
            #import pdb; pdb.set_trace()
            rowDiff = np.bitwise_xor(self.mutatingRow,self.C[n,:])
            columnDiff = np.bitwise_xor(self.mutatingColumn,self.C[:,n])
            
            #save time graph changes for time development
            self.logGraphDev(n, rowDiff, columnDiff)
            
    def unmutateSpecies(self,n,log=True):
        newColumn = self.C[:,n].copy()
        newRow = self.C[n,:].copy()
        self.C[:,n] = self.mutatingColumn
        self.C[n,:] = self.mutatingRow
    
        if(log):
            rowDiff = self.mutatingRow - newRow
            columnDiff = self.mutatingColumn - newColumn
            
            #save time graph changes for time development
            self.unlogGraphDev(n, rowDiff, columnDiff)
            
    """___________________________________________________________"""        
    """Auto-Run"""
    def run_random_phase(self):
        #todo lambda1
        while self.lambdaLog[-1]==0:
            self.iterate()
        
    def run_growth_phase(self):
        self.iterate()
        while (self.s1Log[-1]<self.s):
            self.iterate()
       
    def run_till_catastrophe(self,thres_p=None,max_steps=10000):
        steps=1
        self.iterate()
        threshold = self.s-1 if thres_p is None else self.s*thres_p
        #iterate until threshold is reached        
        while (self.s1Log[-1]>=threshold and steps<max_steps):
            self.iterate()
            steps+=1
        return steps
            
    def run_till_recovery(self,max_steps=10000):
        steps=1
        self.iterate()
        #iterate until fully recovered again
        while (self.s1Log[-1]<self.s-1 and steps<max_steps):
            self.iterate()
            steps+=1
        return steps
            
    def run_till_catastrophe_and_recovery(self,thres_p=None,max_steps=10000):
        steps_c=self.run_till_catastrophe(thres_p=thres_p,max_steps=max_steps)
        self.run_till_recovery(max_steps=max_steps-steps_c)
                    
    """___________________________________________________________"""
    """Stats"""
   
    def getLength(self):
        return len(self.s1Log)         
    
    def countEdges(self):
        '''counts number of edges'''
        return (np.sum(self.C))
       
    def getS1(self):
        '''returns the number of species that is not zero'''   
        return len(np.where(self.x>0.001)[0])
       
    def getDegrees(self,C):
        degrees = np.zeros(self.s)
        for i in range(self.s):
            count1 = np.count_nonzero(C[:,i])
            count2 = np.count_nonzero(C[i,:])
            degrees[i] = count1+count2
        return degrees
          
    def complexity(self,C):
        """return offdiagonal complexity [Claussen 2007]"""
        c = np.zeros((self.s, self.s))
        degrees = self.getDegrees(C) #array of degrees
        #ks = G.degree()  # degree sequence
        for i in range(self.s):
            for j in range(self.s):
                if(C[i,j]==1 and degrees[i]<=degrees[j]):
                    c[int(degrees[i]), int(degrees[j])] += 1
        atilde = np.zeros(self.s)
        for i in range(self.s):
            atilde[i] = np.sum(np.diag(c, i))
        a = atilde / np.sum(atilde)
        compl = st.entropy(a)
        return compl
    
    def runAnalysis(self,thresholds=[]):
        self.analysis=jkA(self,thresholds)        
 
    """__________________________________________________________"""    
    """Logging"""
    
    def initialize_Log(self):
         # List that tracks the number of links in the graph for every step
        self.linkLog = list()
        # List that tracks the biggest population in the graph for every step
        self.maxLog = list()
        # List that tracks the number of species with population bigger than zero
        self.s1Log = list()
        #List the (offdiagonal) complexity for each step
        self.complLog=list()
        #List of the lambda1 values for each step
        self.lambdaLog=list()
		#List of the collapse probabilities for each step
        self.pCollapseLog=list()
        
        #Log the development of the graph, i.e. link changes
        #id for each link from node i to node j: s*i+j
        self.graphHist=dict()

        #Log the development of the populations
        if(self.logPopulation):
            self.popHist=dict()        
            for i in range(self.s):
                self.popHist[i]=list()
        else:
            self.popHist=None
            
    def logGraphDev(self, mutatingspecies, row_diff, column_diff):
        '''log edge changes'''
    
        rowchanges = np.nonzero(row_diff)[0]
        columnchanges = np.nonzero(column_diff)[0]
        
        for i in range(len(rowchanges)):
            changed_id=(self.s*mutatingspecies+rowchanges[i])
            if changed_id not in self.graphHist: self.graphHist[changed_id] = [np.int16(self.timestep)]
            else: self.graphHist[changed_id].append(np.int16(self.timestep))
            
        for i in range(len(columnchanges)):
            changed_id=(self.s*columnchanges[i]+mutatingspecies)
            if changed_id not in self.graphHist: self.graphHist[changed_id] = [np.int16(self.timestep)]
            else: self.graphHist[changed_id].append(np.int16(self.timestep))
        
    def unlogGraphDev(self, mutatingspecies, row_diff, column_diff,printit=False):
        '''unlog edge changes'''
    
        rowchanges = np.nonzero(row_diff)[0]
        columnchanges = np.nonzero(column_diff)[0]
        
        for i in range(len(rowchanges)):
            changed_id=(self.s*mutatingspecies+rowchanges[i])
            del self.graphHist[changed_id][-1]
        for i in range(len(columnchanges)):
            changed_id=(self.s*columnchanges[i]+mutatingspecies)
            del self.graphHist[changed_id][-1]
       
    def logPopDev(self, x, x_old):
        '''log population changes'''
        changes=np.nonzero(x-x_old)
        for i in range(len(changes[0])):
            changed_id=changes[0][i]
            self.popHist[changed_id].append(x[changed_id])
            
    def unlogPopDev(self, x, x_old):
        '''unlog population changes'''
        changes=np.nonzero(x-x_old)
        for i in range(len(changes[0])):
            changed_id=changes[0][i]
            del self.popHist[changed_id][-1]
    
    def log_iteration(self):
        # Log the number of links in the graph
        self.linkLog.append(np.int16(self.countEdges()))
        
        # Log the relative size of the biggest population in the graph
        self.maxLog.append(np.amax(self.x))
        
        # Log the number of species with population bigger than zero
        self.s1Log.append(np.int16(self.getS1()))
        
        #log lambda1
        self.lambdaLog.append(np.float16(self.lambda1.real))
        
        #Log complexity
        #self.complLog.append(complexity(C))
        
        self.countlogs+=1
        
    def unlog_iteration(self):
        # Log the number of links in the graph
        del self.linkLog[-1]
        
        # Log the relative size of the biggest population in the graph
        del self.maxLog[-1]
        
        # Log the number of species with population bigger than zero
        del self.s1Log[-1]
        
        #Log complexity
        #self.complLog.append(complexity(C))
        self.countlogs-=1

    def plotLogs(self,plotly=True,plotS1=True,plotLinks=True,plotComplexity=False,plotLambda1=True,plotPCollapse=True,plotAna=None):
        # Display the evolution of the number of links throughout the simulation
        #plotfig=pl.figure(figsize=(25, 25))
        #pl.subplot(4,4,1)
        #pl.plot(linkLog, label="Run 1")
        #pl.xlabel('n')
        #pl.ylabel('links')
        #pl.title('Number of Links in Graph')
        #pl.legend()
        if(plotly):
            traces=list()
            if(plotS1): traces.append(go.Scatter(y=self.s1Log,name="s1",mode="lines",line=dict(shape="hv")))
            if(plotLinks): traces.append(go.Scatter(y=self.linkLog,name="links",mode="lines",line=dict(shape="hv")))
            if(plotComplexity): traces.append(go.Scatter(y=self.complLog,name="Offdiagonal complexity",yaxis='y2',mode="lines",line=dict(shape="hv")))
            if(plotLambda1): traces.append(go.Scatter(y=self.lambdaLog,name="lambda1",yaxis='y2',mode="lines",line=dict(shape="hv")))
			if(plotPCollapse): traces.append(go.Scatter(y=self.pCollapseLog,name="pCollapse",yaxis='y2',mode="lines",line=dict(shape="hv")))
              
            layout = go.Layout(title=f"Evolution of Network {self.id} (s={self.s},m={self.m}%)",
                               xaxis=dict(title="Time step"))
            
            if(plotLambda1 or plotComplexity):
                layout.update(yaxis2=dict(
                    overlaying='y',
                    side='right'))
                
            if(plotAna is not None):
                shapelist=list()
                
                for i in plotAna:                    
                    shapelist.append({
                        'type': 'rect',
                        # x-reference is assigned to the x-values
                        'xref': 'x',
                        # y-reference is assigned to the plot paper [0,1]
                        'yref': 'paper',
                        'x0': str(self.analysis.tTresh[plotAna[i]]),
                        'y0': 0,
                        'x1': str(self.analysis.tRecov[plotAna[i]]),
                        'y1': 1,
                        'fillcolor': '#a6a6a6',
                        'opacity': 0.3,
                        'layer':"below",
                        'line': {
                            'width': 0,
                        }})
    
                    shapelist.append({
                        'type': 'line',
                        'x0': str(0),
                        'y0': str(self.analysis.thresholds[plotAna[i]]),
                        'x1': str(len(self.s1Log)),
                        'y1': str(self.analysis.thresholds[plotAna[i]]),
                        'line': {
                            'color': 'rgb(255, 77, 77)',
                            'width': 2,
                            'dash': 'dashdot',
                        }})
                
                layout.update({'shapes': shapelist})
                        
            fig = dict(data=traces, layout=layout)
            plot(fig)
    
        else:
            ax1=pl.subplot()
            if(plotS1): ax1.plot(self.s1Log, color='b',label="s1")
            if(plotLinks): ax1.plot(self.linkLog, color='r', label="links")
            ax2=ax1.twinx()
            if(plotComplexity): ax2.plot(self.complLog,color='y', label="Offdiagonal complexity")
            if(plotLambda1): ax2.plot(self.lambdaLog, color="g", label="lambda1")
            ax1.legend()
            #ax2.legend()
            pl.show()
            
            return pl
        
    def plotNetwork(matrix):
        pfe, lambda1=jkNetwork.getPFE(matrix)
        G=nx.DiGraph(np.transpose(matrix))
        core, coreparts, pheriphery = jkT.findCoreAndPheriphery(G, matrix, pfe, lambda1)     
                
        nodesizes=(1+pfe*5)*20
        nodesizedict={key:val for key,val in enumerate(nodesizes)}
        
        nx.set_node_attributes(G, name="size", values=nodesizedict)
        
        origpos={k: (v*4) for k, v in nx.spring_layout(G.subgraph(core)).items()}     
        #import pdb; pdb.set_trace()
        pos=nx.spring_layout(G,k=1/np.sqrt(len(G.nodes)),pos=origpos,fixed=core)
                            
        edgedict=[dict(ax=pos[edge[0]][0], ay=pos[edge[0]][1], axref='x', ayref='y',
                        x=pos[edge[1]][0], y=pos[edge[1]][1], xref='x', yref='y',
                        opacity=0.7,arrowwidth=[0.4,0.4,1,1,2][(3 if edge[0] in core else 0)+(1 if edge[1] in core else 0)+(2 if edge[0] in pheriphery else 0)],
                        arrowcolor=["grey","grey","blue","blue","red"][(3 if edge[0] in core else 0)+(1 if edge[1] in core else 0)+(2 if edge[0] in pheriphery else 0)],
                        startstandoff=G.nodes[edge[0]]["size"]/2,standoff=G.nodes[edge[1]]["size"]/2) for edge in G.edges()]
        
    
        tracecolors=[["yellow","blue","red"][(2 if node in core else 0)+(1 if node in pheriphery else 0)]for node in G.nodes()]

        node_trace = go.Scatter(
            x=[],
            y=[],
            text=[],
            mode='markers',
            hoverinfo='text',
            marker=dict(
                showscale=False,
                color=tracecolors,
                size=nodesizes,
                line=dict(width=1)))
                       
        for node in G.nodes():
            x, y = pos[node]
            node_trace['x'] += tuple([x])
            node_trace['y'] += tuple([y])
            node_trace['text']+=tuple([node])
            
        #import pdb; pdb.set_trace()
            
        fig = go.Figure(data=[node_trace],
             layout=go.Layout(
                title='<br>Jain-Krishna Network Graph',
                titlefont=dict(size=16),
                showlegend=False,
                hovermode='closest',
                margin=dict(b=20,l=5,r=5,t=40),
                annotations = edgedict,
                xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                yaxis=dict(showgrid=False, zeroline=False, showticklabels=False)))

        plot(fig)

    """___________________________________________________________"""
    """Store and Load Functionality"""
    
    def saveCurrentNetwork(self,optionalName=None):
        if(optionalName is None): filename = time.strftime("%y_%m_%d %H_%M_%S") + ".network"
        else: filename = optionalName + ".network"      
        network_file = open('Networks/'+filename, 'wb')
        pickle.dump(self.C, network_file)
        network_file.close()
        
    def loadNetwork(self,filename):
        filename = filename + ".network"       
        network_file = open('Networks/'+filename, 'rb')
        self.C = pickle.load(network_file)
        network_file.close()
    
    def saveCurrentLog(self,optionalName=None,withGraphHist=True):
        log_data = dict()
        
        log_data["Links"]=self.linkLog
        log_data["S1"] = self.s1Log
        log_data["Lambda1"] = self.lambdaLog
        #log_data["Max"]= self.maxLog
        if(withGraphHist): log_data["GraphHist"] = self.graphHist
        if(self.logPopulation): log_data["PopHist"] = self.popHist
        
        if(optionalName is None): filename = time.strftime("%y_%m_%d %H_%M_%S") + ".log"
        else: filename = optionalName + ".log"
        log_file = open('Logs/'+filename, 'wb')
        pickle.dump(log_data, log_file)
        log_file.close()
        
    def loadLog(self, filename):        
        filename = filename + ".log"       
        log_file = open('Logs/'+filename, 'rb')
        log_data = pickle.load(log_file)
        log_file.close()
        
        self.linkLog = log_data["Links"]
        self.s1Log = log_data["S1"]
        self.lambdaLog = log_data["Lambda1"]
        self.maxLog = log_data["Max"]
        self.graphHist = log_data["GraphHist"]
        if(self.logPopulation): self.popHist = log_data["PopHist"]
        
    def saveCompletely(self,optionalName=None):
        if(optionalName is not None): filename = optionalName + ".jk"
        if(self.id is not None): filename = f"jkN {self.id}" + ".jk"      
        else: filename = time.strftime("%y_%m_%d %H_%M_%S") + ".jk" 
        jk_file = open('jkNetworks/'+filename, 'wb')
        pickle.dump(self, jk_file)
        jk_file.close()
        
    def loadCompletely(filename=None,network_id=None):        
        if(network_id is not None): filename = f"jkN {network_id}" + ".jk" 
        else: filename = filename + ".jk"       
        jk_file = open('jkNetworks/'+filename, 'rb')
        jkLoaded = pickle.load(jk_file)
        jk_file.close()
        return jkLoaded
        
    def saveDynamicGraph(self,optionalName=None):
        if not self.popHist:
            raise Exception("No population history exists")
        
        #prepare data for export to DynamicGraph
        nodeHist = dict()
        for i in range(self.s):
            nodeHist[i]=[-1] #all nodes start to exist before start of simulation
        
        if(optionalName is None): filename = time.strftime("%y_%m_%d %H_%M_%S") + ".gexf"
        else: filename = optionalName + ".gexf"
        
        dynGraph = DynamicGraph(nodeHist, self.graphHist, self.popHist)
        f = open("Graphs/" + filename, "w+")
        f.write(dynGraph.makeGraph())
        f.close()
        
    """Restore"""
        
    def restoreGraph(self,time):
        #import pdb; pdb.set_trace()
        edges=[(int(changed_id/self.s),changed_id % self.s) for changed_id,logs in self.graphHist.items() if len(list(filter(lambda x: x <= time, logs))) % 2 != 0]
        C=np.zeros((self.s,self.s),dtype=bool)
        for edge in edges: C[edge]=True
        return C
        
    def restorePopHist(self):
        self.popHist = dict()
        for t in range(self.getLength()):
            self.popHist[t],evalue = jkNetwork.getPFE(self.restoreGraph(t))
