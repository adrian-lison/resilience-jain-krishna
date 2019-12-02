# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 15:50:03 2019

@author: Sascha
"""

import numpy as np
import scipy as sp
import scipy.linalg as la
import scipy.stats as st
import networkx as nx

class jkTopology:

    def __init__(self,jkn):
        
        self.C=jkn.C
        self.pfe = jkn.PFE
        self.lambda1 = jkn.lambda1
        self.s = jkn.s
        self.m=jkn.m
        
        #print("Identify Core")
        G=nx.DiGraph(np.transpose(self.C))
        
        self.core_nodes, self.core_bs_nodes, self.periphery_nodes = jkTopology.findCoreAndPheriphery(
                G, self.C, self.pfe, self.lambda1) #bs_ basic subgraphs 

        self.core=G.subgraph(self.core_nodes)
        #import pdb; pdb.set_trace()
        self.coreC=self.C[np.ix_(list(self.core_nodes),list(self.core_nodes))]

        #print("Calculate entire network statistics")
        #measures for the entire network
        self.nw_edges,self.nw_comple,self.nw_between_central_mean, self.nw_between_central_max, self.nw_between_central_skew=self.topo_analysis(G,self.C)
        
        #print("Calculate overall core statistics")        
        #measures for the core:
        self.core_edges, self.core_comple, self.core_between_central_mean, self.core_between_central_max, self.core_between_central_skew=self.topo_analysis(self.core, self.coreC) 
        
        #print("Calculate average basic subgraph of core statistics")
        #measures for each of the subgraphs in the core
        
        self.core_path_mult=list() #connectivity of each of the subgraphs in the core
        self.core_path_mult_complete=list()
        self.cyclelen=list()
        for core_bs_nodes_i in self.core_bs_nodes:
            #self.core_bs[i]=G.subgraph(core_bs_nodes[i])
            core_bs_mat = self.C[np.ix_(list(core_bs_nodes_i),list(core_bs_nodes_i))]
            core_bs_graph = nx.DiGraph(np.transpose(core_bs_mat))
            core_bs_graph = nx.relabel_nodes(core_bs_graph,dict(zip(core_bs_graph,core_bs_nodes_i)))
            cycles = list(nx.simple_cycles(core_bs_graph))
            
            import pdb; pdb.set_trace()
            self.cyclelen.append(len(cycles))
            self.core_path_mult.append(jkTopology.path_mult(core_bs_mat))
            self.core_path_mult_complete.append(jkTopology.path_mult(core_bs_mat,divide_complete=True))  
            
        self.av_core_path_mult=np.average(self.core_path_mult)
        self.av_core_path_mult_complete=np.average(self.core_path_mult_complete)
        self.av_core_cyclelen=np.average(self.cyclelen)
        
    def findCoreAndPheriphery(G, matrix, pfe, lambda1):
        #import pdb; pdb.set_trace()
        pfe_i=np.where(pfe!=0)[0]
        pfegraph=G.subgraph(pfe_i)     
        sccs = [scc for scc in nx.strongly_connected_components(pfegraph) if len(scc)>1]
        sccs_core = list()
        for scc in sccs:
            scc_mat=matrix[np.ix_(list(scc),list(scc))]
            if np.round(jkTopology.getPFE(scc_mat)[1],5)==np.round(lambda1,5):
                sccs_core.append(scc)
        corenodes=set(set().union(*sccs_core))
        pheripherynodes=set(pfegraph).difference(corenodes)
        #import pdb; pdb.set_trace()
        return corenodes,sccs_core,pheripherynodes

    def topo_analysis(self,G, C): #G: nettworkx graph, C: adjecency matrix
        """creates topological measurs for networks given as networkx"""       
        num_edges = np.sum(C)
        clausen_comple=self.complexity(C) #graph complexity according to clausen
        between_central_nodes=nx.betweenness_centrality(G)
        #import pdb; pdb.set_trace()
        between_central_mean = sum(between_central_nodes.values())/len(between_central_nodes)
        between_central_max = max(between_central_nodes.values())
        between_central_skew = st.skew([i for i in between_central_nodes.values()])
        
        return num_edges, clausen_comple, between_central_mean, between_central_max, between_central_skew
                   
    def complexity(self,C):
        """return offdiagonal complexity [Claussen 2007]"""
        s=C.shape[0]
        c = np.zeros((2*s, 2*s))
        degrees = np.count_nonzero(C,axis=0)+np.count_nonzero(C,axis=1) #array of degrees
        #import pdb; pdb.set_trace()
        for i in range(s):
            for j in range(s):
                if(C[i,j]==1 and degrees[i]<=degrees[j]):
                    #import pdb; pdb.set_trace()
                    c[int(degrees[i]), int(degrees[j])] += 1
        atilde = np.zeros(s)
        for i in range(s):
            atilde[i] = np.sum(np.diag(c, i))
        a = atilde / np.sum(atilde)
        compl = st.entropy(a)
        return compl
    
    def path_mult(C,divide_complete=False): #pathway mulztiplicity within core subgraph
        numC=C.astype(np.int64)
        s=np.shape(C)[0] #number of nodes
        A=np.zeros((s,s))
        g=0
        tempA=C.copy()
        fac=0
        #import pdb; pdb.set_trace()
        for n in range(2,s+1):
            # A = lA + C^n
            tempA = np.dot(tempA,numC)
            A += tempA
            if(divide_complete):
                fac = (s-1)**(n-1) - fac #sum([x for x in map(lambda i: (-1)**(i+1)*(s-1)**(n-i),np.array(range(1,n)))])
                g += fac
            #print(f"A:{A[0,0]}, g:{g}")
         
        if(not divide_complete): g=1   
        return np.trace(A)/(g*s)
      
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
        raise Exception("No eigenvector without negative entries.") #no such eigen vector found
           
        
        