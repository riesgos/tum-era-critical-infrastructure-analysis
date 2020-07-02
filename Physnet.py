# -*- coding: utf-8 -*-
"""
Created on Thu Jul  2 13:54:57 2020

@author: hfrv2
"""

# -*- coding: utf-8 -*-
"""
Python module for network simulation with physics based models for computing loads.
    - other functions:
        - Network update: removes failed/isolated source and consumer nodes, and updates surviving component capacities, 
        for recomputing the loads during the simulation of cascading effects
        - Is edge in path?: returns true or false, depending whether an edge is in a given path.

Created on Wed Aug 14 14:49:53 2019

@author: hfrv2
"""

import numpy as np
import networkx as nx
import Constants as cons
import random



def evaluate_system_loads(G):

    loads_initial={}
    edge_loads_initial={}   
    # Determine shortest path between all source and target nodes and set
    # the nr of shortest paths passing through each node in network as
    # initialize loads with 1 for assigning a minimum non zero capacity to each component
    #loads_initial={no:1 for no in G.nodes()}
    #edge_loads_initial={ed:1 for ed in G.edges()}
    # evaluate inital origin-destination betweenness centrality
    # determine initial number of connecting paths from all source node to
    # each supply node
    #print('\t start load evaluation')
    #path_dict=nx.shortest_path(G,weight=cons.WEIGHT)
    #i_s=0
    param_k=min(30,len(G.nodes()))
    loads_initial=OD_node_betweenness_centrality(G,s_nodes,c_nodes,weight=cons.WEIGHT,normalized=False,k=param_k)
    #edge_loads_initial=nx.edge_betweenness_centrality(G,weight=cons.WEIGHT,normalized=False,k=param_k)
    edge_loads_initial=OD_edge_betweenness_centrality(G,s_nodes,c_nodes,weight=cons.WEIGHT,normalized=False,k=param_k)
    #for s_node in s_nodes:
        #print('\t source node Nr. '+str(i_s))
        #i_s=i_s+1
        #for c_node in c_nodes: 
            #try:
                #a_path=path_dict[s_node][c_node] 
            #except:
                #continue
            #for node in G.nodes():
                #if node in a_path:
                    #loads_initial[node]+=1
            #for edge in G.edges():
                #if is_edge_in_path(edge,a_path):
                    #edge_loads_initial[edge]+=1
    
    # calculate initial capacity of each node using alpha factors
    node_attribs={no:{cons.LOAD:loads_initial[no]} for no in G.nodes()}
    edge_attribs={ed:{cons.LOAD:edge_loads_initial[ed]} for ed in G.edges()}
    nx.set_node_attributes(G,node_attribs)
    nx.set_edge_attributes(G,edge_attribs)

#   CASCADING EFFECTS
#   Updates component state vector with failures due to nodes disconnection
#   and overloading

