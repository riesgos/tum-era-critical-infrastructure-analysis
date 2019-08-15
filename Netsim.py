# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 14:49:53 2019

@author: hfrv2
"""

import numpy as np
import networkx as nx
import Constants as cons
import Sysrel as sr

def direct_hazard_action(G): 
        # Evaluate Component Fragilities
        for node in G.nodes():
            #generate random value between 0 and 1
            r_value=np.random.rand()
            #query probability of failure
            n_pof=G.nodes[node][cons.NODE_POF]
            if r_value<n_pof:
                G.nodes[node][cons.DAMAGE]=1

# UPDATE SOURCE AND CONSUMER NODES; AND EDGE WEIGHTS
# The damage increment DELTADAMAGE is reset to zero after updating the damage of the components
def update_network(G,s_nodes,c_nodes):
    # first, reduce node capacities
    for node in G.nodes():
        # if damage is close to 1 (e.g. by hazard action), set capacity to zero
        if G.node[node][cons.DAMAGE]>1.0-cons.MIN_DAMAGE:
            if node in s_nodes:
                s_nodes.remove(node)
            if node in c_nodes:
                c_nodes.remove(node)
            G.node[node][cons.CAPACITY]=0.0
        #otherwise, only reduce capacity
        else:
            G.node[node][cons.CAPACITY]-=G.node[node][cons.DELTADAMAGE]*G.node[node][cons.CAPACITY]
            #we already applied the damage increment
            G.node[node][cons.DELTADAMAGE]=0
    # now, increase edge cost and reduce their capacities
    for edge in G.edges():
        G.edges[edge][cons.CAPACITY]-=G.edges[edge][cons.DELTADAMAGE]*G.edges[edge][cons.CAPACITY]
        G.edges[edge][cons.WEIGHT]/=max(1-G.edges[edge][cons.DELTADAMAGE],cons.EPS)
        #we already applied the damage increment
        G.edges[edge][cons.DELTADAMAGE]=0


#   CASCADING EFFECTS
#   Updates component state vector with failures due to nodes disconnection
#   and overloading
def simulate_cascading_effects(G,s_nodes,c_nodes,max_iteration):
    iteration_casc=0
    component_state=0
    component_state_upd=1
    # new failures occur
    while iteration_casc<max_iteration and component_state_upd!=component_state: 
       
        component_state=component_state_upd 
        
        # DISCONNECTION FAILURE
        # assess perturbed network
        sr.evaluate_system_loads(G,s_nodes,c_nodes)
        n_loads=nx.get_node_attributes(G,cons.LOAD)
        e_loads=nx.get_edge_attributes(G,cons.LOAD)
        n_caps=nx.get_node_attributes(G,cons.CAPACITY)
        e_caps=nx.get_edge_attributes(G,cons.CAPACITY)
        for node in G.nodes():
            #if new node load exceeds its capacity
            if n_caps[node]>cons.EPS:
                n_ratio=n_loads[node]/n_caps[node]
                if n_ratio>1:
                    #reduce capacity and store de damage increment
                    n_state=1-G.nodes[node][cons.DAMAGE]
                    #update the damage level
                    G.nodes[node][cons.DAMAGE]=1-(1/n_ratio)*n_state
                    #store the damage increment
                    if n_state<cons.EPS:
                        G.nodes[node][cons.DELTADAMAGE]=n_state
                    else:
                        G.nodes[node][cons.DELTADAMAGE]=G.nodes[node][cons.DAMAGE]-(1-n_state)
        
        for edge in G.edges():
            #if new node load exceeds its capacity
            if e_caps[edge]>cons.EPS:
                e_ratio=e_loads[edge]/e_caps[edge]
                if e_ratio>1:
                    #reduce capacity and store de damage increment
                    e_state=1-G.nodes[edge][cons.DAMAGE]
                    #update the damage level
                    G.nodes[edge][cons.DAMAGE]=1-(1/e_ratio)*e_state
                    #store the damage increment
                    if e_state<cons.EPS:
                        G.nodes[edge][cons.DELTADAMAGE]=e_state
                    else:
                        G.nodes[edge][cons.DELTADAMAGE]=G.edges[edge][cons.DAMAGE]-(1-e_state)
        update_network(G,s_nodes,c_nodes)
        component_state_upd={cons.NODES:nx.get_node_attributes(G,cons.DAMAGE),cons.EDGES:nx.get_edge_attributes(G,cons.DAMAGE)}
        iteration_casc+=1
    #return component_state_upd
        
# estimate affectation to consumer areas

def set_state_consumers(ExposureConsumerAreas,Graph):
    areas_damage=[]
    for i in range(0,len(ExposureConsumerAreas[cons.FEATURES])):
        key_a_name=ExposureConsumerAreas[cons.FEATURES][i][cons.PROPERTIES][cons.NAME]
        supply_edge_dam=[Graph.edges[ed][cons.DAMAGE] for ed in Graph.edges() if Graph.edges[ed][cons.TO]==key_a_name]
        areas_damage.append(np.mean(supply_edge_dam))

# IS AN EDGE IN A PATH?
def is_edge_in_path(edge,path):
    G=nx.Graph()
    G.add_path(path)
    if edge in G.edges():
        return True
    else:
        return False
