# -*- coding: utf-8 -*-
"""
Python module for network simulation. Contains the main steps of the Monte Carlo Simulation:
    - Direct Hazard Action: applies a random damage state to the system, given probabilities of failures assigned to the nodes by the damage web service
    - Cascading Effects: simulates systemic failures due to overloading. Affects nodes and lines
    - State of consumer areas: estimates the affectation to the consumer areas, based on the damage level of the supplier lines
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

def direct_hazard_action(G): 
        # Evaluate Component Fragilities
        for node in G.nodes():
            n_pof=G.nodes[node][cons.NODE_POF]
            # if component is certainly damaged
            if n_pof>1-cons.EPS:
                G.nodes[node][cons.NODE_DAMAGE]=1
            # if component is certainly working
            elif n_pof<cons.EPS:
                G.nodes[node][cons.NODE_DAMAGE]=0
            else:
                #generate random value between 0 and 1
                r_value=np.random.rand()
                #query probability of failure
                if r_value<n_pof:
                    G.nodes[node][cons.NODE_DAMAGE]=1

# UPDATE SOURCE AND CONSUMER NODES; AND EDGE WEIGHTS
# The damage increment DELTADAMAGE is reset to zero after updating the damage of the components
def update_network(G,s_nodes,c_nodes):
    # first, reduce node capacities
    for node in G.nodes():
        # if damage is larger than the critical damage (e.g. by hazard action), set capacity to zero
        if G.node[node][cons.NODE_DAMAGE]>cons.CRIT_DAMAGE:
            if node in s_nodes:
                s_nodes.remove(node)
            if node in c_nodes:
                c_nodes.remove(node)
            G.node[node][cons.CAPACITY]=0.0
            # damaged node also makes its edges unavailable
            for edge in G.edges(node):
               G.edges[edge][cons.WEIGHT]=1/cons.EPS 
                
        #otherwise, only reduce capacity
        else:
            G.node[node][cons.CAPACITY]-=G.node[node][cons.NODE_DELTADAMAGE]*G.node[node][cons.CAPACITY]
            #we already applied the damage increment
            G.node[node][cons.NODE_DELTADAMAGE]=0.0
    # now, increase edge cost and reduce their capacities
    for edge in G.edges():
        # if damage is larger than the critical damage, set capacity to zero and weight arbitrarily large
        if G.edges[edge][cons.LINE_DAMAGE]>cons.CRIT_DAMAGE:
            G.edges[edge][cons.CAPACITY]=0.0
            G.edges[edge][cons.WEIGHT]=1/cons.EPS
        else:
            G.edges[edge][cons.CAPACITY]-=G.edges[edge][cons.LINE_DELTADAMAGE]*G.edges[edge][cons.CAPACITY]
            G.edges[edge][cons.WEIGHT]/=max(1-G.edges[edge][cons.LINE_DELTADAMAGE],cons.EPS)
        #we already applied the damage increment
        G.edges[edge][cons.LINE_DELTADAMAGE]=0.0


def evaluate_system_loads(G,s_nodes,c_nodes):

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
def simulate_cascading_effects(G,s_nodes,c_nodes,max_iteration):
    iteration_casc=0
    component_state=0
    component_state_upd=1
    # new failures occur
    while iteration_casc<max_iteration and component_state_upd!=component_state: 
       
        component_state=component_state_upd 
        
        # DISCONNECTION FAILURE
        # assess perturbed network
        evaluate_system_loads(G,s_nodes,c_nodes)
        n_loads=nx.get_node_attributes(G,cons.LOAD)
        e_loads=nx.get_edge_attributes(G,cons.LOAD)
        n_caps=nx.get_node_attributes(G,cons.CAPACITY)
        e_caps=nx.get_edge_attributes(G,cons.CAPACITY)
        for node in G.nodes():
            #if new node load exceeds its capacity
            if n_caps[node]>cons.EPS:
                n_ratio=n_loads[node]/n_caps[node]
                if n_ratio>1:
                    #print('Casc effects nodes')
                    #reduce capacity and store de damage increment
                    n_state=1-G.nodes[node][cons.NODE_DAMAGE]
                    #update the damage level
                    G.nodes[node][cons.NODE_DAMAGE]=1-(1/n_ratio)*n_state
                    #store the damage increment
                    if n_state<cons.MIN_DAMAGE:
                        G.nodes[node][cons.NODE_DELTADAMAGE]=n_state
                    else:
                        G.nodes[node][cons.NODE_DELTADAMAGE]=G.nodes[node][cons.NODE_DAMAGE]-(1-n_state)
        
        for edge in G.edges():
            #if new node load exceeds its capacity
            if e_caps[edge]>cons.EPS:
                e_ratio=e_loads[edge]/e_caps[edge]
                if e_ratio>1:
                    #print('Casc effects edges')
                    #reduce capacity and store de damage increment
                    e_state=1-G.edges[edge][cons.LINE_DAMAGE]
                    #update the damage level
                    G.edges[edge][cons.LINE_DAMAGE]=1-(1/e_ratio)*e_state
                    #store the damage increment
                    if e_state<cons.EPS:
                        G.edges[edge][cons.LINE_DELTADAMAGE]=e_state
                    else:
                        G.edges[edge][cons.LINE_DELTADAMAGE]=G.edges[edge][cons.LINE_DAMAGE]-(1-e_state)
        update_network(G,s_nodes,c_nodes)
        component_state_upd={cons.NODES:nx.get_node_attributes(G,cons.NODE_DAMAGE),cons.EDGES:nx.get_edge_attributes(G,cons.LINE_DAMAGE)}
        iteration_casc+=1
    #return component_state_upd
        
# estimate affectation to consumer areas

def set_state_consumers(ExposureConsumerAreas,Graph):
    areas_damage=[]
    for i in range(0,len(ExposureConsumerAreas[cons.FEATURES])):
        key_a_name=ExposureConsumerAreas[cons.FEATURES][i][cons.PROPERTIES][cons.AREA_NAME]
        if len(Graph.edges(key_a_name))>0:
            supply_edge_dam=[Graph.edges[edge][cons.LINE_DAMAGE] for edge in Graph.edges(key_a_name)]
            if np.mean(supply_edge_dam)>cons.CRIT_DAMAGE:
                areas_damage.append(1)# 1 means has a blackout
            else:
                areas_damage.append(0)# 0 means no blackout
        else:#it is an isolated consumer area (no service)
            areas_damage.append(1)
    return areas_damage

# IS AN EDGE IN A PATH?
def is_edge_in_path(edge,path):
    G=nx.Graph()
    G.add_path(path)
    if edge in G.edges():
        return True
    else:
        return False
    
def OD_node_betweenness_centrality(G, s_nodes, c_nodes, k=None, normalized=True, weight=None,
                                seed=None):
    r"""Compute betweenness centrality for edges, considering paths between source and consumer nodes.
    
    This function is a modified version of networkx package function edge_betweenness_centrality()

    Betweenness centrality of an edge $e$ is the sum of the
    fraction of all-pairs shortest paths that pass through $e$

    .. math::

       c_B(e) =\sum_{s\in SV,t \in CV} \frac{\sigma(s, t|e)}{\sigma(s, t)}

    where $SV$ is the set of source nodes, $CV$ is the set of consumer nodes,
    $\sigma(s, t)$ is the number of shortest $(s, t)$-paths, and $\sigma(s, t|e)$ 
    is the number of     those paths passing through edge $e$ [2]_.

    Parameters
    ----------
    G : graph
      A NetworkX graph.

    k : int, optional (default=None)
      If k is not None use k node samples to estimate betweenness.
      The value of k <= n where n is the number of nodes in the graph.
      Higher values give better approximation.
    
    s_nodes: set of source nodes
    
    c_nodes: set of consumer nodes

    normalized : bool, optional
      If True the betweenness values are normalized by $2/(n(n-1))$
      for graphs, and $1/(n(n-1))$ for directed graphs where $n$
      is the number of nodes in G.

    weight : None or string, optional (default=None)
      If None, all edge weights are considered equal.
      Otherwise holds the name of the edge attribute used as weight.

    seed : integer, random_state, or None (default)
        Indicator of random number generation state.
        See :ref:`Randomness<randomness>`.
        Note that this is only used if k is not None.

    Returns
    -------
    edges : dictionary
       Dictionary of edges with betweenness centrality as the value.


    For weighted graphs the edge weights must be greater than zero.
    Zero edge weights can produce an infinite number of equal length
    paths between pairs of nodes.

    """
    betweenness = dict.fromkeys(G.nodes, 1.0)  # b[v]=0 for v in G
    # b[e]=0 for e in G.edges()
    if seed is None:
        random.seed(seed)
    betweenness.update(dict.fromkeys(G.nodes(), 1.0))
    if k is None:
        sample_s_nodes = s_nodes
        sample_c_nodes = c_nodes
    else:
        sample_s_nodes = random.sample(s_nodes, min(k,len(s_nodes)))
        sample_c_nodes = random.sample(c_nodes, min(k,len(c_nodes)))
    for s in sample_s_nodes:
        for c in sample_c_nodes:
            if c!=s:
                #shortest paths
                if weight is None:  # use BFS
                    sp = nx.shortest_path(G, source=s,target=c)
                else:  # use Dijkstra's algorithm
                    sp = nx.shortest_path(G,source=s,target=c,weight=weight)
                    # accumulation
                for i in range(0,len(sp)):
                    betweenness[sp[i]] +=1
    return betweenness
    
    
def OD_edge_betweenness_centrality(G, s_nodes, c_nodes, k=None, normalized=True, weight=None,
                                seed=None):
    r"""Compute betweenness centrality for edges, considering paths between source and consumer nodes.
    
    This function is a modified version of networkx package function edge_betweenness_centrality()

    Betweenness centrality of an edge $e$ is the sum of the
    fraction of all-pairs shortest paths that pass through $e$

    .. math::

       c_B(e) =\sum_{s\in SV,t \in CV} \frac{\sigma(s, t|e)}{\sigma(s, t)}

    where $SV$ is the set of source nodes, $CV$ is the set of consumer nodes,
    $\sigma(s, t)$ is the number of shortest $(s, t)$-paths, and $\sigma(s, t|e)$ 
    is the number of     those paths passing through edge $e$ [2]_.

    Parameters
    ----------
    G : graph
      A NetworkX graph.

    k : int, optional (default=None)
      If k is not None use k node samples to estimate betweenness.
      The value of k <= n where n is the number of nodes in the graph.
      Higher values give better approximation.
    
    s_nodes: set of source nodes
    
    c_nodes: set of consumer nodes

    normalized : bool, optional
      If True the betweenness values are normalized by $2/(n(n-1))$
      for graphs, and $1/(n(n-1))$ for directed graphs where $n$
      is the number of nodes in G.

    weight : None or string, optional (default=None)
      If None, all edge weights are considered equal.
      Otherwise holds the name of the edge attribute used as weight.

    seed : integer, random_state, or None (default)
        Indicator of random number generation state.
        See :ref:`Randomness<randomness>`.
        Note that this is only used if k is not None.

    Returns
    -------
    edges : dictionary
       Dictionary of edges with betweenness centrality as the value.


    For weighted graphs the edge weights must be greater than zero.
    Zero edge weights can produce an infinite number of equal length
    paths between pairs of nodes.

    """
    betweenness = dict.fromkeys(G.edges, 1.0)  # b[v]=0 for v in G
    # b[e]=0 for e in G.edges()
    if seed is None:
        random.seed(seed)
    betweenness.update(dict.fromkeys(G.edges(), 1.0))
    if k is None:
        sample_s_nodes = s_nodes
        sample_c_nodes = c_nodes
    else:
        sample_s_nodes = random.sample(s_nodes, min(k,len(s_nodes)))
        sample_c_nodes = random.sample(c_nodes, min(k,len(c_nodes)))
    for s in sample_s_nodes:
        for c in sample_c_nodes:
            if c!=s:
                #shortest paths
                if weight is None:  # use BFS
                    sp = nx.shortest_path(G, source=s,target=c)
                else:  # use Dijkstra's algorithm
                    sp = nx.shortest_path(G,source=s,target=c,weight=weight)
                # accumulation
                for i in range(0,len(sp)-1):
                    try:
                        betweenness[(sp[i],sp[i+1])] +=1
                    except:
                        betweenness[(sp[i+1],sp[i])] +=1
    return betweenness
