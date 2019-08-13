# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 18:12:00 2019

@author: hfrv2
"""
## Systemic Reliability Analysis considering of interdependent lifeline systems considering cascading effects
# Amelie Hoffmann, 25.09.2018 (in Matlab only)
# Hugo Rosero, 01.04.2019 (in Python only)
# Hugo Rosero, 10.07.2019 (first version finished)
# Hugo Rosero, 30.07.2019 (modified with WPS standards)


import numpy as np
import scipy as sp
import networkx as nx
from scipy.stats import lognorm
import csv
import os
import copy
import json
import gdal
from geoserver.wps import process

# Note: Edges and Lines are used as sysnonyms through the code
@process(
  title='System Reliability Analysis',
  description='Performs a reliability analysis on critical infrastructure networks',
  inputs={
    'networks': (list, 'list of names of the infrastructure networks to be considered'),
    'NodeFilenames': (list, 'list of locations of the geojson files with the nodes of the networks. The order coincides with the order in network list'),
    'LineFilenames': (list, 'list of locations of the geojson files with the lines of the networks. The order coincides with the order in network list'),
    'areas': (list, 'list of locations of the geojson files with the consumer areas of the networks. The order coincides with the order in network list'),
    'Hazards':(list,'list of names of the hazard intensities that affect the network'),
    'Intensities':(list,'list of locations of the tiff images of the hazard intensities that affect the network. The order coincides with the order in hazards list'),
    'Fragilities':(dict,'dict of dicts of locations of the csv files containing the parameters of fragility functions. Keys correspond to network list, and each item is a dict with keys from hazard list'),
    'InterFilenames':(dict,'dict of dicts of locations of the csv files containing the interdependencies among networks. Keys correspond to network list, and each item is a dict with keys from network list, except the key of the item (i.e. no diagonal-like keys in this data structure)'),
    'NodeProps': (list, 'list of locations of the csv files containing additional node properties. The order coincides with the order in network list'),
    'LineProps': (list, 'list of locations of the csv files containing additional line properties. The order coincides with the order in network list')
  },
  outputs={
    'nodes_out': (list, 'list of locations of the geojson files with the nodes of the networks, modified with new output properties'),
    'lines_out': (list, 'list of locations of the geojson files with the lines of the networks, modified with new output properties'),
    'areas_out': (list, 'list of locations of the geojson files with the distribution areas of the networks, modified with new output properties')
  }
)

def run(networks=None,NodeFilenames=None,LineFilenames=None, InterFilenames=None,Fragfilenames=None,
        hazards=None,TIFFfilenames=None,Outputfilenames=None, NodeProps=None,LineProps=None):


    ##### ----------------------------- define constants ----------------------------------------########

    MIN_DAMAGE=0.1 #minimum acceptable damage. Effects from smaller damage levels are neglected
    EPS=1e-100 #epsilon constant for avoiding divisions by zero
    NODE_TYPE='TYPE'
    DAMAGE='DAM'
    ISSLAVE='SLV'
    DELTADAMAGE='DDAM'
    DAMAGESTATS1='MDAM'#damage statistics
    DAMAGESTATS2='SDAM'#damage statistics
    WEIGHT='WEIGHT'
    EDGES='edges'#keyword for edges in dictionary
    NODES='nodes'#keyword for nodes in dictionary
    NODE_TYPES=['GENERATION','SUBSTATION','CONSUMER','PATHPOINT']#node types: 1- source, 2- transmissison, 3- consumer, 4- path point (for edges)

    NODE_NAME="Name" #keyword for node name. It might change depending of the keywords of input files
    LINE_NAME='Name' #keyword for line name (optional). It might change depending of the keywords of input files
    LOAD="LOAD"
    CAPACITY="CAP"

    #common keywords in GEOJSON files:
    TYPE='type'
    NAME='name'
    CRS='crs'
    FEATURES="features"
    PROPERTIES="properties"
    GEOMETRY="geometry"
    COORDINATES="coordinates"
    FROM="ORIGIN"
    TO="DESTINY"
    LENGTH="LENGTH"
    
    if hazards==None:
        hazards = ['LH'] # (default) list of considered hazards (if one hazard has more than one intensity, then use indexed labels (e.g. EQ1, EQ2 instead EQ for earthquake)
    IM_ATTR=['IM_'+hz for hz in hazards]# default list of intensity measures
    FR_PAR=['FR_'+hz for hz in hazards]# default list of fragility parameters
    FR_BOUNDS=['FRB_'+hz for hz in hazards]# default list of fragility bounds
    
    if networks==None:
        #networks = ['POWER','WATER']# List of considered networks
        networks = ['POWER']# List of considered networks (default)

    ##### ----------------------------- location of input and output files (defaults)----------------------------------------########
    #folder location
    folder_prefix = os.path.dirname(os.path.realpath(__file__))
    if NodeFilenames==None or LineFilenames==None:
        #filenames in csv (if needed)
        #filenames={n:{NODES:folder_prefix+'\\nodelist_{}.csv'.format(networks.index(n)+1),EDGES:folder_prefix+'\\edgelist_{}.csv'.format(networks.index(n)+1)} for n in networks}
        filenames={n:{NODES:folder_prefix+'\\E1_Nodes_{}.geojson'.format(networks.index(n)+1),EDGES:folder_prefix+'\\E1_Lines_{}.geojson'.format(networks.index(n)+1)} for n in networks}
    else:
        filenames={networks[i]:{NODES:NodeFilenames[i],EDGES:LineFilenames[i]} for i in range(0,len(networks))}
    if InterFilenames==None:
        #filenames for interdependence measures (when 2 or more netowrks considered)
        InterFilenames={m:{n:folder_prefix+'\\interface_{}_{}.csv'.format(networks.index(m)+1,networks.index(n)+1) for n in networks if n!=m} for m in networks}

    if Fragfilenames==None:
        #filenames for fragility function parameters
        Fragfilenames={m:{hz:folder_prefix+'\\E1_component_fragilities_{}_{}.csv'.format(networks.index(m)+1,hazards.index(hz)+1) for hz in hazards} for m in networks}
    if TIFFfilenames==None:
        #filenames for tiff output files from hazard models (intensities)
        TIFFfilenames={h:folder_prefix+'\\E1_Lahar_HMax_S4.tif' for h in hazards}
    else:
        TIFFfilenames={hazards[i]:TIFFfilenames[i] for i in range(0,len(hazards))}
    #filenames for outputs
    NodeOutput=[folder_prefix+'\\E1_NodesOut_{}.geojson'.format(networks.index(n)+1) for n in networks]
    LineOutput=[folder_prefix+'\\E1_LinesOut_{}.geojson'.format(networks.index(n)+1) for n in networks]
    AreaOutput=[folder_prefix+'\\E1_AreasOut_{}.geojson'.format(networks.index(n)+1) for n in networks]
    #Outputfilenames={n:{NODES:folder_prefix+'\\E1_NodesOut_{}.geojson'.format(networks.index(n)+1),EDGES:folder_prefix+'\\E1_LinesOut_{}.geojson'.format(networks.index(n)+1)} for n in networks}


    ##### ----------------------------- simulation settings ----------------------------------------########

    nr_int = 2*sp.misc.comb(len(networks),2,exact=True) # max number of interdependence matrices
    alpha = {net:0.5 for net in networks}    # tolerance parameter between 0 and 1 for each network
    mcs=20      # number of iterations in the monte carlo simulation
    max_iteration=10   # max nr of iterations if no covergence is reached
    seed=677


    ##### ----------------------------- initialization of variables ---------------------------------########
 


    s_nodes={n:None for n in networks}# source nodes
    t_nodes={n:None for n in networks}# transmission nodes
    c_nodes={n:None for n in networks}# consumer nodes
    paths_upd={n:None for n in networks}# path list
    connectivity={n:np.zeros(mcs) for n in networks}    # total connectivity loss in system
    damage_propagation={n:[[] for x in range(0,mcs)] for n in networks}    # component damage propagation in each cycle



    ##### ----------------------------- Auxiliary functions called in the code ---------------------------------########

    # IMPORT TIFF IMAGES AS ARRAYS, INCLUDING COORD TRANSFORM PARAMS
    def import_tiff_to_array(filename):
        ds = gdal.Open(filename)
        # GetGeoTransform gives the (x, y) origin of the top left pixel,
        # the x and y resolution of the pixels, and the rotation of the raster.
        # The information is returned as tuple:
        # (TL x, X resolution, X rotation, TL y, Y rotation, y resolution)
        TiffTransform = ds.GetGeoTransform()
        # Read the raster as an array
        TiffArray = ds.ReadAsArray()
        return TiffTransform,TiffArray


    # IMPORT JSON FILES AND CREATE DICTIONARY
    def import_json_to_dict(filename):
        data={}
        with open(filename, encoding="utf8") as f:
            s=f.read()
            s = s.replace('ï»¿','')#for some reason, this weird char sequence might appear when opening
            data=json.loads(s)
        return data
        
    # IMPORT CSV TO MATRIX
    def import_csv_to_matrix(filename):
        RawList=[]
        with open(filename, 'r') as csvFile:
            reader = csv.reader(csvFile)
            for row in reader:
                RawList.append(row)
        csvFile.close()
        header=list(RawList[0])
        C=RawList[1:len(RawList)]
        for i in range(0,len(C)):
            for j in range(0,len(C[i])):
                #format numeric data as int or float instead of string
                try:
                    C[i][j]=int(C[i][j])
                except:
                    try:
                        C[i][j]=float(C[i][j])
                    except:
                        continue
        return header, C

    # DESCRIBE NETWORK BY READING CSV FILES
    def network_description(n,filenames):
        C={}
        header={}
        header[EDGES],C[EDGES]=import_csv_to_matrix(filenames[n][EDGES])
        header[NODES],C[NODES]=import_csv_to_matrix(filenames[n][NODES])
        return header,C

    # DEFINE NODE TO NODE INTERDEPENDENCIES
    def interdependence_description(m,n,filenames,Gs):
        #check if networks have interdependencies
        try:
            header,Inter_list= import_csv_to_matrix(filenames[m][n])
            #return if there is no interdependencies
        except:
            return
        n_inter=len(Inter_list[:][0])
        for i in range(0,n_inter):
            mr_node=Inter_list[i][0]
            sl_node=Inter_list[i][1]
            Gs[m].nodes[sl_node][ISSLAVE][n][mr_node]=float(Inter_list[i][2])

    # DEFINE FRAGILITY PARAMETERS
    def fragility_params(filename,G,f_par):
        try:
            Fheader,Frag_pars= import_csv_to_matrix(filename)
            #return if there are no fragilities
        except:
            #print('No fragilities found')
            return
        for node in G.nodes():
        
            G.nodes[node][f_par]=[
                    [Frag_pars[i][2],Frag_pars[i][3]]
                    for i in range(0,len(Frag_pars)) 
                    if Frag_pars[i][0]==1+NODE_TYPES.index(G.nodes[node][NODE_TYPE])]
        # for edges, we assume that the path points are of the same type
        # e.g. a power transmission line has the same type of towers
        # or a pipeline segment has the same diameter and materials
        for edge in G.edges():
            #store the fragility parameters for each coordinate point in the edge
            G.edges[edge][f_par]=[[
                    [Frag_pars[i][2],Frag_pars[i][3]]
                    for i in range(0,len(Frag_pars)) 
                    if Frag_pars[i][0]==NODE_TYPES[3]
                    ]for poin in G.edges[edge][COORDINATES]]

    # IS AN EDGE IN A PATH?
    def is_edge_in_path(edge,path):
        G=nx.Graph()
        G.add_path(path)
        if edge in G.edges():
           return True
        else:
           return False

    # Initialize the graphs with data from CSV files
    def load_from_CSV(networks, Gs, hazards, filenames):
        headers={}
        Gs={net:nx.Graph()for net in networks}
        edgelist={n:None for n in networks}
        nodelist={n:None for n in networks}
        edgelength={n:None for n in networks}
        nodelength={n:None for n in networks}
   
        for n in networks:
            headers[n],C=network_description(n,filenames)
            # determine edgelist and nodelist of network based on  data given in
            # network_description
            edgelist[n]=C[EDGES]
            nodelist[n]=C[NODES]
            edgelength[n]=len(edgelist[n])
            nodelength[n]=len(nodelist[n])
            #Define minimum attributes for nodes and edges
            NodeBunch={nodelist[n][nod][0]:{
                NODE_TYPE:C[NODES][nod][1],
                DAMAGE:0,
                COORDINATES:(0,0),
                ISSLAVE:{n:{} for n in networks},
                DELTADAMAGE:0,
                DAMAGESTATS1:0.0,
                DAMAGESTATS2:0.0
                } for nod in range(0,nodelength[n])}
            EdgeBunch=[(edgelist[n][x][0],edgelist[n][x][1],{
                WEIGHT:float(edgelist[n][x][headers[n][EDGES].index(WEIGHT)]),
                DAMAGE:0,
                COORDINATES:[],
                DELTADAMAGE:0,
                DAMAGESTATS1:0.0,
                DAMAGESTATS2:0.0
                }) for x in range(0,edgelength[n])]
            Gs[n].add_edges_from(EdgeBunch)
            del EdgeBunch
        
            nx.set_node_attributes(Gs[n],NodeBunch)
        
            #Add user-defined attributes (input file shall include header)
            edge_attributes={(edgelist[n][x][0],edgelist[n][x][1]):{
                     headers[n][EDGES][y]:C[EDGES][x][y] for y in range(0,len(edgelist[n][0]))
                     } for x in range(0,edgelength[n])}
            nx.set_edge_attributes(Gs[n],edge_attributes)
    
            node_attributes={nodelist[n][x][0]:{
                headers[n][NODES][y]:C[NODES][x][y] for y in range(0,len(nodelist[n][0]))
                } for x in range(0,nodelength[n])}
            nx.set_node_attributes(Gs[n],node_attributes)
            del C
        return Gs, headers

    # Initialize the graphs with data from GEOJSON files
    def load_from_JSON(networks,hazards, filenames,weight=LENGTH):
        Gs={net:nx.Graph()for net in networks}
        headers={n:{NODES:None,EDGES:None} for n in networks}
        for n in networks:
            dataEdges=import_json_to_dict(filenames[n][EDGES])
            dataNodes=import_json_to_dict(filenames[n][NODES])
            headers[n]={NODES:{k:dataNodes[k] for k in dataNodes.keys() if k!=FEATURES},EDGES:{k:dataEdges[k] for k in dataEdges.keys() if k!=FEATURES}}
            EdgesFea=dataEdges[FEATURES]
            NodesFea=dataNodes[FEATURES]        
            #Define minimum attributes for nodes and edges
            NodeBunch={NodesFea[i][PROPERTIES][NODE_NAME]:{
                NODE_TYPE:NodesFea[i][PROPERTIES][NODE_TYPE],
                DAMAGE:0,
                COORDINATES:NodesFea[i][GEOMETRY][COORDINATES],
                ISSLAVE:{n:{} for n in networks},
                DELTADAMAGE:0,
                DAMAGESTATS1:0.0,
                DAMAGESTATS2:0.0
                } for i in range(0,len(NodesFea))}
            EdgeBunch=[(EdgesFea[i][PROPERTIES][FROM],EdgesFea[i][PROPERTIES][TO],{
                WEIGHT:float(EdgesFea[i][PROPERTIES][weight]),
                DAMAGE:0,
                COORDINATES:EdgesFea[i][GEOMETRY][COORDINATES][0],
                DELTADAMAGE:0,
                DAMAGESTATS1:0.0,
                DAMAGESTATS2:0.0
                }) for i in range(0,len(EdgesFea))]
            Gs[n].add_edges_from(EdgeBunch)
            del EdgeBunch
        
            nx.set_node_attributes(Gs[n],NodeBunch)
        
            #create extra node attributes (only from 'properties' key)
            for attr in NodesFea[0][PROPERTIES].keys():
                if attr==NODE_NAME:
                    continue
                NodeAttr={NodesFea[i][PROPERTIES][NODE_NAME]:{attr:NodesFea[i][PROPERTIES][attr]} for i in range(0,len(NodesFea))}
                nx.set_node_attributes(Gs[n],NodeAttr)
            #create extra edge attributes (only from 'properties' key)
            for attr in EdgesFea[0][PROPERTIES].keys():
                EdgeAttr={(EdgesFea[i][PROPERTIES][FROM],EdgesFea[i][PROPERTIES][TO]):{attr:EdgesFea[i][PROPERTIES][attr]} for i in range(0,len(EdgesFea))}
                nx.set_node_attributes(Gs[n],EdgeAttr)
        return Gs, headers

    # Create the graph with interdependencies. By default, the source data are geojson files
    def load_network_data(networks, hazards, filenames, InterFilenames,IType="geojson"):
        node_types={n:None for n in networks}
    
        if IType!="geojson":
            Graphs,headers=load_from_CSV(networks, hazards, filenames)
        else:
            Graphs,headers=load_from_JSON(networks, hazards, filenames)
    
        #Add default intensity measures
        for n in networks:
            node_IMs={nod:{
                hz:None for hz in IM_ATTR
                } for nod in Graphs[n].nodes()}
            nx.set_node_attributes(Graphs[n],node_IMs)
    
            edge_IMs={ed:{
                hz:{} for hz in IM_ATTR
                } for ed in Graphs[n].edges()}
            nx.set_edge_attributes(Graphs[n],edge_IMs)
        
            #Add default fragility parameters
            node_FRs={nod:{
                fr:[] for fr in FR_PAR
                } for nod in Graphs[n].nodes()}
            nx.set_node_attributes(Graphs[n],node_FRs)
    
            edge_FRs={ed:{
                fr:[] for fr in FR_PAR
                } for ed in Graphs[n].edges()}
            nx.set_edge_attributes(Graphs[n],edge_FRs)
        
            #Add default fragility bounds
            node_FRBs={nod:{
                fr:[] for fr in FR_BOUNDS
                } for nod in Graphs[n].nodes()}
            nx.set_node_attributes(Graphs[n],node_FRBs)
    
            edge_FRBs={ed:{
                fr:[] for fr in FR_BOUNDS
                } for ed in Graphs[n].edges()}
            nx.set_edge_attributes(Graphs[n],edge_FRBs)
        
            # create vector of all pairs of source and consumption nodes in network
            Types=nx.get_node_attributes(Graphs[n],NODE_TYPE)
            node_types[n]={typ:[nod for nod in Graphs[n].nodes() if Types[nod]==typ] for typ in NODE_TYPES}
            del Types

        # Assign fragilitiy functions
        for m in networks:
            for hz in hazards:
                fragility_params(Fragfilenames[m][hz],Graphs[m],FR_PAR[hazards.index(hz)])

        # Load interdependence structure
        for m in networks:
            for n in networks:
                if n!=m:
                    interdependence_description(m,n,InterFilenames,Graphs) 
        return Graphs,node_types,headers

    # evaluation of system loads
    # load is computed as the number of shortest paths that pass through the component (node or edge)
    def evaluate_system_loads(G,s_nodes,c_nodes):
        loads_initial={}
        edge_loads_initial={}   
        # Determine shortest path between all source and target nodes and set
        # the nr of shortest paths passing through each node in network as
        # initialize loads with 1 for assigning a minimum non zero capacity to each component
        loads_initial={no:1 for no in G.nodes()}
        edge_loads_initial={ed:1 for ed in G.edges()}
        # evaluate inital origin-destination betweenness centrality
        # determine initial number of connecting paths from all source node to
        # each supply node
        for s_node in s_nodes:
            for c_node in c_nodes:
                path_list=nx.all_shortest_paths(G,s_node,c_node,weight=WEIGHT)
                path_list=list(path_list)               
                if len(path_list)>0:
                    for a_path in path_list:
                        for node in G.nodes():
                            if node in a_path:
                                loads_initial[node]+=1
                        for edge in G.edges():
                            if is_edge_in_path(edge,a_path):
                                edge_loads_initial[edge]+=1
    
        # calculate initial capacity of each node using alpha factors
        node_attribs={no:{LOAD:loads_initial[no]} for no in G.nodes()}
        edge_attribs={ed:{LOAD:edge_loads_initial[ed]} for ed in G.edges()}
        nx.set_node_attributes(G,node_attribs)
        nx.set_edge_attributes(G,edge_attribs)

    # ASSIGN HAZARD INTENSITIES
    # Reads the tiff file and assigns to the points of interest the intensity value
    # if the point is outside the extent of the tiff file, the intensity is zero
    def set_Hazard_Intensities(hazard,G): 
        h=hazards.index(hazard)
        HTransf, HArray= import_tiff_to_array(TIFFfilenames[hazard])
        for node in G.nodes():
            NodeLoc=G.nodes[node][COORDINATES]
            x_index = int((NodeLoc[0] - HTransf[0]) / HTransf[1])
            y_index = int((NodeLoc[1] - HTransf[3]) / HTransf[5])
            #Assign non-zero IM to locations inside the extent of the hazard array
            try:
                G.nodes[node][IM_ATTR[h]]=HArray[y_index, x_index]
                #otherwise, set intensity to zero
            except:
                G.nodes[node][IM_ATTR[h]]=0.0
        for edge in G.edges():
            EdgeLocs=G.edges[edge][COORDINATES]
            G.edges[edge][IM_ATTR[h]]=[]
            for poin in EdgeLocs:
                x_index = int((poin[0] - HTransf[0]) / HTransf[1])
                y_index = int((poin[1] - HTransf[3]) / HTransf[5])
                #Assign non-zero IM to locations inside the extent of the hazard array
                try:
                    G.edges[edge][IM_ATTR[h]].append(HArray[y_index, x_index])
                    #otherwise, set intensity to zero
                except:
                    G.edges[edge][IM_ATTR[h]].append(0.0)

    # SET FRAGILITY BOUNDS
    # Computes the bounds between discrete damage states given an intensity measure. 
    # All the information is read and modified in the graph G   
    def set_Fragility_Bounds(hazard,G):
        h=hazards.index(hazard)
        #case for nodes
        for node in G.nodes():
            # get the intensity measures for the node, and the corresponding 
            # parameters of fragility functions
            im_h=G.nodes[node][IM_ATTR[h]]
            Fparams=G.nodes[node][FR_PAR[h]]
            #number of fragility functions for this intensity measure
            fr_states=len(Fparams)
            if fr_states>0:
                #damage states
                d_states=range(0,len(Fparams))
                #means
                lloc=[np.log(Fparams[ds][0]) for ds in d_states]
                #standard deviations
                lstd=[Fparams[ds][1] for ds in d_states]
                # compute the bounds and store in graph
                G.nodes[node][FR_BOUNDS[h]]=[
                        lognorm.cdf(im_h,s=lstd[ds],loc=lloc[ds]) 
                        for ds in d_states
                        ]
            # if the component has no fragility for this intensity measure
            else:
                continue
        #case for edges
        for edge in G.edges():
            # get the intensity measures for the node, and the corresponding 
            # parameters of fragility functions
            im_h=G.edges[edge][IM_ATTR[h]]
            PFparams=G.edges[edge][FR_PAR[h]]
            #number of path points
            n_pts=len(im_h)
            # if there are points
            if n_pts>0:
                for poin in range(0,n_pts):
                    #damage state parameters
                    Fparams=PFparams[poin]
                    if len(Fparams)>0:
                        d_states=range(0,len(Fparams))
                        #means
                        lloc=[np.log(Fparams[ds][0]) for ds in d_states]
                        #standard deviations
                        lstd=[Fparams[ds][1] for ds in d_states]
                        # compute the bounds and store in graph
                        G.edges[edge][FR_BOUNDS[h]][poin]=[
                            lognorm.cdf(im_h,s=lstd[ds],loc=lloc[ds]) 
                            for ds in d_states
                            ]
            else:
                continue

    #MODIFY FRAGILITY PARAMETERS OF A NODE
    # For dynamic exposure, the parameters of fragility functions are modified, based on the present damage level of the node
    def  modify_node_fragilities(G,node,damage,h,hazards):
        i_myhazard=hazards.index(h)
        #Here we assume that the hazards happen sequentially
        for hz in range(i_myhazard+1,len(hazards)):
            Frag_pars=G.nodes[node][FR_PAR[hz]]
            G.nodes[node][FR_PAR[hz]]=[
                    [Frag_pars[i][0]*(1-damage),Frag_pars[i][1]*(1-damage)]
                    for i in range(0,len(Frag_pars))]
        
    #MODIFY FRAGILITY PARAMETERS OF AN EDGE
    # For dynamic exposure, the parameters of fragility functions are modified, based on the present damage level of the edges
    def  modify_edge_fragilities(G,edge,damages,h,hazards):
        i_myhazard=hazards.index(h)
        #Here we assume that the hazards happen sequentially
        for hz in range(i_myhazard+1,len(hazards)):
            Frag_pars=G.edges[edge][FR_PAR[hz]]
            p_count=range(0,len(G.edges[edge][COORDINATES]))
            G.edges[edge][FR_PAR[hz]]=[
                [
                        [Frag_pars[i][poin][0]*(1-damages[poin]),
                         Frag_pars[i][poin][1]*(1-damages[poin])]
                        for i in range(0,len(Frag_pars[poin]))]
                for poin in p_count]
                    
                
    #   DIRECT HAZARD ACTION OPERATION
    #   Simulates direct perturbation action on network components by all
    #   hazards


    #   ATTENTION: For more than two hazards and networks additional component fragility
    #   functions must be added. 

    #   input: hazard list, network information as Graph object
    #   the following attributes of Graph object are modified: damage state and fragility parameters of nodes and edges
    ##
    def direct_hazard_action(hazards,G): 
        # Evaluate Component Fragilities
        for h in hazards:
            i_h=hazards.index(h)
            #refresh the fragility functions for the known IM
            set_Fragility_Bounds(h,G)
            for node in G.nodes():
                #generate random value between 0 and 1
                r_value=np.random.rand()
                #set the bounds of fragility functions conditioned to the IM
                n_bounds=G.nodes[node][FR_BOUNDS[i_h]]
                ln_bounds=len(n_bounds)
                #read the current state of the node
                n_state=1-G.node[node][DAMAGE]
                #save initial state
                n_state0=n_state
                # if there is damage due to the hazard
                if (ln_bounds>0) and (n_state>0.0) and (n_bounds[0]>r_value):
                    ds=0
                    while ds<ln_bounds and n_bounds[ds]>r_value:
                        ds+=1
                    #reduce the serviceability of the node
                    n_state*=(1-ds/ln_bounds)
                
                    # components with damage level larger than 90% are considered failed
                    if n_state<0.1:
                        n_state=0.0
                    else:
                        #modify node fragilities for the upcoming hazards, given the damage increase
                        modify_node_fragilities(G,node,1-n_state,h,hazards)
                    #save the new damage level of the node
                    G.node[node][DAMAGE]=1.0-n_state
                    #save the damage relative change
                    G.node[node][DELTADAMAGE]=1-n_state/n_state0
            
            #the same for edges    
            for edge in G.edges():

                #set the bounds of fragility functions conditioned to the IM
                n_bounds=G.edges[edge][FR_BOUNDS[i_h]]
                ledge_pts=len(n_bounds)
                ln_bounds=[len(n_bounds[pt]) for pt in range(0,ledge_pts)]
                #generate random vector between 0 and 1
                r_values=[sp.rand() for x in range(0,ledge_pts)]
                #current serviceability of the edge
                n_state=1-G.edges[edge][DAMAGE]
                #store initial serviceability
                n_state0=n_state
                # if there is damage due to the hazard
                if (ledge_pts>0) and (n_state>0.0) and any([n_bounds[x][0]<r_values[x] for x in range(0,ledge_pts)]):
                    ds=[0 for x in range(0,ledge_pts)]
                    #increase damage level to each point location
                    for poin in range(0,ledge_pts):
                        while n_bounds[poin][ds[poin]]<r_values[poin] and ds[poin]<ln_bounds:
                            ds[poin]+=1
                    #reduce the serviceability of the line
                    n_state*=(1-max(ds)/ln_bounds)
                    #modify edge fragilities for the upcoming hazards, given the damage increase
                    modify_edge_fragilities(G,edge,[ds[x]/ln_bounds for x in range(0,ledge_pts)],h,hazards)
                    #save the new damage level of the line
                    # components with damage level larger than 90% are considered failed
                    if n_state<0.1:
                        n_state=0.0
                    G.edges[edge][DAMAGE]=1.0-n_state
                    G.edges[edge][DELTADAMAGE]=1-n_state/n_state0

    # UPDATE SOURCE AND CONSUMER NODES; AND EDGE WEIGHTS
    # The damage increment DELTADAMAGE is reset to zero after updating the damage of the components
    def update_network(G,s_nodes,c_nodes):
        # first, reduce node capacities
        for node in G.nodes():         
            G.node[node][CAPACITY]-=G.node[node][DELTADAMAGE]*G.node[node][CAPACITY]
            #we already applied the damage increment
            G.node[node][DELTADAMAGE]=0
            # damaged source nodes do not produce resource anymore
            if node in s_nodes and G.node[node][DAMAGE]>1.0-MIN_DAMAGE:
                s_nodes.remove(node)
                # damaged consumer nodes cannot receive resources anymore
            if node in c_nodes and G.node[node][DAMAGE]>1.0-MIN_DAMAGE:
                c_nodes.remove(node)
        
        # now, increase edge cost and reduce their capacities
        for edge in G.edges():
            G.edges[edge][CAPACITY]-=G.edges[edge][DELTADAMAGE]*G.edges[edge][CAPACITY]
            G.edges[edge][WEIGHT]/=max(1-G.edges[edge][DELTADAMAGE],EPS)
            #we already applied the damage increment
            G.edges[edge][DELTADAMAGE]=0

    #   INTERNAL DAMAGE PROPAGATION
    #   Updates component state vector with failures due to nodes disconnection
    #   and overloading

    def internal_damage_propagation(component_state_upd,G,s_nodes,c_nodes,max_iteration):
        iteration_idp=0
        component_state=0
        # new failures occur
        while iteration_idp<max_iteration and component_state_upd!=component_state: 
       
            component_state=component_state_upd 
        
            # DISCONNECTION FAILURE
            # assess perturbed network
            evaluate_system_loads(G,s_nodes,c_nodes)
            n_loads=nx.get_node_attributes(G,LOAD)
            e_loads=nx.get_edge_attributes(G,LOAD)
            n_caps=nx.get_node_attributes(G,CAPACITY)
            e_caps=nx.get_edge_attributes(G,CAPACITY)
            for node in G.nodes():
                #if new node load exceeds its capacity
                if n_caps[node]>EPS:
                    n_ratio=n_loads[node]/n_caps[node]
                    if n_ratio>1:
                        #reduce capacity and store de damage increment
                        n_state=1-G.nodes[node][DAMAGE]
                        #update the damage level
                        G.nodes[node][DAMAGE]=1-(1/n_ratio)*n_state
                        #store the damage increment
                        if n_state<EPS:
                            G.nodes[node][DELTADAMAGE]=n_state
                        else:
                            G.nodes[node][DELTADAMAGE]=G.nodes[node][DAMAGE]-(1-n_state)
        
            for edge in G.edges():
                #if new node load exceeds its capacity
                if e_caps[edge]>EPS:
                    e_ratio=e_loads[edge]/e_caps[edge]
                    if e_ratio>1:
                        #reduce capacity and store de damage increment
                        e_state=1-G.nodes[edge][DAMAGE]
                        #update the damage level
                        G.nodes[edge][DAMAGE]=1-(1/e_ratio)*e_state
                        #store the damage increment
                        if e_state<EPS:
                            G.nodes[edge][DELTADAMAGE]=e_state
                        else:
                            G.nodes[edge][DELTADAMAGE]=G.edges[edge][DAMAGE]-(1-e_state)
            update_network(G,s_nodes,c_nodes)
            component_state_upd={NODES:nx.get_node_attributes(G,DAMAGE),EDGES:nx.get_edge_attributes(G,DAMAGE)}
            iteration_idp+=1


    #   NETWORK INTERDEPENDENT PROPAGATION OF DAMAGES
    #   updates component state vectors of both networks with failures due to
    #   failed components in the master network
    #   this function still requires revision

    def external_damage_propagation(G,network,component_state_upd,s_nodes,c_nodes):
        #component_state_upd=component_state    
        for n in networks:
            for node in G[n].nodes():
                # if node is not completely damaged
                if G[n].nodes[node][DAMAGE]<1.0-MIN_DAMAGE:
                    # get the list of master nodes for this node
                    mstr_nodes=G[n].nodes[node][ISSLAVE]
                #run through the other networks
                networks2=copy.deepcopy(networks)
                networks2.remove(n)
                for m in networks2:
                    # check if there are master nodes
                    try:
                        m_mstr_nodes=mstr_nodes[m]
                        for m_node in m_mstr_nodes:
                            #get damage state of master node
                            m_damage=G[m].node[m_node][DAMAGE]
                            # if damage level is significant and interdependence randomly affects slave node
                            if m_damage>MIN_DAMAGE and sp.rand()<m_mstr_nodes[m_node]:
                                #we assume that the damage increase in slave node is proportional to the damage level of master node
                                G[n].node[node][DELTADAMAGE]=m_damage
                    # ignore if there are no master nodes for this network
                    except:
                        continue
        update_network(G[n],s_nodes[n],c_nodes[n])
        #obtain the updated states of the system
        component_state_upd[n]={NODES:nx.get_node_attributes(G[n],DAMAGE),EDGES:nx.get_edge_attributes(G[n],DAMAGE)}

    # SAVE ANALYSIS IN GEOJSON FILE
    # creates a geojson file with the information of the input geojson files (resp. CSV), plus the attributes created during the analysis (such as fragility parameters and damage level)
    def save_graph_to_JSON(G,headers,NodeOut, LineOut):
        NodeOutDict={k: headers[NODES][k] for k in headers[NODES].keys()}
        NodeOutDict[FEATURES]=[]
        for node in G.nodes():
            TheNode=G.nodes[node]
            NodeFea={'type':'Feature',GEOMETRY:{},PROPERTIES:{}}
            NodeFea[GEOMETRY]={'type':'Point',COORDINATES:TheNode[COORDINATES]}
            NodeFea[PROPERTIES]={k:TheNode[k] for k in TheNode.keys() if k!=COORDINATES}
            NodeFea[PROPERTIES][NAME]=node
            NodeOutDict[FEATURES].append(NodeFea)
        with open(NodeOut, 'w') as outfile:
            json.dump(NodeOutDict, outfile)
        EdgeOutDict={k: headers[EDGES][k] for k in headers[EDGES].keys()}
        EdgeOutDict[FEATURES]=[]
        for edge in G.edges():
            TheEdge=G.edges[edge]
            EdgeFea={'type':'Feature'}
            EdgeFea[GEOMETRY]={'type':'MultiLineString',COORDINATES:[TheEdge[COORDINATES]]}
            EdgeFea[PROPERTIES]={k:TheEdge[k] for k in TheEdge.keys() if k!=COORDINATES}
            EdgeOutDict[FEATURES].append(EdgeFea)   
        with open(LineOut, 'w') as outfile2:
            json.dump(EdgeOutDict, outfile2)   



    ## Notes
    # Matlab version: considering node failures only
    # Python version: extended to line failures

    ########## -------------------------------------------------------------------- #################
    ########## --------------------------MAIN SCRIPT STARTS HERE------------------- #################        
    ########## -------------------------------------------------------------------- #################

    ##### ----------------------------- Load network data ---------------------------------########

    Graphs,node_types,headers=load_network_data(networks,hazards,filenames,InterFilenames)
    s_nodes0={net:node_types[net][NODE_TYPES[0]] for net in networks}
    c_nodes0={net:node_types[net][NODE_TYPES[2]] for net in networks}

    ##### --------------------- Assess unperturbed system and capacities -------------------########


    for n in networks:
        evaluate_system_loads(Graphs[n],s_nodes0[n],c_nodes0[n])
        node_caps={no:{CAPACITY:(1+alpha[n])*Graphs[n].node[no][LOAD]} for no in Graphs[n].nodes()}
        edge_caps={ed:{CAPACITY:(1+alpha[n])*Graphs[n].edges[ed][LOAD]} for ed in Graphs[n].edges()}
        nx.set_node_attributes(Graphs[n],node_caps)
        nx.set_edge_attributes(Graphs[n],edge_caps)

    ##### --------------------- Assign hazard intensities to the networks -------------------########    
    
    for n in networks:
        for h in hazards:
            set_Hazard_Intensities(h,Graphs[n])
        
    ##### ----------------------------- Monte Carlo Simulation ----------------------------########
  
 
    component_state_dha={net:[] for net in networks }
    component_state_idp={net:[] for net in networks }
    component_state_edp={net:[] for net in networks }
    n_loads0=nx.get_node_attributes(Graphs[n],LOAD)
    e_loads0=nx.get_edge_attributes(Graphs[n],LOAD)

    #store damage samples in list
    G_damages=[]
    #initialize graph with damage as only attribute
    DamGraph={net:nx.Graph()for net in networks}

    for n in networks:
        DamGraph[n].add_edges_from(Graphs[n].edges(),DAMAGE=0)
        nx.set_node_attributes(DamGraph[n],0,DAMAGE)
    

 
    for i in range(0,mcs):
        myflag=False
        i_component_state_dha={}
        i_component_state_idp={}
        i_component_state_edp={}
        iteration_edp=[0,0]
        #modify the networks in a copy of them
        NetG=copy.deepcopy(Graphs)
        s_nodes=copy.deepcopy(s_nodes0)
        c_nodes=copy.deepcopy(c_nodes0)
        nr_its=0
        G_damages.append(copy.deepcopy(DamGraph))
        ## Direct Hazard Action
        for n in networks:
            # Simulate effects of hazard action on components
            direct_hazard_action(hazards,NetG[n])
        
            # store coponent damage states
            i_component_state_dha={NODES:nx.get_node_attributes(NetG[n],DAMAGE),EDGES:nx.get_edge_attributes(NetG[n],DAMAGE)}
            component_state_dha[n].append(i_component_state_dha)
        
            # Update network description (weights and capacities)
            update_network(NetG[n],s_nodes[n],c_nodes[n])    
       

        # at least one component has significant damage (>10% damage)
        if any(
            [max(component_state_dha[n][i][NODES].values())>MIN_DAMAGE for n in networks]
            ) or any(
            [max(component_state_dha[n][i][EDGES].values())>MIN_DAMAGE for n in networks]
            ):
            myflag=True
            ## Damage Propagation
            component_state_edp=component_state_dha
            # at least one new component has failed
            while component_state_idp!=component_state_edp:
            
                component_state_idp=component_state_edp
            
                ## Internal Damage Propagation
            
                for n in networks:
                    i_component_state_idp={NODES:nx.get_node_attributes(NetG[n],DAMAGE),EDGES:nx.get_edge_attributes(NetG[n],DAMAGE)}
                    # if any component in this network has failed
                    if max(i_component_state_idp[NODES].values())>MIN_DAMAGE or max(i_component_state_idp[EDGES].values())>MIN_DAMAGE:   
                    
                        # if there are surviving supply and consumtion nodes
                        if len(s_nodes[n])>0 and len(c_nodes[n])>0:
                        
                            # start idp cycle
                            internal_damage_propagation(i_component_state_idp,NetG[n],s_nodes[n],c_nodes[n],max_iteration)
                        #else:
                            # loss of functionality, i.e. no more flow in network
                        update_network(NetG[n],s_nodes[n],c_nodes[n])                         
                  
                    # no disconnection or cascading failures in this network              
            
                    # Save Results
                    component_state_idp[n][i]=i_component_state_idp
                     

                ## External Damage Propagation
                i_component_state_edp=i_component_state_idp
                external_damage_propagation(NetG,networks,i_component_state_edp,s_nodes,c_nodes)
                nr_its+=1
            
                # if max iterations reached
                if nr_its==max_iteration: 
                    print('max. number of iterations reached in EDP \n')
                    #force the results to be equal
                    component_state_edp[n]=component_state_idp[n]

            # Save Results

            ## Evaluate System Performance
            for n in networks:
                # there are surviving supply and consumption nodes
                if len(s_nodes[n])>0 and len(c_nodes[n])>0:
                    # determine shortest path between surviving source and consumption nodes
                    evaluate_system_loads(NetG[n],s_nodes[n],c_nodes[n])
                    # calculate connectivity measure
                    n_loads=nx.get_node_attributes(NetG[n],LOAD)
                    e_loads=nx.get_edge_attributes(NetG[n],LOAD)

                    connectivity[n][i]= 1/len(c_nodes0[n])*sum([n_loads[nod]/n_loads0[nod] for nod in c_nodes0[n]])
                else:# no sources, no consumers. No flow
                    connectivity[n][i]=0.0
 
        else:
            # No component failures in any of the networks
            for n in networks:
                connectivity[n][i]=1
            # total affected people: none, as no failed components
        #save output data of this sample
        for n in networks:
            nx.set_node_attributes(G_damages[i][n],nx.get_node_attributes(NetG[n],DAMAGE),DAMAGE)
            nx.set_edge_attributes(G_damages[i][n],nx.get_edge_attributes(NetG[n],DAMAGE),DAMAGE)
    
        print('Iteration: #{}'.format(i))

    ##### ----------------------------- Post Processing ----------------------------########

    for n in networks:
        for node in Graphs[n].nodes():
            SampleDamages=[]
            for i in range(0,mcs):
                SampleDamages.append(G_damages[i][n].node[node][DAMAGE])
            Graphs[n].node[node][DAMAGESTATS1]=np.mean(SampleDamages)
            Graphs[n].node[node][DAMAGESTATS2]=np.std(SampleDamages)
        for edge in Graphs[n].edges():
            SampleDamages=[]
            for i in range(0,mcs):
                SampleDamages.append(G_damages[i][n].edges[edge][DAMAGE])
            Graphs[n].edges[edge][DAMAGESTATS1]=np.mean(SampleDamages)
            Graphs[n].edges[edge][DAMAGESTATS2]=np.std(SampleDamages)        
        save_graph_to_JSON(Graphs[n],headers[n],NodeOutput[n],LineOutput[n])
    print('finished')    
    return NodeOutput,LineOutput, AreaOutput

    
    
    

