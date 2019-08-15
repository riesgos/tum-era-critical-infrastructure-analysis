# -*- coding: utf-8 -*-
"""
Python module containing constants that are called in the functions of Netsim and Sysrel modules.

Created on Wed Aug 14 11:58:25 2019

@author: hfrv2
"""

##### ----------------------------- define constants ----------------------------------------########

# useful constants for the functions
MIN_DAMAGE=0.1 #minimum acceptable damage. Effects from smaller damage levels are neglected
EPS=1e-100 #epsilon constant for avoiding divisions by zero
EDGES='edges'#keyword for edges in dictionary
NODES='nodes'#keyword for nodes in dictionary

# general geojson file
TYPE='type'
NAME='name'
CRS='crs'
FEATURES="features"
PROPERTIES="properties"
GEOMETRY="geometry"
COORDINATES="coordinates"

# load and capacity
LOAD="LOAD"#this property is added in code
CAPACITY="CAP"#this property is added in code

#Exposure nodes geojson file
NODE_TYPE="taxonomy"
NODE_NAME="Name_Node" #keyword for node name
NODE_POF="ProbFailure"#this property should come in the output from damage web service
NODE_DAMAGE='DAM'#this property is added in code 
NODE_DELTADAMAGE='DDAM'#this property is added in code

#Exposure Lines geojson file
LINE_NAME='Name' #keyword for line name (optional)
FROM="FROM"
TO="TO"
LENGTH="length"
LINE_TYPE="taxonomy"
VOLTAGE='voltage'
WEIGHT='WEIGHT'#this property is added in code
LINE_DAMAGE='DAM'#this property is added in code
LINE_DELTADAMAGE='DDAM'#this property is added in code

#Exposure areas geojson file
AREA_NAME='Name' #keyword for area name (must coincide with name of a consumer node)
AREA_POF="Prob_Disruption" #this property is added in code
AREA_POPULATION="population"

# Network fragility json file
DATA="data"
META="meta"
TAXONOMY="taxonomy"
TAXONOMIES="taxonomies"
SOURCE="source"
CONSUMER="consumer"