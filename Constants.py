# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 11:58:25 2019

@author: hfrv2
"""

##### ----------------------------- define constants ----------------------------------------########

MIN_DAMAGE=0.1 #minimum acceptable damage. Effects from smaller damage levels are neglected
EPS=1e-100 #epsilon constant for avoiding divisions by zero
NODE_TYPE='TYPE'
LINE_TYPE='TYPE'
DAMAGE='DAM'
ISSLAVE='SLV'
DELTADAMAGE='DDAM'
DAMAGESTATS1='MDAM'#damage statistics
DAMAGESTATS2='SDAM'#damage statistics
WEIGHT='WEIGHT'
VOLTAGE='VOLTAGE'
EDGES='edges'#keyword for edges in dictionary
NODES='nodes'#keyword for nodes in dictionary
NODE_TYPES=['GENERATION','SUBSTATION','CONSUMER','PATHPOINT']#node types: 1- source, 2- transmissison, 3- consumer, 4- path point (for edges)

NODE_NAME="Name" #keyword for node name. It might change depending of the keywords of input files
LINE_NAME='Name' #keyword for line name (optional). It might change depending of the keywords of input files
LOAD="LOAD"
CAPACITY="CAP"
NODE_POF="ProbFailure"
AREA_POF="ProbDisruption"
AREA_HOUSEHOLDS="NumberHouseholds"
DATA="data"
META="meta"
TAXONOMY="taxonomy"
TAXONOMIES="taxonomies"
SOURCE="source"
CONSUMER="consumer"
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

