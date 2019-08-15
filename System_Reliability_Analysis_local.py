# -*- coding: utf-8 -*-
"""

Local execution test file of the content of the WPS that performs a reliability analysis on one critical infrastructure network with Monte Carlo Simulation.

This service uses the probability of failure at nodes for estimating the affectation to population 
(i.e. people affected by supply disruption).

Input data: 
    -Nodal exposure and damage (probability of failure at each node)
    -Network fragility (which taxonomies are source and consumer nodes?)
    -Line exposure (how are the nodes connected?)
    -Consumer Areas exposure (how many people is supplied by which lines?)
The input data is processed with a function that creates dicts from json files

Output data:
    -Consumer Areas damage (probability of each area to have supply disruption)
    -Random sample of total affected population (raw data for performing statistics and plots)
Created on Mon Mar 11 18:12:00 2019
Created on Wed Aug 14 11:15:58 2019

@author: hfrv2
"""

import Sysrel as sr 
import json
import os

########## -------------------------------------------------------------------- #################
########## --------------------------MAIN FUNCTION----------------------------- #################        
########## -------------------------------------------------------------------- #################
def main():
    ##### ----------------------------- Load network data -------------------------------########
    Graph,source_nodes,consumer_nodes=sr.load_network_data(DamageNodes,ExposureLines,NetworkFragility)
    ##### --------------------- Assess unperturbed system and capacities ----------------########
    sr.evaluate_system_loads(Graph,source_nodes,consumer_nodes)
    alpha=1.5#safety factor (>=1.0) for estimating capacity based on initial loads 
    sr.assign_initial_capacities(Graph,alpha)
    ##### ----------------------------- Monte Carlo Simulation --------------------------########
    nmcs=20 #number of samples
    # obtain samples of affected areas
    SampleAreas=sr.run_Monte_Carlo_simulation(Graph,source_nodes,consumer_nodes,ExposureConsumerAreas,nmcs)
    ##### ----------------------------- Post Processing ---------------------------------########
    DamageConsumerAreas,SampleDamageNetwork=sr.compute_output(SampleAreas,ExposureConsumerAreas,nmcs)
        
    return DamageConsumerAreas,SampleDamageNetwork
    
# IMPORT JSON FILES AND CREATE DICTIONARY
def import_json_to_dict(filename):
    data={}
    with open(filename, encoding="utf8") as f:
        s=f.read()
        s = s.replace('ï»¿','')#for some reason, this weird char sequence might appear when opening
        data=json.loads(s)
    return data

# Local execution
if __name__ == '__main__':
    ##### ----------------------------- location of input files----------------------------------------########
    #folder location
    folder_prefix = os.path.dirname(os.path.realpath(__file__))
    # Exposure data from Ecuador
    DamageNodes=import_json_to_dict(folder_prefix+'\\E1_EPN_ExposureNodes.geojson')
    ExposureLines=import_json_to_dict(folder_prefix+'\\E1_EPN_ExposureLines.geojson')
    ExposureConsumerAreas=import_json_to_dict(folder_prefix+'\\E1_EPN_ExposureConsumerAreas.geojson')
    NetworkFragility=import_json_to_dict(folder_prefix+'\\NetworkFragility.json')
    DamageConsumerAreas,SampleDamageNetwork = main()