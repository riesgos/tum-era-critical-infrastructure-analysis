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
# import matplotlib.pyplot as plt
import numpy as np

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

# SAVE ANALYSIS IN GEOJSON FILE
# creates a geojson file with the information of the input geojson files (resp. CSV), plus the attributes created during the analysis (such as fragility parameters and damage level)
def save_to_JSON(DataOutDict,filename):
    with open(filename, 'w') as outfile:
        json.dump(DataOutDict, outfile)
 

# Local execution
if __name__ == '__main__':
    ##### ----------------------------- location of input files----------------------------------------########
    #folder location
    folder_prefix = os.path.dirname(os.path.realpath(__file__))
    # Exposure data from Ecuador
    DamageNodes=import_json_to_dict(os.path.join(folder_prefix, 'E1_EPN_ExposureNodes_withDamage.geojson'))
    ExposureLines=import_json_to_dict(os.path.join(folder_prefix, 'E1_EPN_ExposureLines.geojson'))
    ExposureConsumerAreas=import_json_to_dict(os.path.join(folder_prefix, 'E1_EPN_ExposureConsumerAreas.geojson'))
    NetworkFragility=import_json_to_dict(os.path.join(folder_prefix, 'NetworkFragility.json'))
    # execute main function
    DamageConsumerAreas,SampleDamageNetwork = main()
    # save consumer areas output as geojson file
    save_to_JSON(DamageConsumerAreas,os.path.join(folder_prefix, 'E1_EPN_ExposureConsumerAreas_withDamage.geojson'))
    # make a histogram with the output vector of total affected population
    # SampleDamageNetwork_1000=[SampleDamageNetwork[i]/1000 for i in range(0,len(SampleDamageNetwork))]
    # plt.hist(SampleDamageNetwork_1000,normed=True,stacked=True)
    # plt.xlabel('Affected population / Población afectada (thousands/miles)')
    # plt.ylabel('Probability / Probabilidad')
    # plt.title('Histogram of affected population / histograma de población afectada \n Scenario/Escenario VEI>=4')
    # plt.grid(True)
    # plt.show()
    # print('mean (thousands): '+str(np.mean(SampleDamageNetwork_1000))+" , Coeff. of Variation: "+str(np.std(SampleDamageNetwork_1000)/np.mean(SampleDamageNetwork_1000)))
