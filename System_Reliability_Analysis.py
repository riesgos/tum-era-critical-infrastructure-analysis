# -*- coding: utf-8 -*-
"""

WPS that performs a reliability analysis on one critical infrastructure network with Monte Carlo Simulation.

This service uses the probability of failure at nodes for estimating the affectation to population 
(i.e. people affected by supply disruption).

Input data: 
    -Nodal exposure and damage (probability of failure at each node)
    -Network fragility (which taxonomies are source and consumer nodes?)
    -Line exposure (how are the nodes connected?)
    -Consumer Areas exposure (how many people is supplied by which lines?)

Output data:
    -Consumer Areas damage (probability of each area to have supply disruption)
    -Random sample of total affected population (raw data for performing statistics and plots)
Created on Mon Mar 11 18:12:00 2019

@author: hfrv2
"""
## Systemic Reliability Analysis of lifeline system considering cascading effects
# Amelie Hoffmann, 25.09.2018 (in Matlab only)
# Hugo Rosero, 01.04.2019 (in Python only)
# Hugo Rosero, 10.07.2019 (first version finished)
# Hugo Rosero, 30.07.2019 (modified with WPS standards)
# Hugo Rosero, 14.08.2019 (simplified for one network, compatible with GFZ web services)


import Sysrel as sr

from geoserver.wps import process

# Note: Edges and Lines are used as sysnonyms through the code
@process(
  title='System Reliability Analysis',
  description='Performs a reliability analysis on one critical infrastructure network with Monte Carlo Simulation',
  inputs={
    'DamageNodes': (dict, '''Output from Damage web-service (e.g. DEUS from GFZ) that includes the information of the exposure 
                    corresponding to the nodes of the network (node, location, taxonomy, etc.), plus the output of that web-service (e.g. probability of 
                    failure given a specific multi-hazard scenario). The dictionary structure follows the geojson standard format'''),
    'ExposureLines': (dict, '''Information of the exposure corresponding to the topology of the network (start and end node of each line, 
                    coordinates of the path for visualization, etc.). The start and end nodes coincide with node names/ids from DamageNodes.
                    The dictionary structure follows the geojson standard format'''),
    'ExposureConsumerAreas': (dict, '''Information of the exposure corresponding to the consumer areas (i.e. extent of distribution networks),
                    containing properties such as id, name, number of households. The names/ids coincide with an end node name/id in the ExposureLines dict.
                    The dictionary structure follows the geojson standard format'''),
    'NetworkFragility':(dict,'''information about the fragility parameters for the network taxonomy. The dictionary structure follows the same json standard format 
                 used for fragility of buildings''')
  },
  outputs={
    'DamageConsumerAreas': (dict, '''Damage in the consumer areas. Contains the information of the exposure, plus properties storing damage metrics, such as 
                            probability of blackout, based on Monte Carlo Simulation. The dictionary structure follows the same geojson standard format as 
                            ExposureConsumerAreas'''),
    'SampleDamageNetwork': (list, '''list of samples of global damage metric values (e.g. total affected population/households), 
                            which are result of Monte Carlo Simulation.'''),
  }
)

def run(DamageNodes=None,ExposureLines=None,ExposureConsumerAreas=None,NetworkFragility=None):
    
    if (DamageNodes==None or ExposureLines==None or ExposureConsumerAreas==None or NetworkFragility==None):
        raise ValueError("Either nodes, lines, area or fragility information are missing. Analysis aborted")

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
        nmcs=100 #number of samples
        # obtain samples of affected areas
        SampleAreas=sr.run_Monte_Carlo_simulation(Graph,source_nodes,consumer_nodes,ExposureConsumerAreas,nmcs)
        ##### ----------------------------- Post Processing ---------------------------------########
        DamageConsumerAreas,SampleDamageNetwork=sr.compute_output(SampleAreas,ExposureConsumerAreas,nmcs)
        
        return DamageConsumerAreas,SampleDamageNetwork
    
    if __name__ == '__main__':
        return main()