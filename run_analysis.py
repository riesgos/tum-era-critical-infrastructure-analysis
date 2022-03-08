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

import argparse
import json
import os

import shapely.geometry
import numpy as np
import matplotlib.pyplot as plt

import shakemap
import fragility
import Sysrel as sr

def run_network_simulation(DamageNodes, ExposureLines, NetworkFragility, ExposureConsumerAreas):
    ##### ----------------------------- Load network data -------------------------------########
    Graph,source_nodes,consumer_nodes=sr.load_network_data(DamageNodes,ExposureLines,NetworkFragility)
    ##### --------------------- Assess unperturbed system and capacities ----------------########
    sr.evaluate_system_loads(Graph,source_nodes,consumer_nodes)
    alpha=1.5#safety factor (>=1.0) for estimating capacity based on initial loads 
    sr.assign_initial_capacities(Graph,alpha)
    ##### ----------------------------- Monte Carlo Simulation --------------------------########
    nmcs=50 #number of samples
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

def make_histogram(SampleDamageNetwork):
    # make a histogram with the output vector of total affected population
    SampleDamageNetwork_1000=[SampleDamageNetwork[i]/1000 for i in range(0,len(SampleDamageNetwork))]
    plt.hist(SampleDamageNetwork_1000,normed=True,stacked=True)
    plt.xlabel('Affected population / Población afectada (thousands/miles)')
    plt.ylabel('Probability / Probabilidad')
    plt.title('Histogram of affected population / histograma de población afectada')
    plt.grid(True)
    plt.show()
    print('mean (thousands): '+str(np.mean(SampleDamageNetwork_1000))+" , Coeff. of Variation: "+str(np.std(SampleDamageNetwork_1000)/np.mean(SampleDamageNetwork_1000)))

def evaluate_ProbFailure_oneIntensity(fragility_file,im_file,Nodes):

    intensity_provider = shakemap.Shakemaps.from_file(im_file).to_intensity_provider()
    fragility_provider = fragility.Fragility.from_file(fragility_file).to_fragility_provider()

    # add the fragility value for the element in the field "ProbFailure"

    for feature in Nodes['features']:
        centroid = shapely.geometry.shape(feature['geometry']).centroid
        lon, lat = centroid.x, centroid.y

        intensity, units = intensity_provider.get_nearest(lon=lon, lat=lat)
        # there is just one damage state for each taxonomy
        damage_state = fragility_provider.get_damage_states_for_taxonomy(feature['properties']['taxonomy'])[0]
        p = damage_state.get_probability_for_intensity(intensity, units)
        feature['properties']['ProbFailure'] = p

def evaluate_ProbFailure_multiIntensity(fragility_file_list,im_file_list,Nodes):
    
    #first intensity (lahars, height) for creating the field ProbFailure in Nodes
    evaluate_ProbFailure_oneIntensity(fragility_file_list[0],im_file_list[0],Nodes)
    # we assume that the fragility and intensity lists have the same order ( for lahars, second intensity and fragility is velocity)
    for i in range(1,len(fragility_file_list)):
        intensity_provider = shakemap.Shakemaps.from_file(im_file_list[i]).to_intensity_provider()
        fragility_provider = fragility.Fragility.from_file(fragility_file_list[i]).to_fragility_provider()
        #print("Matched frag and haz")
        # add the fragility value for the element in the field "ProbFailure"
        for feature in Nodes['features']:
            centroid = shapely.geometry.shape(feature['geometry']).centroid
            lon, lat = centroid.x, centroid.y
            intensity, units = intensity_provider.get_nearest(lon=lon, lat=lat)
            # there is just one damage state for each taxonomy
            damage_state = fragility_provider.get_damage_states_for_taxonomy(feature['properties']['taxonomy'])[0]
            p = damage_state.get_probability_for_intensity(intensity, units)
            q = feature['properties']['ProbFailure']
            #print(feature['properties']['Name Node']+": ("+str(lon)+", "+str(lat)+"), "+" "+str(intensity)+" "+str(q)+" "+str(p))
            feature['properties']['ProbFailure'] =max(p,q)

    
def main():
    ##### ----------------------------- location of input files----------------------------------------########

    argparser = argparse.ArgumentParser(
        description='Script to compute the probability of disruption given a shakemap')
    argparser.add_argument(
        '--intensity_file', nargs='*',
        help='File with the hazard intensities, for example a shakemap')
    argparser.add_argument(
        '--country',
        help='Country for which the simulation should be done. Supported: chile, ecuador')
    argparser.add_argument(
        '--hazard',
        help='Hazard for chosing the fragility functions. Supported: earthquake, lahar')
    argparser.add_argument(
        '--output_file',
        help='Name of the output file for the consumer areas with damage.')

    args = argparser.parse_args()

    prefixes_by_hazard = {
        'earthquake': 'EQ',
        'lahar': ['LH_maxheight','LH_maxvelocity']
    }
    
    prefixes_by_country = {
        'chile': 'C1',
        'ecuador': 'E1',       
        'peru':'P1',
    }

    if args.country not in prefixes_by_country.keys():
        raise Exception('{0} is not a supported country'.format(args.country))

    country_prefix = prefixes_by_country[args.country]

    if args.hazard not in prefixes_by_hazard.keys():
        raise Exception('{0} is not a supported hazard'.format(args.hazard))

    fragility_file_prefix = prefixes_by_hazard[args.hazard]
    im_file_list=args.intensity_file
        
    #folder location
    folder_prefix = os.path.dirname(os.path.realpath(__file__))
    # Exposure data 
    DamageNodes=import_json_to_dict(os.path.join(folder_prefix, country_prefix + '_EPN_ExposureNodes.geojson'))
    ExposureLines=import_json_to_dict(os.path.join(folder_prefix, country_prefix + '_EPN_ExposureLines.geojson'))
    ExposureConsumerAreas=import_json_to_dict(os.path.join(folder_prefix, country_prefix + '_EPN_ExposureConsumerAreas.geojson'))

    # if the hazard uses more than one intensity measure
    if args.hazard in ['lahar']:
        fragility_files = [os.path.join(folder_prefix, ffp + '_NetworkFragility.json') for ffp in fragility_file_prefix]
        evaluate_ProbFailure_multiIntensity(fragility_files,im_file_list,DamageNodes)
        NetworkFragility=import_json_to_dict(fragility_files[0])# sources and terminals do not depend on the intensity measure
    else: #hazard with one single intensity measure
        
        if isinstance(im_file_list,list):
            im_file_list=im_file_list[0]
        if isinstance(fragility_file_prefix,list):
            fragility_file_prefix=fragility_file_prefix[0]   
        fragility_file = os.path.join(folder_prefix, fragility_file_prefix + '_NetworkFragility.json')
        evaluate_ProbFailure_oneIntensity(fragility_file,im_file_list,DamageNodes)
        NetworkFragility=import_json_to_dict(fragility_file)
    
    
    

    # execute main function
    DamageConsumerAreas,SampleDamageNetwork = run_network_simulation(DamageNodes, ExposureLines, NetworkFragility, ExposureConsumerAreas)

    if args.output_file is None:
        output_filename = country_prefix + '_EPN_ExposureConsumerAreas_withDamage.geojson'
    else:
        output_filename = args.output_file
    # save consumer areas output as geojson file
    save_to_JSON(DamageConsumerAreas, os.path.join(folder_prefix, output_filename))
    #make_histogram(SampleDamageNetwork)

if __name__ == '__main__':
    main()
