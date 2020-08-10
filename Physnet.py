# -*- coding: utf-8 -*-
"""
Created on Thu Jul  2 13:54:57 2020

@author: hfrv2
"""

# -*- coding: utf-8 -*-
"""
Python module for network simulation with physics based models for computing loads.


Created on Thu Jul  2 13:54:57 2020

@author: hfrv2
"""

#import numpy as np
import networkx as nx
import Constants as cons
import random
import os, sys
import pandas as pd
import geopandas as gpd
import pycrs
import numpy
from shapely.geometry import Point, LineString, Polygon
import pypsa
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
plt.style.use("bmh")
%matplotlib inline
import math

def time_stamps(day):
    # generates csv file with timestamps for power flow analysis
    # Generate time snapshots for analysis in Pypsa
    snapshots=pd.DataFrame(columns=['name', 'weightings'])
    dates = pd.date_range(start=day, periods = 24, freq='H').strftime('%d/%m/%Y %H:%M')
    #print(dates)

    snapshots['name'] = dates
    snapshots['weightings'] = 1
    snapshots.to_csv(pypsafp +'/snapshots.csv', index=False)
    return dates

def parser(date_string):
    # Function to read of a specific format from a csv or excel file
    return datetime.strptime(date_string, "%Y %m %d").strftime("%Y %m %d %H:%M" )
   
def parser_load(date_string):
    # Function to read of a specific format from a csv or excel file
    return datetime.strptime(year+'Dec'+date_string, '%Y%b%d %H')  

day='2017-12-28'
year = '2017'
dates=time_stamps(day)    
pypsafp = 'valparaiso-with-load-gen-trafos' # path to save input data for Pypsa
if not os.path.exists(pypsafp):
    os.mkdir(pypsafp)

busesfp = 'red_chile/buses.csv'
busesexposurefp = 'red_chile/buses_elec_exposure.csv'
allnodesfp = 'red_chile/all_nodes.csv'
linesexposurefp = 'red_chile/C1_EPN_ExposureLines_modified.shp'
linesfp = 'red_chile/lines.csv'
powerplantsconsumersfp = 'red_chile/power_plants_consumers.csv'
loadsfp = 'red_chile/loads.csv'
transformersfp = 'red_chile/transformers.csv'

buses = pd.read_csv(busesfp)
busesexposure = pd.read_csv(busesexposurefp)
allnodes = pd.read_csv(allnodesfp)
lines = pd.read_csv(linesfp)
linesexposure = gpd.read_file(linesexposurefp)
powerplantsconsumers = pd.read_csv(powerplantsconsumersfp)
loads = pd.read_csv(loadsfp)
transformers = pd.read_csv(transformersfp)

#print(allnodes)

folderfp = 'network'
networkfp = 'network/network_valparaiso' # path to save output shape and GeoJSON files

if not os.path.exists(folderfp):
    os.mkdir(folderfp)

historicalfp = 'historical_data_chile' # path to databases of load and generation
historicalgen = pd.read_csv(historicalfp + '/SIC_DIC.csv') # historical generation database

data_generationfp = 'red_chile/Generacion_bruta_horaria_SIC_2017/Diciembre/GeneracionBrutaHorariaSIC_1480923703947556098.xlsx'
data_generation = pd.read_excel(data_generationfp, sep='\s+', parse_dates={'datetime': [1, 2, 3]}, index_col='datetime', date_parser=parser)

data_linesfp = 'red_chile/Measurements_lines/OP171228.xls' # path to data of measurements in transmission lines
data_lines = pd.read_excel(data_linesfp, sheet_name='Tx', header=3) # data of measurements in transmission lines

data_loadfp = 'red_chile/Load_SIC/SIC_DIC.csv' # path to data of hourly load for every month

resultsfp = 'network/results' # path to save results of electrical analysis
if not os.path.exists(resultsfp):
    os.mkdir(resultsfp)
        
def build_buses():
    # Creates file for the buses including buses of electrical analysis, power plants and consumers of exposure file and loads GEOJson and shape formats and creates indiviual shape files
    # Code to convert DataFrame to GeodataFrame adapted from https://www.kaggle.com/learn-forum/111860
    buses_df_geometry = [Point(xy) for xy in zip(buses.x, buses.y)] # Buses as points in space with coordinates of buses
    buses_gdf = gpd.GeoDataFrame(buses, geometry=buses_df_geometry) # Dataframe of buses
    #print(buses_gdf.head())

     # include power plants and consumers
    pp_consumers_point = gpd.GeoDataFrame(columns=powerplantsconsumers.columns.values)
    pp_consumers_geometry_point = []

    pp_consumers_xy = gpd.GeoDataFrame(columns=powerplantsconsumers.columns.values)
    pp_consumers_geometry_xy = []

    for index2, row2 in powerplantsconsumers.iterrows():
        if(pd.isnull(row2['x'])):
            # part adapted from https://stackoverflow.com/questions/50155703/pandas-column-convert-string-to-shapely-point-using-map-function
            pieces = row2['wkt_geom'].split() # splits each record into a list of 3
            x = float(pieces[1].lstrip('(')) # latitude
            y = float(pieces[2].rstrip(')')) # longitude 
            pp_consumers_geometry_point.insert(len(pp_consumers_geometry_point), Point(x,y))
            pp_consumers_point.loc[len(pp_consumers_point)] = row2
        else:          
            pp_consumers_geometry_xy.insert(len(pp_consumers_geometry_xy), Point(row2['x'], row2['y']))
            pp_consumers_xy.loc[len(pp_consumers_xy)] = row2 
    pp_consumers_xy_gdf = gpd.GeoDataFrame(pp_consumers_xy, geometry=pp_consumers_geometry_xy)
    pp_consumers_point_gdf = gpd.GeoDataFrame(pp_consumers_point, geometry=pp_consumers_geometry_point)
    #print(pp_consumers_point_gdf)
    #print(pp_consumers_xy_gdf)

    # include loads
    loads_gdf = gpd.GeoDataFrame(columns = numpy.append(loads.columns.values,'geometry'))
    for index, row in loads.iterrows():
        if(len(buses_gdf.loc[buses_gdf['name'] == row['bus']].index) == 0):
            print('The bus: ', row['bus'], 'does not exist, please add it to the network')
        else:
            idxload = buses_gdf.loc[buses_gdf['name'] == row['bus']].index[0]
            loads_gdf.loc[len(loads_gdf)] = [row['name'], row['bus'], row['Name Node'], 
                                             row['FID'], row['taxonomy'], row['Assumption'], 
                                             row['p_factor'], row['Load Name Database'],
                                             row['Point of Connection'], row['User Type'],
                                             buses_gdf['geometry'][idxload]]
    #print(loads_gdf)

    # crs definition buses
    buses_gdf.crs = "EPSG:32719"
    #print(buses_gdf.crs)
    # re-projection
    buses_gdf = buses_gdf.to_crs("EPSG:4326")    
    # Write the data into the Shapefile
    buses_gdf.to_file(networkfp + '_buses.shp')

    # crs definition power plants and consumers with coordinates y x and y
    pp_consumers_xy_gdf.crs = "EPSG:32719"
    #print(pp_consumers_xy_gdf.crs)
    # re-projection
    pp_consumers_xy_gdf = pp_consumers_xy_gdf.to_crs("EPSG:4326")
    # Write the data into the Shapefile
    #pp_consumers_xy_gdf.to_file(networkfp + '_consumers_generators.shp')

    # crs definition power plants and consumers with coordinates by points
    pp_consumers_point_gdf.crs = "EPSG:4326"
    #print(pp_consumers_point_gdf.crs)
    # Write the data into the Shapefile

    loads_gdf.crs = "EPSG:32719"
    #print(loads_gdf.crs)
    # re-projection
    loads_gdf = loads_gdf.to_crs("EPSG:4326") 
    # Write the data into the Shapefile
    #loads_gdf.to_file(networkfp + '_loads.shp')

    # Join buses of electrical analysis, power plants and consumers of exposure file and loads and export as GeoJSON file
    pp_consumers = pd.concat([pp_consumers_xy_gdf, pp_consumers_point_gdf])
    pp_consumers.to_file(networkfp + '_pp_consumers.shp')
    all_buses_gdf = pd.concat([buses_gdf, pp_consumers_xy_gdf, pp_consumers_point_gdf, loads_gdf])
    del all_buses_gdf['wkt_geom']
    #all_buses_gdf = all_buses_gdf.to_crs("EPSG:32719")
    all_buses_gdf['x'] = all_buses_gdf.geometry.x
    all_buses_gdf['y'] = all_buses_gdf.geometry.y
    all_buses_gdf.to_file(networkfp + '_all_buses.shp')
    all_buses_gdf.to_file(networkfp + '_all_buses.geojson', driver='GeoJSON')
    #print(all_buses_gdf)
    
def straight_lines():
# Create straigth lines from buses coordinates
    buses_gdf = gpd.read_file(networkfp + '_buses.shp')
    lines_gdf = gpd.GeoDataFrame(columns=['name', 'bus0', 'bus1', 'x', 'r', 's_nom', 'Nemotecnico', 'Nombre', 'Nombre Linea', 'Nombre Tramo', 'geometry'])
    for index, row in lines.iterrows():
        #print('row', row.values)
        p0 = row['bus0']
        p1 = row['bus1']

        #print(lines_gdf)
        if(len(buses_gdf.loc[buses_gdf['name'] == p0].index) == 0):
            print("There is no bus:", p0, "Please add the bus to the network")
            sys.exit( 1 )
        else:
            indexp0 = buses_gdf.loc[buses_gdf['name'] == p0].index[0]

            p0geom = buses_gdf['geometry'][indexp0] 

        if(len(buses_gdf.loc[buses_gdf['name'] == p1].index) == 0):
            print("There is no bus:", p1, "Please add the bus to the network")
            sys.exit( 1 )
        else:
            indexp1 = buses_gdf.loc[buses_gdf['name'] == p1].index[0]

            p1geom = buses_gdf['geometry'][indexp1] 
        line = LineString([p0geom, p1geom])

        lines_gdf.loc[index] = [row['name'], row['bus0'], row['bus1'], row['x'], row['r'], row['s_nom'], row['Nemotecnico'], row['Nombre'], row['Nombre Linea'], row['Nombre Tramo'], line]

    #print(lines_gdf)
    # crs definition straight lines
    lines_gdf.crs = "EPSG:32719"
    #print(lines_gdf.crs)
    # re-projection
    lines_gdf = lines_gdf.to_crs("EPSG:4326")    
    # Write the data into the Shapefile
    lines_gdf.to_file(networkfp + '_straigt_lines.shp')        
    
def build_lines():
    # include lines considered in the electrical analysis, the distribution lines of exposure file, additional lines, the lines corresponding to transformers and the conection between generators and substations and Generate a GeoJSON file
    # Generate lines with real path from exposure file with numbers of the "bus0" and "bus1" columns of electrical analysis
    linesexposure_gdf = gpd.GeoDataFrame(columns=numpy.append(linesexposure.columns.values, lines.columns.values))
    pp_consumers = gpd.read_file(networkfp + '_pp_consumers.shp')
    buses_gdf = gpd.read_file(networkfp + '_buses.shp')
    for index, row in lines.iterrows():
        #print('row', row.values)
        p0 = row['bus0']
        p1 = row['bus1']

        p0idx = buses.loc[buses['name'] == p0].index[0]
        p1idx = buses.loc[buses['name'] == p1].index[0]
        p0exposure = buses['FID'][p0idx]
        p1exposure = buses['FID'][p1idx]

        name = row['name']
        lineexists = False
        for index2, row2 in linesexposure.iterrows():     
            if((row2['FROM_ID'] == p0exposure and row2['TO_ID'] == p1exposure) or (row2['FROM_ID'] == p1exposure and row2['TO_ID'] == p0exposure)):           
                lineexists = True
                idxline = index2
                break
        if (lineexists == False):
            print("The line: ", name, " is not in the  Exposure File, please add the line to the network")
            #sys.exit( 1 )
        else:
            linesexposure_gdf.loc[len(linesexposure_gdf)] = numpy.append(linesexposure.iloc[idxline], row)
            # Update FROM and TO values
            fromidx = buses.loc[buses['FID'] == row2['FROM_ID']].index[0]
            toidx = buses.loc[buses['FID'] == row2['TO_ID']].index[0]
            fromname = buses['Name Node'][fromidx]
            toname = buses['Name Node'][toidx]
            linesexposure_gdf.loc[(len(linesexposure_gdf)-1), 'FROM'] = fromname
            linesexposure_gdf.loc[(len(linesexposure_gdf)-1), 'TO'] = toname

    # Distribution lines and additional lines to the lines considered in the electrical analysis
    distlines_gdf = gpd.GeoDataFrame(columns=numpy.append(linesexposure.columns.values, lines.columns.values))
    for index, row in linesexposure.iterrows():
        fidline = row['FID']
        lineincluded = False
        for index2, row2 in linesexposure_gdf.iterrows():
            if(row['FID'] == row2['FID']):
                lineincluded = True
        if(lineincluded == False):
            distlines_gdf.loc[len(distlines_gdf)] = row
            # Update FROM and TO values
            fromidx = allnodes.loc[allnodes['FID'] == row['FROM_ID']].index[0]
            toidx = allnodes.loc[allnodes['FID'] == row['TO_ID']].index[0]
            fromname = allnodes['Name Node'][fromidx]
            toname = allnodes['Name Node'][toidx]
            distlines_gdf.loc[(len(distlines_gdf)-1), 'FROM'] = fromname
            distlines_gdf.loc[(len(distlines_gdf)-1), 'TO'] = toname

    #print(len(linesexposure_gdf))
    #print(len(distlines_gdf))

    # include transformers as lines between primary and secondary side nodes
    transformers_gdf = gpd.GeoDataFrame(columns = numpy.append(transformers.columns.values,['FROM','TO','FROM_ID','TO_ID','Reactance','geometry', 'length_m', 'length']))

    for index, row in transformers.iterrows():
        indexbus0 = buses_gdf.loc[buses_gdf['name'] == row['bus0']].index[0]
        indexbus1 = buses_gdf.loc[buses_gdf['name'] == row['bus1']].index[0]

        bus0 = buses_gdf['geometry'][indexbus0]
        bus1 = buses_gdf['geometry'][indexbus1]

        line = LineString([bus0, bus1])
        #length = line.length
        length = math.acos(math.sin(math.radians(bus0.y))*math.sin(math.radians(bus1.y))+math.cos(math.radians(bus0.y))*math.cos(math.radians(bus1.y))*math.cos(math.radians(bus1.x)-math.radians(bus0.x)))*6371           
        transformers_gdf.loc[len(transformers_gdf)] = [row['name'], row['bus0'], row['bus1'], row['s_nom'], row['x'], row['Nombre'], row['Nemotecnico'], row['FID'], row['taxonomy'], buses_gdf['Name Node'][indexbus0], buses_gdf['Name Node'][indexbus1], buses_gdf['FID'][indexbus0], buses_gdf['FID'][indexbus1], row['x'], line, length*1000, length]
    #print(transformers_gdf)

    # include connection between generators and substations
    generators_lines_gdf = gpd.GeoDataFrame(columns = ['FID', 'FROM', 'TO', 'taxonomy', 'FROM_ID', 'TO_ID', 'Reactance', 'name', 'Nombre', 'geometry', 'length_m', 'length'])
    for index, row in pp_consumers.iterrows():
        if((row['taxonomy'] == "Power Plant") or (row['taxonomy'] == "Power Plant Equivalent")):
            indexbus1 = buses_gdf.loc[buses_gdf['name'] == row['bus']].index[0]            
            bus1 = buses_gdf['geometry'][indexbus1]
            bus0 = row.geometry
            #print(row['Name Node'])
            line = LineString([bus0, bus1])
            
            #length = math.acos(math.sin(math.radians(bus0.y))*math.sin(math.radians(bus1.y))+math.cos(math.radians(bus0.y))*math.cos(math.radians(bus1.y))*math.cos(math.radians(bus1.x)-math.radians(bus0.x)))*6371
            length = line.length
            #print(length)
            generators_lines_gdf.loc[len(generators_lines_gdf)] = [row['FID'], row['Name Node'], buses['Name Node'][indexbus1], 'generation', row['FID'], buses['FID'][indexbus1], row['Synchronou'], row['name'], row['Name Node'], line, length*1000, length]
    
    # Delete lines of generators which are already included in the exposure file
    #print('index', generators_lines_gdf.loc[generators_lines_gdf['FID'] == '267'].index[0])
    generators_lines_gdf.drop(generators_lines_gdf.loc[generators_lines_gdf['FID'] == '267'].index[0], axis=0, inplace=True)
    generators_lines_gdf.drop(generators_lines_gdf.loc[generators_lines_gdf['FID'] == '268'].index[0], axis=0, inplace=True)
    generators_lines_gdf.drop(generators_lines_gdf.loc[generators_lines_gdf['FID'] == '269'].index[0], axis=0, inplace=True)
    generators_lines_gdf.drop(generators_lines_gdf.loc[generators_lines_gdf['FID'] == '270'].index[0], axis=0, inplace=True)
    generators_lines_gdf.drop(generators_lines_gdf.loc[generators_lines_gdf['FID'] == '271'].index[0], axis=0, inplace=True)

    # crs definition lines considered in the electrical analysis
    linesexposure_gdf.crs = "EPSG:4326"
    # Write the data into the Shapefile
    linesexposure_gdf.to_file(networkfp + '_lines_elec_exposure.shp')

    # crs definition Distribution lines
    distlines_gdf.crs = "EPSG:4326"
    # Write the data into the Shapefile
    distlines_gdf.to_file(networkfp + '_lines_distribution.shp')

    # crs definition transformers
    transformers_gdf.crs = "EPSG:4326"
    #print(transformers_gdf.crs) 
    # Write the data into the Shapefile
    transformers_gdf.to_file(networkfp + '_transformers.shp')

    # crs definition lines generators   
    generators_lines_gdf.crs = "EPSG:4326"       
    #print(generators_lines_gdf.crs) 
    # Write the data into the Shapefile
    generators_lines_gdf.to_file(networkfp + '_lines_generators.shp')

    # Join lines considered in the electrical analysis, the distribution lines of exposure file and the lines corresponding to transformers and the conection between generators and substations
    all_lines = pd.concat([linesexposure_gdf, distlines_gdf, transformers_gdf, generators_lines_gdf])
    # Reassing FID
    all_lines['FID'] = numpy.arange(len(all_lines))
    # Write the data into the Shapefile
    #all_lines = all_lines.to_crs("EPSG:32719")
    all_lines.to_file(networkfp + '_all_lines.shp')
    all_lines.to_file(networkfp + '_all_lines.geojson', driver='GeoJSON')

def pypsa_network_files():
    # Generate network csv files for Pyspsa from exposure file including buses, generators, lines, loads and transformers
    
    # buses
    all_buses_gdf = gpd.read_file(networkfp + '_all_buses.shp')
    busesp = all_buses_gdf.loc[(all_buses_gdf['taxonomy'] == 'Boundary Substation') | (all_buses_gdf['taxonomy'] == 'Substation') | (all_buses_gdf['taxonomy'] =='Boundary Tap') | (all_buses_gdf['taxonomy'] == 'Tap')]
    busesp[['name', 'v_nom', 'x', 'y']].to_csv(pypsafp + '/buses.csv', index = False)
    #print(all_buses_gdf)
    # generators
    generatorsp = all_buses_gdf.loc[(all_buses_gdf['taxonomy'] == 'Power Plant') | (all_buses_gdf['taxonomy'] == 'Power Plant Equivalent')]
    generatorsp[['name', 'bus', 'p_nom', 'carrier']].to_csv(pypsafp + '/generators.csv', index=False)

    # lines
    all_lines = gpd.read_file(networkfp + '_all_lines.shp')

    linesp = all_lines.loc[all_lines['taxonomy'] == 'Transmission'].reset_index(drop=True)
    # Drop S/E GNL QUINTERO - S/E CENTRAL QUINTERO FID 69 because it was not considered in the electrical analysis
    linesp.drop(linesp.loc[linesp['FID'] == 69].index[0], axis=0 , inplace=True)
    linesp[['name', 'bus0', 'bus1', 'x', 'r', 's_nom']].to_csv(pypsafp + '/lines.csv', index=False)

    # loads
    loadsp = all_buses_gdf.loc[all_buses_gdf['taxonomy'] == 'Load']
    loadsp[['name', 'bus']].to_csv(pypsafp + '/loads.csv', index = False)

    # transformers
    transformersp = all_lines.loc[all_lines['taxonomy'] == 'Transformer']
    transformersp[['name', 'bus0', 'bus1', 's_nom', 'x']].to_csv(pypsafp + '/transformers.csv', index=False)
 
def gen_time_series():
    # Generate file 'generators-p_set.csv' with the generation in time for each power plant
    # Read excel file with data of hourly generation of power plants for every month
    data_generation['Hour Corrected'] = data_generation['Hour']-1

    data_generation.index += pd.TimedeltaIndex(data_generation['Hour Corrected'], unit='h')
    #print(data_generation.loc['2017-12-01'])
    #print(data_generation.dtypes)
    
    pp_consumers = gpd.read_file(networkfp + '_pp_consumers.shp')
    generators = pp_consumers.loc[(pp_consumers['taxonomy'] == 'Power Plant') | (pp_consumers['taxonomy'] == 'Power Plant Equivalent')]
    generation = pd.DataFrame(columns=numpy.append('name', generators['name']))
    generation['name'] = dates
    # list with the hours in a day from 1 to 24
    hours = list(range(1,25))
    for index, row in generators.iterrows():
    #    print(row['Assumption'])
        if(row['Assumption'] == 0): # Generation taken from the database of gross generation of the power plants
            #print(row['Plant Name'])
            generation[row['name']] = data_generation.loc[data_generation['Central'] == row['Plant Name']].loc[day, 'Generacion_MWh'].values
        if(row['Assumption'] == 1): # Generation calculated as a percentage of generation capacity given by 'p_factor'
            generation[row['name']] = row['p_nom']*row['p_factor']
            #print(row['name'])
        if(row['Assumption'] == 2): # Generation calculated as a percentage of generation capacity given by 'p_factor'. Corresponds to equivalent line with no power measurement available
            generation[row['name']] = row['p_nom']*row['p_factor']
        if(row['Assumption'] == 3): # Generation taken from the power measuremt in line. Corresponds to an equivalent line
            # Select the measurement in the line and choose active power, field corresponding to 'Transferencia por líneas en MW'
            p_measurement = data_lines.loc[(data_lines['Barra i'] == row['Barra i']) & (data_lines['Barra j'] == row['Barra j'])].iloc[0]
            #print(row['Barra i'],row['Barra j'])
            #print(p_measurement)
            #p_measurementac = p_measurement.iloc[0]
            p_measurementdf = pd.DataFrame(p_measurement) 
            generation[row['name']] = p_measurementdf.loc[hours].values * row['p_factor']
        if(row['Assumption'] == 4): # Solar power plant with no data: 85% of capacity during day and 0 during night
            generation.loc[0:5, row['name']] = 0
            generation.loc[6:18, row['name']] = row['p_nom']*row['p_factor']
            generation.loc[19:24, row['name']] = 0    
    # Save file
    generation.to_csv(pypsafp +'/generators-p_set.csv', index=False)
 
    #print(generation)
  
 def load_time_series(year):
    # Generate file 'loads-p_set.csv' with the power consumed in time for each load
    # Translate month names, convert day and hour numbers to integer
    #year = '2017'
    data_load = pd.read_csv(data_loadfp) # Read excel file with data of hourly load for every month
    data_load['DIANUM'] = data_load['DIANUM'].astype(int)
    data_load['HoraDia'] = data_load['HoraDia'].astype(int)-1
    data_load.loc[data_load['MESN'] == 'DIC', 'MESN'] = 'Dec'
    data_load['year'] = year

    dates_l = data_load['year'].map(str) + data_load['MESN'].map(str) + data_load['DIANUM'].map(str) + data_load['HoraDia'].map(str)
    datestime = pd.to_datetime(dates_l.values, format='%Y%b%d%H')
    data_load.index = datestime

    #print(data_load)
    
    #Create dataframe with time series of load
    consumption = pd.DataFrame(columns=numpy.append('name',loads['name']))
    consumption['name'] = dates

    for index, row in loads.iterrows():    
        if(row['Assumption'] > 0):
            consumption[row['name']] = row['Assumption']*row['p_factor']
        if(row['Assumption'] == 0):          
            names_loads = row['Load Name Database'].split('--') # loads in the database corresponding to the analized load in node
            number_loads = len(names_loads)        
            if(number_loads == 1):
                consumption[row['name']] = data_load.loc[(data_load['Barra'] == row['Load Name Database']) & 
                                                         (data_load['PuntoConexion'] == row['Point of Connection']) & 
                                                         (data_load['Tipo2'] == row['User Type'])].loc[day, 'Total'].values

            if(number_loads > 1):
                points_connection = row['Point of Connection'].split('--')
                user_types = row['User Type'].split('--')
                consumption_node = numpy.zeros(24)    
                for i in range(number_loads):        
                    consumption_node += data_load.loc[(data_load['Barra'] == names_loads[i]) & 
                                                      (data_load['PuntoConexion'] == points_connection[i]) & 
                                                      (data_load['Tipo2'] == user_types[i])].loc[day, 'Total'].values
                    consumption[row['name']] = consumption_node             
    #consumption
    consumption.to_csv(pypsafp +'/loads-p_set.csv', index=False)

def pypsa_network():
    # Generate Pyspsa netkork and run a linear power flow
    # Returns updated network with the results of the power flow
    
    network = pypsa.Network("valparaiso-with-load-gen-trafos")

#    network.loads.q_set = 0.7
    network.srid = 32719
    # network.plot() # uncomment to plot network
#     print("Buses:") 
#     print(network.buses)
#     print("Generators:")
#     print(network.generators)
#     print("Loads:")
#     print(network.loads)
#     print("Power in Generators:")
#     print(network.generators_t.p_set)
#     print("Loads:")
#     print(network.loads_t.p_set)
#     print("Transformers:")
#     print(network.transformers)

    #Do a Newton-Raphson power flow
    
    #print(network.buses)
    network.lpf()
    #network.pf(use_seed=True)
    #network.pf()

    print('Power in Lines - MW:')
    print(network.lines_t.p0)

    print('Angle of Voltage in Buses - Degrees:')
    print(network.buses_t.v_ang*180/numpy.pi)

    print('Magnitude of Voltage in Buses - Per Unit:')
    print(network.buses_t.v_mag_pu)

    
# Uncomment below to plot line loading. Also uncoment line #all_buses_gdf = all_buses_gdf.to_crs("EPSG:32719") in function
# build_buses() and line #all_lines = all_lines.to_crs("EPSG:32719") in function build_lines()
#     loading = (network.lines_t.p0.abs().mean().sort_index())
    
#     fig,ax = plt.subplots(
#     figsize=(15,14),
#     subplot_kw={"projection": ccrs.PlateCarree()})
#     #subplot_kw={"projection": ccrs.Orthographic()})    
#     #subplot_kw={"projection": ccrs.Mercator()})
#     #subplot_kw={"projection": ccrs.epsg(32719)})
   
#     network.plot(ax=ax,
#             bus_colors='gray',
#             branch_components=["Line"],
#             line_widths=3,
#             line_colors=loading,
#             line_cmap=plt.cm.viridis,
#             color_geomap=True,
#             bus_sizes=0.0001)
#     plt.rcParams['axes.titlesize'] = 30
#     ax.axis('off');
#     ax.set(title='Line Loading Valparaiso (Chile)')
#     #plt.savefig(resultsfp + '/line_loading.png')

    network.export_to_csv_folder("results_valparaiso-pf")
            
    return network

def lines_measur_results(network):
# returns and save a dataframe wiht:
# 1) Available measurements of power in the lines considered in the electrical analysis in the file
# 'red_chile/Measurements_lines/OP171228.xls' (variable: data_lines).
# 2) the results obtained for the power in the lines in the power flow analysis.
# 3) log10 error between measurements and results.

# network is the Pypsa network after running the power flow
    
    #lines

    results_p0 = numpy.transpose(network.lines_t.p0)
    results_p0.columns = dates + ' P0'

    results_p1 = numpy.transpose(network.lines_t.p1)
    results_p1.columns = dates + ' P1'
    results = pd.concat([results_p0, results_p1], axis=1)
    #print(results)
    errors = pd.DataFrame(columns=['FID','P0 Error Max', 'P0 Error Mean', 'P1 Error Max', 'P1 Error Mean'])
    
    measurements_empty = pd.DataFrame(columns=(dates + ' P0 Measurement').append(dates + ' P1 Measurement').append(dates + ' P0 Error').append(dates + ' P1 Error'))
    measurements = pd.concat([results, measurements_empty, errors], axis=1)
    
    #measurements.set_index = dates
    hours = list(range(1,25))
    
    # Assign available measurements for the power in lines, according to the information given in the file:
    # 'red_chile/lines.csv' (variable 'lines')
    for index, row in lines.iterrows():
        if(~numpy.isnan(row['Type'])):        
            names_bus0 = row['Barra i'].split('--') # bus 0 in the database crresponding to the analized line
            names_bus1 = row['Barra j'].split('--') # bus 1 in the database corresponding to the analized line
            rows = row['Barra j'].split('--')
            number_lines = len(names_bus0)

            names_results = row['Result Line'].split('--')
            direction_results = row['Result Direction'].split('--')

            number_results = len(names_results)

            if(number_lines == 1): # One measurement for the line
                # Select the measurement in the line and choose active power, field corresponding to 'Transferencia por líneas en MW'
                #print(row['Barra i'],row['Barra j'])            
                p_measurement = data_lines.loc[(data_lines['Barra i'] == row['Barra i']) & (data_lines['Barra j'] == row['Barra j'])].iloc[0]

                #print(p_measurement)
                p_measurementdf = pd.DataFrame(p_measurement)
                measurement =  numpy.transpose(p_measurementdf.loc[hours].values).astype(float)           
                #generation[row['name']] = p_measurementdf.loc[hours].values * row['p_factor']            
                if(number_results == 1): # Measurement compared to result in one direction
                    #measurements.set_index = names_results
                    result = measurements.loc[names_results, dates + ' ' + direction_results].values.astype(float)
                    #print(measurement)
                    #print(result)                
                    measurements.loc[names_results, dates + ' ' + direction_results + ' Measurement'] = measurement
                    errors = abs(numpy.log10(abs(measurement)+1)-numpy.log10(abs(result)+1))
                    measurements.loc[names_results, dates + ' ' + direction_results + ' Error'] = errors
                    measurements.loc[names_results, direction_results[0] + ' Error Max'] = errors.max()
                    measurements.loc[names_results, direction_results[0] + ' Error Mean'] = errors.mean()

            if(number_lines  > 1): # Measurements in two directions
                if(row['Used'] != -1): # use measurement for the index given in column 'Used'
                    #print(names_results)
                    #print(direction_results)
                    bus0 = names_bus0[int(row['Used'])]
                    bus1 = names_bus1[int(row['Used'])]
                    #print(bus0,bus1)          
                    p_measurement = data_lines.loc[(data_lines['Barra i'] == bus0) & (data_lines['Barra j'] == bus1)].iloc[0]
                    p_measurementdf = pd.DataFrame(p_measurement)
                    #print(p_measurementdf)
                    measurement =  numpy.transpose(p_measurementdf.loc[hours].values).astype(float)
                    #print(measurement[0])
                    #print(measurements)
                    result = measurements.loc[names_results[0], dates + ' ' + direction_results[0]].values.astype(float)
    
                    if(number_results == 1): # Measurement compared to result in one direction 
                        #print(names_results[0])
                        #print(direction_results[0])
                        measurements.loc[names_results, dates + ' ' + direction_results + ' Measurement'] = measurement
                        #print(measurement)
                        #print(result)
                        errors = abs(numpy.log10(abs(measurement)+1)-numpy.log10(abs(result)+1))
                        measurements.loc[names_results, dates + ' ' + direction_results + ' Error'] = errors
                        measurements.loc[names_results, direction_results[0] + ' Error Max'] = errors.max()
                        measurements.loc[names_results, direction_results[0] +  ' Error Mean'] = errors.mean()

   
                if(row['Used'] == -1): # use all measurements
                    for i in range(number_lines):
                        #print(names_results)
                        #print(direction_results)
                        bus0 = names_bus0[i]
                        bus1 = names_bus1[i]
                        #print(bus0, bus1)
                        p_measurement = data_lines.loc[(data_lines['Barra i'] == bus0) & (data_lines['Barra j'] == bus1)].iloc[0]
                        p_measurementdf = pd.DataFrame(p_measurement)
                        measurement =  numpy.transpose(p_measurementdf.loc[hours].values).astype(float)[0]
                        #print('aca', measurement)
                        result = measurements.loc[names_results[i], dates + ' ' + direction_results[i]].values.astype(float)
                        if(row['Type'] == 1): # two measurements, one for each direction
                            #print(type(names_results[i]))
                            errors = abs(numpy.log10(abs(measurement)+1)-numpy.log10(abs(result)+1))
                            line = names_results[i]
                            measurements.loc[line, dates + ' ' + direction_results[i] + ' Measurement'] = measurement
                            measurements.loc[line, dates + ' ' + direction_results[i] + ' Error'] = errors
                            measurements.loc[line, direction_results[i] + ' Error Max'] = errors.max()
                            measurements.loc[line, direction_results[i] +  ' Error Mean'] = errors.mean()

    #measurements.to_csv(resultsfp + '/measur_results.csv', index=True)
    #measurements
    return measurements          
          
def errors_summary(measurements):
# Create plots of summary of log10 errors of the power in lines of results of load flow analysis with respect to the
# measurements in both directions for the period of analysis

# measurements is the dataframe with measurements, results of power flow and error for the lines generated with the function
# lines_measur_results(network):

    lines_measurements_p0 = []
    lines_measurements_p1 = []
    for index, row in lines.iterrows():
        if(~numpy.isnan(row['Type'])):
            names_results = row['Result Line'].split('--')
            direction_results = row['Result Direction'].split('--')

            number_results = len(names_results)
            #print(names_results)
            #print(direction_results)
            if(number_results == 1):
                if(direction_results[0] == 'P0'):
                    lines_measurements_p0.append(names_results[0])
                    #print('P0', lines_measurements_p0)
                if(direction_results[0] == 'P1'):
                    lines_measurements_p1.append(names_results[0])
            if(number_results > 1):
                for i in range(number_results):
                    if(direction_results[i] == 'P0'):
                        lines_measurements_p0.append(names_results[i])
                    if(direction_results[i] == 'P1'):
                        lines_measurements_p1.append(names_results[i])
    #print(lines_measurements_p0)
    #print(lines_measurements_p1)

    fig, ax = plt.subplots(1,1, figsize=(15,6))
    measurements.loc[lines_measurements_p0, dates + ' P0 Error'].T.plot(ax=ax)
    ax.set(title='Log10 Error of Power in Lines from Bus P0 to P1', xlabel='Date')
    plt.savefig(resultsfp + '/Log10 Error of Power in Lines from Bus P0 to P1.png')
    
    fig, ax1 = plt.subplots(1,1, figsize=(15,4))
    measurements.loc[lines_measurements_p1, dates + ' P1 Error'].T.plot(ax=ax1)
    ax1.set(title='Log10 Error of Power in Lines from Bus P1 to P0', xlabel='Date')
    plt.savefig(resultsfp + '/Log10 Error of Power in Lines from Bus P1 to P0.png')
          
# Merging adapted from https://stackoverflow.com/questions/19125091/pandas-merge-how-to-avoid-duplicating-columns
def update_lines(measurements):
# Adds errors to line file, saves the updated file and returns path of the file
    
# measurements is the dataframe with measurements, results of power flow and error for the lines generated with the function
# lines_measur_results(network):

    
    lines = gpd.read_file(networkfp + '_all_lines.geojson', driver='GeoJSON')
    lines_error = pd.concat([lines, measurements], axis=1)

    cols_additional = measurements.columns.difference(lines.columns)
    lines_error = pd.merge(lines, measurements, left_on='name', right_index=True, how='left', sort=False)
    lines_error.loc[4, 'name']
    lines_error.to_file(resultsfp + '/lines_power_flow.geojson', driver='GeoJSON')
    #lines_error
    #lines.name
    #measurements
    return resultsfp + '/lines_power_flow.geojson'    
          
#def evaluate_system_loads(G):
def evaluate_system_loads():
    # returns the path of the geojson file that contains the information of the lines of the network, the results of the
    # power in the lines obtained in the DC power flow analysis, the available measurements of the power and the error. 
    # Prints results of the power flow analysis and plots the errors.
    build_buses()
    straight_lines()
    build_lines()
    pypsa_network_files()
    gen_time_series()
    load_time_series(year)
    network = pypsa_network()
    lines_measur_results(network)
    measurements = lines_measur_results(network)
    errors_summary(measurements)
    path = update_lines(measurements)
    return path
    #return G
     
    
