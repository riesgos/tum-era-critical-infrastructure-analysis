# System_Reliability
WPS for performing the reliability of infrastructure networks

File content:

-System_Reliability_Analysis: contains the WPS standard and the main function in the run function (wrapper)
-System_Reliability_Analysis_local: the main function that runs locally (i.e. without interaction with other WPS), for testing the functions that run within the WPS.
-Sysrel: Python module, contains functions called by the main function
-Netsim: Python module, contains functions called by the Sysrel module
-Constants: constants used by the previously mentioned modules
-geojson files for exposure in Chile: nodes and lines (areas pending), and in Ecuador: nodes, lines and areas
-geojson files of nodes with probability of failure (mimicking output from damage WPS). Used for testing the code locally
-json file for network fragility taxonomy
