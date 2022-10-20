#!/usr/bin/env python
# coding: utf-8

#######################################################################################
# This Python script was created by Maria Russo (mrr32@cam.ac.uk); 2022
# This Python script will:
#     1) read aircraft data (netcdf file)
#     2) write aircraft data in a consistent format to Aerocom flight data

#######################################################################################

#############
#   0. LOAD REQUIRED APIs and assign variables with info passed from rose/cylc
#import iris
from iris import time, cube
import cis
from cis.data_io.ungridded_data import UngriddedData, Metadata
from cis.time_util import PartialDateTime
import cf_units
import copy
import os
import numpy as np

# PLEASE EDIT BELOW ACCORDING TO SPECIFIC DATA TO PROCESS #################################################
# Directory containing raw aircraft data
inputdir='/home/mrr32/data/CIS_TESTS/Flight_raw/FAAM/ACSIS/' 
# Directory to write processed aircraft data to
outdir='/home/mrr32/data/CIS_TESTS/Flights/'       
# The campaign name is added as a string variable to the processed netcdf file
campaign_name='CONSTRAIN'   
# Variable names for standard coordinates to read from raw aircraft files
varnames=['Time','LAT_GIN','LON_GIN','ALT_GIN','PS_RVSM'] 
# Position of date characters in the filename string
date_position=slice(10,18) 
# Plugin to be used by cis.read_data
aircraft_plugin='FAAM'
###########################################################################################################


#############~~~~~~~~~~~~~~~~~~~~~
# import the correct plugin for the chosen aircraft data 
import_plugin = "from flight_py_tools import " + aircraft_plugin
exec(import_plugin)
from patch_cis_string_data import UngriddedData_string
from patch_cis_string_data import write_coordinate_list, add_data_to_file

# Find out days for which aircraft data exists
inputdir=inputdir
files=sorted(os.listdir(inputdir))
dates=[file[date_position] for file in files] 
# Remove duplicate dates (if more than one flight per day)
dates=list(dict.fromkeys(dates))
for date in dates:
    trackfile=inputdir + '*' + date + '*.nc'
    rdata=cis.read_data_list(trackfile,varnames,product=aircraft_plugin)

    # Convert FAAM time coordinates from seconds to days and from int to float
    # First set out new units to convert to
    new_time_units = cf_units.Unit("days since 1600-01-01 00:00:00", calendar=rdata.coord("time").units.calendar)

    for n in range(len(varnames)):
        rdata[n].coord("time").data = rdata[n].coord("time").data.astype(np.float64)  # convert time coord
        rdata[n].coord("time").convert_units(new_time_units)

    campaign_data=np.array([campaign_name] * len(rdata[0].data))

    campaign_metadata = Metadata(name='campaign', standard_name=None, long_name='campaign', history='', units='unkown')
    new_coords = rdata[0]._coords
    
    campaign_var = UngriddedData_string(data=campaign_data, metadata=campaign_metadata, coords=new_coords)
        
    # Replace first variable with campaign_var
    rdata[0]=campaign_var

    # Write to output file
    outfile=outdir+'core_faam_'+date+'_'+campaign_name+'.nc'
    rdata[0].save_data(outfile)


