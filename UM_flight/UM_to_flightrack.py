#!/usr/bin/env python
# coding: utf-8

#######################################################################################
# This Python script was created by Maria Russo (mrr32@cam.ac.uk); 2022
# the script will perform the following steps:
#     1) read air pressure and campaign name from flight track files 
#     2) read model variables and Heaviside functions on p levels from hourly pp files
#     3) collocate model variable onto flight track
#     4) write daily netcdf file containing model variables colocated onto flight track
#     5) read daily netcdf files and write monthly netcdf files
#     6) if archiving of hourly pp files is set to false delete hourly pp files before postproc

# This script can be run interactively within a UM suite or offline for postprocessing of UM hourly files
# This is controlled by the 'batch_job' logical argument: set to False for offline and True for running within a UM suite 

# How to call the script on the command line: 
# python3 UM_to_flightrack.py -i 'UM_inputdir' -t 'trackdir' -d 'YYYYMM' -r 'runid' -p 'pp_stream' -b 'batch_job' -k 'keep_hourly_files' -o 'outdir'
# where:
# 'UM_inputdir' = directory containing the UM hourly pp or fieldsfiles
# 'trackdir' = directory containing input flight track netcdf files
# 'YYYYMM' = date (year-month) pointing to specific month to be processed
# 'runid' = UM run id (5 characters) eg cm185
# 'pp_stream' = the UM pp_stream containing the required hourly files (1 character) eg f
# 'outdir' = directory to write monthly UM output files colocated onto flight track (optional): 
#            1) for offline jobs default=current directory
#            2) for batch running, output is written to $DATAM for subsequent archiving.
#               If outdir is present output is also copied to outdir (default = files only written to $DATAM)
# 'batch_job' = logical: set to False for running outside UM suite' (optional; default=True)
# 'keep_hourly_files' = logical: set to False if you want to delete hourly UM input files (optional; default=True)
#######################################################################################


#############
#   LOAD REQUIRED APIs and read variables from command line arguments
from datetime import datetime #, timedelta
from dateutil.relativedelta import relativedelta
import numpy as np
import hashlib
import argparse
import cf
from iris import time, cube
import cis
import cf_units
import copy
import os
import sys

def generate_campaign_code(campaign_name):
    import hashlib

    # Make sure all letters are capitalised
    campaign_name = campaign_name.upper()
    campaign_code = int(hashlib.sha1(campaign_name.encode("utf-8")).hexdigest(), 16) % (10 ** 8)

    return campaign_code

def cis_from_cf(cfvar):
    import numpy as np
    import cf
    import iris
    from iris import time, cube
    import cis
    from cis import data_io
    from cis.data_io import gridded_data
    from cis.data_io.gridded_data import GriddedData

    # If cfvar is a cf.fieldlist extract the first field; if not leave as it is
    try:
        len(cfvar)
        cfvar=cfvar[0]
    except TypeError:
        cfvar=cfvar

    # Find number of dimension coordinates in cfvar
    n_dim = np.shape(cfvar.dimension_coordinates())[0]
    # Initialise list of dimension coordinates: this will be filled with coords in iris construct)
    coords_and_dims=[]
    # Loop through dimension coordinates
    for nd in range(n_dim):
        string='dimensioncoordinate'+str(nd)
        dim_array=cfvar.dimension_coordinate(string).array
        dim_name=cfvar.dimension_coordinate(string).standard_name
        dim_units=cfvar.dimension_coordinate(string).units
        iris_coord=iris.coords.DimCoord(dim_array, standard_name=dim_name, units=dim_units)
        iris_dim=(iris_coord,nd)
        coords_and_dims.append(iris_dim)

    # Create CIS gridded data
    data=cfvar.data.array
    s_name=None
    if cfvar.has_property('standard_name'):
        s_name=cfvar.get_property('standard_name')
    l_name=None
    if cfvar.has_property('long_name'):
        l_name=cfvar.get_property('long_name')
    v_name=None
    if cfvar.has_property('um_stash_source'):
        v_name=cfvar.get_property('um_stash_source')
    units=None
    if cfvar.has_property('units'):
        units=cfvar.get_property('units')

    cisvar=GriddedData(data=data, standard_name=s_name, long_name=l_name, var_name=v_name,
                       units=units, dim_coords_and_dims=coords_and_dims)

    return cisvar

# Set defaults for optional arguments 
outdir=None
batch_job=True
keep_hourly=True
# Create the parser with required arguments
parser=argparse.ArgumentParser()
parser.add_argument('-i','--indir',required=True,type=str,help='Input directory containing hourly pp files')
parser.add_argument('-t','--trackdir',required=True,type=str,help='Directory with input files to colocate onto')
parser.add_argument('-d','--date',required=True,type=str,help='Date tag to identify files from one UM cycle')
parser.add_argument('-r','--runid',required=True,type=str,help='UM job id')
parser.add_argument('-p','--ppstream',required=True,type=str,help='ppstream containing hourly data')
parser.add_argument('-b','--batch_job',type=str,help='logical set to True for running within UM suite')
parser.add_argument('-k','--keep_hourly',type=str,help='logical to archive hourly UM fieldfiles')
parser.add_argument('-o','--outdir',type=str,nargs='*',help='Output directory for model output on flight track')
# Parse the arguments
args=parser.parse_args()
inputdir = args.indir
trackdir = args.trackdir
cycle_date = args.date
jobid = args.runid
ppstream = args.ppstream
batch_job = args.batch_job
keep_hourly = args.keep_hourly
outdir = args.outdir

# Test parsing has worked correctly
if inputdir[-1] != '/':
   inputdir=inputdir+'/'
print('UM data dir = ' + inputdir)

if trackdir[-1] != '/':
   trackdir=trackdir+'/'
print('Flight input dir = ' + trackdir)

print('UM cycle date = ' + cycle_date)

print('UM job id = ' + jobid)

print('pp stream to read from = ' + ppstream)

if batch_job == 'True' or batch_job == 'true' or batch_job == 'TRUE' or batch_job == 'T':
    offline=False
    print('This script is running within a UM suite')
elif batch_job == 'False' or batch_job == 'false' or batch_job == 'FALSE' or batch_job == 'F':
    offline=True
    print('This script is running offline using hourly pp files')
else:
    offline=False

if keep_hourly == 'True' or keep_hourly == 'true' or keep_hourly == 'TRUE' or keep_hourly == 'T':
    delete_ff=False
    print('Hourly ppstream will not be deleted')
elif keep_hourly == 'False' or keep_hourly == 'false' or keep_hourly == 'FALSE' or keep_hourly == 'F':
    delete_ff=True
    print('Hourly ppstream will be deleted')
else:
    delete_ff=False

if offline == True:
    # Running from command line in the terminal (outside UM suite)
    if outdir == None:
        outdir = './'
        print ('Colocated files are written to current directory')
    else:
        outdir=outdir[0]
        if outdir[-1] != '/':
            outdir=outdir+'/'
        print ('Colocated files are written to: ',outdir)
else: 
    # Running within UM suite: always write output to $DATAM
    outdir=inputdir
    if outdir != [""] and outdir != None:
        # If outdir directory is provided, write additional output there
        additional_outdir=outdir[0]
        if additional_outdir[-1] != '/':
            additional_outdir=additional_outdir+'/'
        print('Colocated files are archived and also copied locally to: ', additional_outdir)
    else:
        print('Colocated files are archived but not copied locally')

######################################
#   0. Initialise variables before daily loop

# Initialise variables (needed to write monthly files at the end)
stash_save=[]
var_save=[]
campaign_history=[]
# Find out days within UM cycle for which file track data exists (so we only read and process UM output for those days)
flight_dates=sorted(os.listdir(trackdir))
flight_dates=[filename for filename in flight_dates if cycle_date in filename]
read_dates=[filename[filename.index(cycle_date):filename.index(cycle_date)+8] for filename in flight_dates]
# Check Daily subdirectory and prepare appropriately
daily_dir = outdir + 'Daily/'
if not os.path.exists(daily_dir):
    # Create one if it doesn't exist
    os.mkdir(daily_dir)
else:
    # Remove existing files (this works if there are existing files and does nothing if there are no files)
    all_files=os.listdir(daily_dir)
    for dfile in all_files:
        os.remove(daily_dir + dfile)
######################################

###  DAILY LOOP ######################
# Loop through all dates for which a flight track exists and process UM output day by day
for date in read_dates:
    #############~~~~~~~~~~~~~~~~~~~~~
    #   1. READ FLIGHT TRACK DATA
    # Define flight track filename for dates within the current UM cycle
    trackfile=trackdir + '*' + date +'*.nc'
    print('Reading '+trackfile)          
    try:
        flight=cis.read_data_list(trackfile,['air_pressure','campaign'])
    except OSError as err:
        # If file does not exists or problems reading it: 
        print("Error: {0}".format(err))
        raise exception
    else:    
        # Convert CIS time units
        new_time_units = cf_units.Unit("days since 1900-01-01", calendar=flight[0].coord("time").units.calendar)
        flight.coord("time").convert_units(new_time_units)
        # Convert campaigns name to 8 digit integer (campaign code) and write into history metadata)
        c_data = copy.deepcopy(flight[1].data)
        c_names = list(dict.fromkeys(c_data))  # remove duplicates from list
        c_codes = [generate_campaign_code(name) for name in c_names]
        campaigns=dict(zip(c_names, c_codes))
        # convert campaign names to campaign codes and list to numpy array
        flight[1].data=np.array([campaigns[data] for data in c_data])
        # Save campaign information to add to monthly files later
        campaign_history.append(str(campaigns))
    #############~~~~~~~~~~~~~~~~~~~~~

    #############=====================
    #   2. READ MODEL DATA
    # Define input filename
    # Specify UM model output so only hourly fields that we want to colocate go into selected pp files

    # The line below only reads UM files for days which also have a flight track file 
    infile=inputdir+jobid+'a.p'+ppstream+date+'*' 
    print('Reading', infile)  # temporary

    # Reading all variables in the pp stream (this produces a cf.Fieldlist)
    try:
        reading_vars=cf.read(infile)
    except OSError as err:
        # If file does not exists or problems reading it:
        print("Error: {0}".format(err))
        raise exception
    else:
        # Extract Heaviside step functions from reading_vars (if heaviside is not there will return an empty list)
        heaviside_51 = [var for var in reading_vars if var.get_property("um_stash_source") == "m01s51i999"]
        heaviside_30 = [var for var in reading_vars if var.get_property("um_stash_source") == "m01s30i301"]
    #############=====================

    #############+++++++++++++++++++++
    #   3. PROCESS AND COLOCATE (ONE STASHCODE AT THE TIME)
    for var in reading_vars:
	# Loop through all stashcode fields in daily variable
        print('Reading ',var.get_property("um_stash_source"))

        # Check if field is a Heaviside function and only process field if not heaviside
        if (var.get_property("um_stash_source") != "m01s51i999" and var.get_property("um_stash_source") != "m01s30i301"):
            print('Processing ',var.get_property("um_stash_source"))

            # For section 51 and 52
            if var.get_property('um_stash_source')[0:6] == 'm01s51' or var.get_property('um_stash_source')[0:6] == 'm01s52':
		# Check that the appropriate Heaviside function has been read
                if len(heaviside_51) == 1:
                    print('Dividing field by Heaviside step function')
                    var = var/heaviside_51[0]
                elif len(heaviside_51) == 0:
                    raise Exception('Heaviside function is required for section 51 or 52: add 51999 to UM output')
                elif len(heaviside_51) > 1:
                    raise Exception('More than one Heaviside function is found for section 51 or 52')

            # For section 30
            if var.get_property('um_stash_source')[0:6] == 'm01s30':
		# Check that the appropriate Heaviside function has been read
                if len(heaviside_30) == 1:
                    print('Dividing field by Heaviside step function')
                    var = var/heaviside_30[0]
                elif len(heaviside_30) == 0:
                    raise Exception('Heaviside function is required for section 30: add 30301 to output')
                elif len(heaviside_30) > 1:
                    raise Exception('More than one Heaviside function is found for section 30')

	    # Move data to cis variable format (cf.field.Field to cis.data_io.gridded_data.GriddedData)
            try:
                cisvar=cis_from_cf(var)
            except BaseException as err:
                # If file does not exists or problems reading it: 
                print("Error: {0}".format(err))
                raise exception
            else:
	        # Convert working_var_cis to same time units as flight data
                cisvar.coord("time").convert_units(new_time_units)
                trackvar=(cisvar.collocated_onto(flight[0]))   #cis.data_io.ungridded_data.UngriddedDataList
                flight[0]=trackvar[0]
	    #############+++++++++++++++++++++

	    #############---------------------
	    #   4. WRITE TEMPORARY DAILY OUTPUT
	    # Output files: one file per day for each variable
            # Define stashcodes for writing monthly files later
            stash=cisvar.var_name[4:6] + cisvar.var_name[7:10]
            if date == read_dates[0]:
                stash_save.append(stash)
                var_save.append(cisvar.var_name)

	    # Define output filename for daily files
            outfile=daily_dir + jobid + '_' + date + '_stash' + stash + '_flight_track.nc' # one file per day

            try:
                flight.save_data(outfile)
            except BaseException as err:
                # If file cannot be written: 
                print("Error: {0}".format(err))
                raise exception

	    #############---------------------

#############@@@@@@@@@@@@@@@@@@@@@
#   5. WRITE MONTHLY OUTPUT
# Read daily files and write one file per month for each variable
# If successful delete daily files

# Check that the additional output directory exists and make one if not
if 'additional_outdir' in locals() :
    if not os.path.exists(additional_outdir):
        try:
            # Create one if it doesn't exist
            os.mkdir(additional_outdir)
        except BaseException as err:
            # If directory cannot be written:
            print("Error: {0}".format(err))
            raise exception

# Check if any daily files exist for that month
all_daily_files=sorted(os.listdir(daily_dir))
all_daily_files=[daily_dir + file for file in all_daily_files]
if len(all_daily_files) > 0:
    # Read daily files for each stash and write monthly file (one monthly file per stashcode)
    for nv in range(len(stash_save)):
	# Define and read daily files
        daily_files=[file for file in all_daily_files if stash_save[nv] in file]
        try:
            monthly_data=cis.read_data_list(daily_files,[var_save[nv],'campaign'])
        except BaseException as err:
            # If file does not exists or problems reading it: 
            print("Error: {0}".format(err))
            raise exception
        else:
            # Add campaign data to history (after tidying up string)
            campaign_string=str(list(dict.fromkeys(campaign_history))) # remove duplicates
            campaign_string=campaign_string.replace('"', "") # remove unwanted characters
            campaign_string=campaign_string.replace('{', "")
            campaign_string=campaign_string.replace('}', "")
            campaign_string=campaign_string.replace("'", "")
            monthly_data[1].add_history(campaign_string)
            # Calculate date for 1 month after cycle date
            next_month = (datetime.strptime(cycle_date, "%Y%m") + relativedelta(months=1)).strftime("%Y%m")
            # Filename follows CMIP6 naming convention
            cmip6_filename= 'atmos_' + jobid + 'a_1h_' + cycle_date + '01-' + next_month + '01_stash'+stash_save[nv] + '.nc'
            monthly_outfile=outdir + cmip6_filename 
            print(nv, 'Writing data to ', monthly_outfile)
            monthly_data.save_data(monthly_outfile)
            if not offline:
                # Running within UM suite: save monthly files to additional directory if one is specified 
                if 'additional_outdir' in locals() :
                    if not os.path.exists(additional_outdir):
                        try:
                            # Create one if it doesn't exist
                            os.mkdir(additional_outdir)
                        except BaseException as err:
                            # If directory cannot be written:
                            print("Error: {0}".format(err))
                            raise exception
                    # Define filenames
                    additional_monthly_outfile=additional_outdir + cmip6_filename
                    print('Also writing data to ', additional_monthly_outfile)
                    monthly_data.save_data(additional_monthly_outfile)

    # Check if monthly_outfile exists and delete Daily output on flight track
    monthly_files=(os.listdir(outdir))
    if cmip6_filename in monthly_files:
        print('Delete daily files')
        for file in all_daily_files:
            os.remove(file)
    else:
        print('Keeping daily files')
#############@@@@@@@@@@@@@@@@@@@@@

#############*********************
#   6. TIDY UP
# Check logical and remove hourly pp stream before archiving
all_files=sorted(os.listdir(inputdir))
previous_month = (datetime.strptime(cycle_date, "%Y%m") + relativedelta(months=-1)).strftime("%Y%m")
ppstream_cycle_files=[file for file in all_files if 'a.p'+ppstream in file and cycle_date in file]
previous_month_file=[file for file in all_files if 'a.p'+ppstream in file and previous_month in file]
if delete_ff == True:
    print('Delete last day of previous month')
    try:
        for dfile in previous_month_file:
            os.remove(inputdir + dfile)
    except OSError as err:
        print("Error: {0}".format(err))

    print('Delete all files for selected ppstream except last day of month: ',ppstream_cycle_files[-1])
    for file in ppstream_cycle_files[:-1]:
        try:
            os.remove(inputdir+file)
        except OSError as err:
            print("Error: {0}".format(err))
else:
    print('Keeping hourly ppstream')
#############*********************

