#!/usr/bin/env python
# coding: utf-8

#######################################################################################
# This Python script was created by Maria Russo (mrr32@cam.ac.uk); 2022
# the script will perform the following steps:
#     0) Initialise, define functions, parse arguments
#     1) read air pressure and campaign name from flight track files 
#     2) read model variables and Heaviside functions on p levels from hourly pp files
#     3) collocate model variable onto flight track
#     4) write daily netcdf file containing model variables colocated onto flight track
#     5) read daily netcdf files and write monthly netcdf files
#     6) if archiving of hourly pp files is set to false delete hourly pp files before postproc

# This script can be run interactively within a UM suite or offline for postprocessing of UM hourly files
# This script can be used with model data produced with both 360_day and gregorian calendars
# This is controlled by an input argument: choose 'postprocessing' for offline and 'batch' for running within a UM suite 
# How to call the script on the command line: 
# python3 UM_to_flightrack.py -i 'UM_inputdir' -t 'trackdir' -d 'YYYYMM' -n 'n_months' -r 'runid' -p 'pp_stream'  -o 'outdir' -m 'method' -c 'True' postprocessing -s 'select_stash'
# where:
# 'UM_inputdir' = directory containing the UM hourly pp or fieldsfiles
# 'trackdir' = directory containing input flight track netcdf files
# 'YYYYMM' = start date for processing
# 'n_months' = number of months to process (optional, integer, default is 1)
# 'runid' = UM run id (5 characters) eg cm185
# 'pp_stream' = the UM pp_stream containing the required hourly files (1 character) eg f
# 'outdir' = directory to write monthly UM output files colocated onto flight track (optional): 
#            1) for offline jobs default=current directory
#            2) for batch running, default is $DATAM
#               If outdir is present output is also copied to outdir
# 'interpolation_method' = choose between 'lin' (linear) and 'nn' (nearest neighbour); (optional; default=lin)
# -c 'True' produces a model climatology for a small set of flights, e.g. from a field campaign. (optional; default='False')
# 'batch' = if selected the script is expected to run within a UM suite; you can also input: 
#     archive_hourly =  set to 'False' if you do not want to archive hourly UM input files (optional; default=True)
# 'postprocessing' = if selected the script is expected to run outside a UM suite; you can also input:
#     select_stash = stash codes of variables you want to process (optional, default = all variables in input file will be processed)

## EXAMPLES:
# python3 UM_to_flightrack.py -i ~/nethome/data/CIS_TESTS/UM_Input/cm020 -t ~/nethome/data/CIS_TESTS/Flights -d 201001 -n 1 -r cm020 -p l -o ~/data/CIS_TESTS/Out_test postprocessing --select_stash 51001 51009
# python3 UM_to_flightrack.py -i ~/data/CIS_TESTS/UM_Input -t ~/data/CIS_TESTS/Flights -d 201003 -r cm020 -p l -o ~/data/CIS_TESTS/Out_test batch --archive_hourly False 
#######################################################################################

######################################
#   0. Initialise, define functions, parse arguments

#   LOAD REQUIRED APIs 
from datetime import datetime 
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

#########################################################################################################
# Required functions below 

# Generate 8 digit integer code from campaign name ------------------------------------------
def generate_campaign_code(campaign_name):
    import hashlib

    # Make sure all letters are capitalised
    campaign_name = campaign_name.upper()
    campaign_code = int(hashlib.sha1(campaign_name.encode("utf-8")).hexdigest(), 16) % (10 ** 8)

    return campaign_code

# Convert cf variable to cis variable -------------------------------------------------------
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

# Function to read, colocate, write output and remove files if required (this works one month at a time) ---------------------------------------
def process_data_monthly(args,datetag):
    ######
    #print(args)
    ######
    inputdir = args.inputdir
    trackdir = args.trackdir
    cycle_date = datetag
    runid = args.runid
    ppstream = args.ppstream
    method = args.method
    climatology = args.climatology
    outdir = args.outdir
    jobtype = args.jobtype

    # Process parsed arguments and print input variables 
    print(' ')
    print('########## SCRIPT INPUT ############')
    print(' ')
    if inputdir[-1] != '/':
       inputdir=inputdir+'/'
    print('UM data dir = ' + inputdir)

    if trackdir[-1] != '/':
       trackdir=trackdir+'/'
    print('Flight input dir = ' + trackdir)

    print('UM cycle date = ' + cycle_date)

    print('UM run id = ' + runid)

    print('pp stream to read from = ' + ppstream)

    print('Interpolation method = ' + method)

    if climatology == 'True' or climatology == 'true' or climatology == 'TRUE' or climatology == 'T':
        multi_year=True
    else:
        multi_year=False

    if jobtype == 'batch':
        print('This script is running within a UM suite')
        offline=False
        archive_hourly = args.archive_hourly

        # Handling of outdir and additional_outdir (specific to batch jobs)
        if 'outdir' in locals() and outdir != "":
            # If outdir directory is provided, write additional output there
            additional_outdir=outdir
            if additional_outdir[-1] != '/':
                additional_outdir=additional_outdir+'/'
            print('Colocated files are archived and also copied locally to: ', additional_outdir)
        else:
            print('Colocated files are archived but not copied locally')
        # Output is always written to $DATAM (to allow for archiving)
        outdir=inputdir

        # Handilng of conditional arguments
        if archive_hourly == 'True' or archive_hourly == 'true' or archive_hourly == 'TRUE' or archive_hourly == 'T':
            delete_ff=False
            print('Hourly ppstream will not be deleted')
        elif archive_hourly == 'False' or archive_hourly == 'false' or archive_hourly == 'FALSE' or archive_hourly == 'F':
            delete_ff=True
            print('Hourly ppstream will be deleted')
        else:
            delete_ff=False

    elif jobtype == 'postprocessing':
        print('This script is running offline using hourly pp files')
        offline=True

        # Handling of outdir (specific to postprocessing)
        if 'outdir' not in locals():
            outdir = './'
            print ('Colocated files are written to current directory')
        else:
            if outdir[-1] != '/':
                outdir=outdir+'/'
            print ('Colocated files are written to: ',outdir)

        # Handling of conditional arguments
        if args.select_stash != None:
            select_stash=args.select_stash
            print('Only processing selected variables from hourly files: ', select_stash)
        else:
            print('Processing all variables from hourly files')

    ######################################
    print(' ')
    print('########## RUNNING SCRIPT ############')
    print(' ')
    # Initialise variables before daily loop 
    stash_save=[]
    var_save=[]
    campaign_history=[]
    first=True 

    # Find out days within UM cycle for which file track data exists (so we only read and process UM output for those days)
    files=sorted(os.listdir(trackdir))
    if multi_year == True:
        # Find date location in filenames
        filename0=files[0]
        date0="".join([n for n in filename0 if n.isdigit()])
        date_start = filename0.index(date0)
        # Select dates to read flight data (with same month as cycle date)
        flight_files=[filename for filename in files if filename[date_start + 4 : date_start + 6] == cycle_date[4:6]]
        # Select dates with same month and day as flight dates (this will produce output for multiple years for each flight)  
        read_dates = [cycle_date[0:4] + filename[date_start + 4 : date_start + 8] for filename in flight_files]
        flight_dates= [filename[date_start : date_start + 8] for filename in flight_files]
    else:
        # Select dates with same year, month and date as flight dates
        flight_files=[filename for filename in files if cycle_date in filename]
        read_dates=[filename[filename.index(cycle_date):filename.index(cycle_date)+8] for filename in flight_files]
        flight_dates=read_dates

    # Check Daily subdirectory and prepare appropriately
    daily_dir = outdir + 'Daily/'
    if not os.path.exists(daily_dir):
        # Create one if it doesn't exist
        os.makedirs(daily_dir)
    else:
        # Remove existing files (this works if there are existing files and does nothing if there are no files)
        all_files=os.listdir(daily_dir)
        for dfile in all_files:
            os.remove(daily_dir + dfile)
    ######################################

    ###  TIME LOOP #######
    for dd in range(len(read_dates)):
        # Loop through all selected dates
        m_date=read_dates[dd]
        f_date=flight_dates[dd]

        # Test if model data exists for specified date
        input_files=sorted(os.listdir(inputdir))
        test_file=[filename for filename in input_files if (runid in filename and m_date in filename and "a.p" + ppstream in filename)]
        if len(test_file) == 0:
            # Do not read or process data if flight data exists but model data does not.
            print('Model data for ',m_date,' does not exist. Skipping this date')
            if first:
                first_date=read_dates[dd+1]
        else:
            #############~~~~~~~~~~~~~~~~~~~~~
            #   1. READ FLIGHT TRACK DATA
            #MRR add one line below
            first_date=read_dates[0]
            first=False
            # Define flight track filename for dates within the current UM cycle
            trackfile=trackdir + '*' + f_date +'*.nc'
            print('Reading ', trackfile)
            try:
                flight=cis.read_data_list(trackfile,['air_pressure','campaign'])
            except OSError as err:
                # If file does not exists or problems reading it: 
                print("Error: {0}".format(err))
                raise Exception
            else:    
                if multi_year == True:
                    # Convert time to start of flight month 
                    start_of_flight_month="days since "+f_date[0:4]+"-"+f_date[4:6]+"-01"
                    new_time_units_1 = cf_units.Unit(start_of_flight_month, calendar=flight[0].coord("time").units.calendar)
                    flight[0].coord("time").convert_units(new_time_units_1)
                    # Save time array
                    saved_time_data=copy.deepcopy(flight[0].coord("time").data)
                    # Now convert time to start of model month
                    start_of_model_month="days since "+m_date[0:4]+"-"+m_date[4:6]+"-01"
                    new_time_units_2 = cf_units.Unit(start_of_model_month, calendar=flight[0].coord("time").units.calendar)
                    flight[0].coord("time").convert_units(new_time_units_2)
                    # Replace time array with previously saved time array
                    flight[0].coord("time").data=saved_time_data

                # Convert CIS time units to common starting point
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
            # Specify UM model output so only hourly fields that we want to colocate go into selected pp files

            # Define input filename 
            infile=inputdir+runid+'a.p'+ppstream+m_date+'*'
            print('Reading', infile)
            if 'select_stash' in locals():
                # Read only selected variables (this produces a cf.Fieldlist)

                try:
                    reading_vars=cf.read(infile,select=['stash_code='+name for name in select_stash])
                except OSError as err:
                    print("Error: {0}".format(err))
                    raise Exception
                else:
                    read_h51=['true' for name in select_stash if name[0:2] == '51' or name[0:2] == '52']
                    read_h30=['true' for name in select_stash if name[0:2] == '30']
                    # Read Heaviside step functions
                    if len(read_h51) >= 1:
                        heaviside_51=cf.read(infile,select='stash_code=51999')[0] #cf.field.Field
                    if len(read_h30) >= 1:
                        heaviside_30=cf.read(infile,select='stash_code=30301')[0] #cf.field.Field
            else:
                # Reading all variables in the pp stream (this produces a cf.Fieldlist)
                try:
                    reading_vars=cf.read(infile)
                except OSError as err:
                    # If file does not exists or problems reading it:
                    print("Error: {0}".format(err))
                    raise Exception
                else:
                    # Extract Heaviside step functions from reading_vars (if heaviside is not there will return an empty list)
                    heaviside_51 = [var for var in reading_vars if var.get_property("um_stash_source") == "m01s51i999"][0]
                    heaviside_30 = [var for var in reading_vars if var.get_property("um_stash_source") == "m01s30i301"][0]
            #############=====================

            #############+++++++++++++++++++++
            #   3. PROCESS AND COLOCATE (ONE STASHCODE AT THE TIME)

            ###  VARIABLES LOOP ####
            for var in reading_vars:
                # Loop through all stashcode fields in daily variable
                print('Reading ',var.get_property("um_stash_source"))

                # Check if field is a Heaviside function and only process field if not heaviside
                if (var.get_property("um_stash_source") != "m01s51i999" and var.get_property("um_stash_source") != "m01s30i301"):
                    print('Processing ',var.get_property("um_stash_source"))

                    # For section 51 and 52
                    if var.get_property('um_stash_source')[0:6] == 'm01s51' or var.get_property('um_stash_source')[0:6] == 'm01s52':
                        # Check that the appropriate Heaviside function has been read
                        if 'heaviside_51' in locals():
                            print('Dividing field by Heaviside step function')
                            var = var/heaviside_51
                        else:
                            raise Exception('Heaviside function is required for section 51 or 52: add 51999 to UM output')

                    # For section 30
                    if var.get_property('um_stash_source')[0:6] == 'm01s30':
                    # Check that the appropriate Heaviside function has been read
                        if 'heaviside_30' in locals():
                            print('Dividing field by Heaviside step function')
                            var = var/heaviside_30
                        else:
                            raise Exception('Heaviside function is required for section 30: add 30301 to output')

                    # Move data to cis variable format (cf.field.Field to cis.data_io.gridded_data.GriddedData)
                    try:
                        cisvar=cis_from_cf(var)
                    except BaseException as err:
                        # If file does not exists or problems reading it: 
                        print("Error: {0}".format(err))
                        raise Exception
                
                    # Convert working_var_cis to same time units as flight data
                    cisvar.coord("time").convert_units(new_time_units)
                    # Collocate
                    trackvar=cisvar.collocated_onto(flight[0], how=method)   #cis.data_io.ungridded_data.UngriddedDataList
                    flight[0]=trackvar[0]
                    #############+++++++++++++++++++++

                    #############---------------------
                    #   4. WRITE TEMPORARY DAILY OUTPUT
                    # Output files: one file per day for each variable
                    # Define stashcodes for writing monthly files later
                    stash=cisvar.var_name[4:6] + cisvar.var_name[7:10]
                    if m_date == first_date: #read_dates[0]:
                        stash_save.append(stash)
                        var_save.append(cisvar.var_name)

                    # Define output filename for daily files
                    outfile=daily_dir + runid + '_' + m_date + '_stash' + stash + '_flight_track.nc' # one file per day

                    try:
                        flight.save_data(outfile)
                    except BaseException as err:
                        # If file cannot be written: 
                        print("Error: {0}".format(err))
                        raise Exception

                    #############---------------------

    #############@@@@@@@@@@@@@@@@@@@@@
    #   5. WRITE MONTHLY OUTPUT
    # Read daily files and write one file per month for each variable
    # If successful delete daily files

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
                raise Exception
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
                cmip6_filename= 'atmos_' + runid + 'a_1h_' + cycle_date + '01-' + next_month + '01_' + method +'_stash'+stash_save[nv] + '.nc'
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
                                raise Exception
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
    # If requested, remove hourly pp stream before archiving
    all_files=sorted(os.listdir(inputdir))
    previous_month = (datetime.strptime(cycle_date, "%Y%m") + relativedelta(months=-1)).strftime("%Y%m")
    ppstream_cycle_files=[file for file in all_files if 'a.p'+ppstream in file and cycle_date in file]
    previous_month_file=[file for file in all_files if 'a.p'+ppstream in file and previous_month in file]
    if jobtype == 'batch' and delete_ff == True:
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

# End of functions
#########################################################################################################


######## MAIN PROGRAMM ############
# Argument handling --------------
# Create the parser
parser=argparse.ArgumentParser()
parser.add_argument('-i','--inputdir',required=True,type=str,help='Input directory containing hourly pp files')
parser.add_argument('-t','--trackdir',required=True,type=str,help='Directory with input files to colocate onto')
parser.add_argument('-d','--cycle_date',required=True,type=str,help='Date tag to identify start period for analysis')
parser.add_argument('-n','--n_months',default=1,type=int,help='Number of months to process')
parser.add_argument('-r','--runid',required=True,type=str,help='UM job id')
parser.add_argument('-p','--ppstream',required=True,type=str,help='ppstream containing hourly data')
parser.add_argument('-m','--method',type=str,choices=['lin','nn'],default='lin',
        help='Interpolation method, choose between linear and nearest neighbour')
parser.add_argument('-c','--climatology',type=str,default='False',
        help='Model calendar, choose between 360_day and gregorian')
parser.add_argument('-o','--outdir',type=str,help='Output directory for model output on flight track')

# Create subparsers for jobtype
subparser = parser.add_subparsers(dest='jobtype')
# Add subparsers
batch=subparser.add_parser('batch')
postproc=subparser.add_parser('postprocessing')
# Add conditional arguments for batch job
batch.add_argument('-a','--archive_hourly',type=str,default='True',help='logical to archive hourly UM fieldfiles')
# Add conditional arguments for postprocessing job
postproc.add_argument('-s','--select_stash',type=str,nargs='+',help='Optional: only process selected stashcodes')

# Parse the arguments
args=parser.parse_args()

# Identify start date and number of months to process
start_date=args.cycle_date
n_months=args.n_months

# Loop through months to process (default is one)
for nm in range(n_months):
    # Calculate date tag (YEARMONTH) for month to be processed
    datetag=(datetime.strptime(start_date, "%Y%m") + relativedelta(months=nm)).strftime("%Y%m")
    # Call function to process UM data for specified month
    process_data_monthly(args,datetag)

##### END MAIN ####################


