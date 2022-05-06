# CF-CIS-Iris_python_tools

<<<<<<< HEAD
This repository contains open source code to process model data. A brief description of the directories is given here (for more details refer to each directory README file).

- *Modules*: contains functions to convert cf variables into 
cis, iris or xarray format. CF reads large files much faster than    
iris or cis. These functions allow you to use CF for fast reading   
but then converting to the appropriate format to interface with 
existing scripts using cis, iris or xarray.
=======
Functions to convert cf variables into cis, iris or xarray format. 
This is useful as it allows to use CF for fast reading of large 
files and then convert to the appropriate format to interface with 
existing scripts using cis, iris or xarray.

The module **convert_CFvar.py** contains three functions:  
- *cis_from_cf*: produces a cis variable from a cf variable  
`from convert_CFvars import cis_from_cf`   
`new_cis_var = cis_from_cf(my_cf_var)`

- *iris_from_cf*: produces a iris variable from a cf variable  
`from convert_CFvars import iris_from_cf`   
`new_iris_var = iris_from_cf(my_cf_var)`

- *xarray_from_cf*: produces an xarray variable from a cf variable  
`from convert_CFvars import xarray_from_cf`   
`new_xarray_var = xarray_from_cf(my_cf_var)`
>>>>>>> de17df27c6af2b8fb9111f1932de37b81458d745

- *UM_flight*: contains python routines and functions supporting  
output of the *Unified Model*, *UKCA* or *UKESM* on specified  
aircraft flighttracks
