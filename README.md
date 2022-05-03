# CF-CIS-Iris_python_tools

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

This currently works for gridded variables only and does not consider
auxilliary coordinates.
