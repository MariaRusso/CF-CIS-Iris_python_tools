# CF-CIS-Iris_python_tools

Functions to convert cf variables into iris or cis format. 
This is useful as it allows to use CF for fast reading of large 
files and then convert to the appropriate format to interface with 
existing scripts (which might require a CIS or Iris variable format).

The module **convert_CFvar.py** contains two functions:  
- *cis_from_cf*: produces a cis variable from a cf variable  
new_cis_var = cis_from_cf(my_cf_var)
<<<<<<< HEAD
- *iris_from_cf*: produces a iris variable from a cf variable. 
=======
- *iris_from_cf*: produces a iris variable from a cf variable.  
>>>>>>> 4aaf705add954fbb9113917b1d10e94a1e0a62dc
new_iris_var = iris_from_cf(my_cf_var)

This currently works for gridded variables only and does not consider
auxilliary coordinates.
