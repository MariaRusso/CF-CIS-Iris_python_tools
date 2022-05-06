# CF-CIS-Iris_python_tools

This repository contains open source code to process model data.  
A brief description of the directories is given here  
(for more details refer to the README file in each directory).

- *Modules*: contains functions to convert cf variables into 
cis, iris or xarray format. CF reads large files much faster than    
iris or cis. These functions allow you to use CF for fast reading   
but then converting to the appropriate format to interface with 
existing scripts using cis, iris or xarray.

- *UM_flight*: contains python routines and functions to enable  
the output of the *Unified Model*, *UKCA* or *UKESM* variables  
on specified aircraft flight tracks
