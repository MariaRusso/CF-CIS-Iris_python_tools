## CF-CIS-Iris_python_tools

This repository contains open source python code to process model data. A brief description of the directories is given here (for more details refer to the README file in each directory).

### Modules

Contains functions to convert *cf* variables into *cis*, *iris* or *xarray* format.  
Reading large files with *cf.read* is much faster than iris or cis (~6 and ~30 times faster respectively). If you decide to use *cf* for fast reading, these functions can then convert your *cf* structure to the required structure needed by existing scripts (*cis*, *iris* or *xarray*).

### UM_flight

Contains python routines and functions to enable the colocation of gridded variables from existing files (.nc, .pp or UM field files) on specified aircraft flight tracks.  
The main python script (*UM_to_flightrack.py*) can also be embedded into a UM suite to produce monthly netcdf files of selected variables on specified flight tracks, which are archived for later use.  
