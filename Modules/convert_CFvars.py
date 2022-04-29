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

############################################################################################

def iris_from_cf(cfvar):
    import numpy as np
    import cf
    import iris
    from iris import time, cube
    from iris.cube import Cube
    #import cis
    #from cis import data_io
    #from cis.data_io import gridded_data
    #from cis.data_io.gridded_data import GriddedData

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

    # Create IRIS gridded data
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

    irisvar=iris.cube.Cube(data, standard_name=s_name, long_name=l_name, var_name=v_name,
            units=units, dim_coords_and_dims=coords_and_dims)

    return irisvar

################################################################################################

def xarray_from_cf(cfvar):
    import numpy as np
    import cf
    import iris
    from iris import time, cube
    from iris.cube import Cube
    import xarray
    from xarray.convert import from_iris

    # If cfvar is a cf.fieldlist extract the first field; if not leave as it is
    try:
        len(cfvar)
        cfvar=cfvar[0]
    except TypeError:
        cfvar=cfvar

    cube=iris_from_cf(cfvar)
    xarray=from_iris(cube)

    return xarray

################################################################################################
