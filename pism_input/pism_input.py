
import os, sys
import numpy as np
from pyproj import Proj
from optparse import OptionParser
import netCDF4
from netCDF4 import Dataset as CDF
import jinja2


def get_path_to_data(output_data_path,dataset,resolution,
                     custom_data_path=None, custom_file_name=None):

    if custom_data_path is None:
        custom_data_path =  os.path.join(output_data_path, dataset)

    if custom_file_name is None:
        custom_file_name = dataset+"_"+str(resolution)+"km.nc"

    path_to_data = os.path.join(custom_data_path,custom_file_name)

    return path_to_data


def write_regrid_submission_file(config, data_path, dataset, inputfile, resolution,
                                 cdo_targetgrid_file, regridded_file, use_conservative_regridding):

    """ This writes a SLURM submission file for the CPU heavy task of regridding.
        Regridding is done via CDO.
        Output: a cdo_remap.sh file in the folder you executed this function.
    """

    # make jinja aware of templates in the pism_input/tools folder
    jinja_env = jinja2.Environment(loader=jinja2.FileSystemLoader(
                searchpath=os.path.join(config.project_root,"tools")))

    scen_template_file = "GENERATED_SCENARIO.SCEN.template"
    scen_template = jinja_env.get_template("cdo_remap.sh.template")
    mapweights = os.path.join(data_path, "mapweights_"+str(resolution)+"km.nc")

    out = scen_template.render(user=config.username,
                               use_conservative_regridding = use_conservative_regridding,
                               targetgrid = cdo_targetgrid_file,
                               inputfile = inputfile,
                               mapweights = mapweights,
                               regridded_file = regridded_file,
                              )

    with open("cdo_remap.sh", 'w') as f:
        f.write(out)
    print "Wrote cdo_remap.sh, submit with sbatch cdo_remap.sh to compute nodes."


def create_grid_for_cdo_remap(path_to_write,use_PISM_grid=True,resolution=15.0):

    """
    Create a netcdf file holding the target grid for cdo in folder path_to_write.
    If use_PISM_grid is true, we create the grid exactly as PISM would do it,
    so PISM will not regrid internally when such grid is used as input.
    Else, and for any resolution, grid spacing will be calculated by equal
    distancing between domain boundaries.
    """

    # hardcoded grids, as inferred from PISM output
    PISM_grid_set={}
    PISM_grid_set[50]=[120,120,-2777500,-2777500,3172500,3172500]
    PISM_grid_set[30]=[200,200,-2787500,-2787500,3182500,3182500]
    PISM_grid_set[20]=[300,300,-2792500,-2792500,3187500,3187500]
    PISM_grid_set[15]=[400,400,-2795000,-2795000,3190000,3190000]
    PISM_grid_set[12]=[500,500,-2796500,-2796500,3191500,3191500]
    PISM_grid_set[10]=[600,600,-2797500,-2797500,3192500,3192500]
    PISM_grid_set[7]=[800,800,-2798750,-2798750,3193750,3193750]
    PISM_grid_set[5]=[1200,1200,-2800000,-2800000,3195000,3195000] # 5km standard Albmap input grid
    PISM_grid_set[3]=[2000,2000,-2801000,-2801000,3196000,3196000]
    PISM_grid_set[2]=[3000,3000,-2801500,-2801500,3196500,3196500]
    PISM_grid_set[1]=[6000,6000,-2802000,-2802000,3197000,3197000]

    ## create output directory if it does not exist.
    if not os.path.exists(path_to_write): os.makedirs(path_to_write)
    nc_outfile = os.path.join(path_to_write,'grid_'+str(int(resolution))+'km.nc')
    grid_spacing = resolution * 1.e3 # convert km -> m

    # define output grid
    de = dn =  grid_spacing # m

    if (use_PISM_grid):

        dxy = int(grid_spacing/1e3)

        e0 = PISM_grid_set[int(dxy)][2]
        n0 = PISM_grid_set[int(dxy)][3]

        M = PISM_grid_set[int(dxy)][0]
        N = PISM_grid_set[int(dxy)][1]

        easting  = np.linspace(e0, e0+(M-1)*de, M)
        northing = np.linspace(n0, n0+(N-1)*dn, N)

        print PISM_grid_set[int(dxy)]

        nc_outfile = nc_outfile.replace('grid_', 'pism_')

    else:

        dxy_albmap=5
        e0 = PISM_grid_set[dxy_albmap][2]
        n0 = PISM_grid_set[dxy_albmap][3]
        e1 = PISM_grid_set[dxy_albmap][4]
        n1 = PISM_grid_set[dxy_albmap][5]
        #x: -2800000 to 3195000
        #y: -2800000 to 3195000

        M = int((e1 - e0)/de) + 1
        N = int((n1 - n0)/dn) + 1

        easting  = np.linspace(e0, e1, M)
        northing = np.linspace(n0, n1, N)


    ee, nn = np.meshgrid(easting,northing)

    print "Grid is created for "+str(int(resolution))+"km resolution:"
    print M,N,easting[0],northing[0],easting[-1],northing[-1],np.diff(easting)[0],np.diff(northing)[0]

    #projection = "+proj=stere +ellps=WGS84 +datum=WGS84 +lon_0=0 +lat_0=-90 +lat_ts=-71 +units=m"
    projection = "+lon_0=0.0 +ellps=WGS84 +datum=WGS84 +lat_ts=-71.0 +proj=stere +x_0=0.0 +units=m +y_0=0.0 +lat_0=-90.0"


    proj = Proj(projection)

    lon,lat = proj(ee,nn,inverse=True)

    nc = netCDF4.Dataset(nc_outfile,'w',format='NETCDF4_CLASSIC')

    nc.createDimension("x", size=easting.shape[0])
    nc.createDimension("y", size=northing.shape[0])

    var = 'x'
    var_out = nc.createVariable(var, 'float64', dimensions=("x"))
    var_out.units = "meters";
    var_out[:] = easting

    var = 'y'
    var_out = nc.createVariable(var, 'float64', dimensions=("y"))
    var_out.units = "meters";
    var_out[:] = northing

    var = 'lon'
    var_out = nc.createVariable(var, 'float64', dimensions=("y","x"))
    var_out.units = "degreesE";
    var_out.bounds = "lon_bounds"
    var_out[:] = lon

    var = 'lat'
    var_out = nc.createVariable(var, 'float64', dimensions=("y","x"))
    var_out.units = "degreesN";
    var_out.bounds = "lat_bounds"
    var_out[:] = lat

    var = 'ocean_kill_mask'
    var_out = nc.createVariable(var, 'b', dimensions=("y","x"))
    var_out.units = "";
    var_out.long_name = "mask specifying fixed calving front locations"
    var_out.comment = "mask is set via topg and thk"
    var_out.grid_mapping = "mapping"
    var_out.coordinates = "lon lat"
    var_out[:] = 0

    mapping = nc.createVariable("mapping",'c')
    mapping.ellipsoid = "WGS84"
    mapping.false_easting = 0.
    mapping.false_northing = 0.
    mapping.grid_mapping_name = "polar_stereographic"
    mapping.latitude_of_projection_origin = -90.
    mapping.standard_parallel = -71.
    mapping.straight_vertical_longitude_from_pole = 0.

    from time import asctime
    historystr = 'Created ' + asctime() + ' Torsten Albrecht (PIK), using a script by Andy Aschwanden (UAF, USA), modified by Maria Martin (PIK) \n'
    nc.history = historystr
    nc.comments = 'Thanks to Uwe Schulzweida (MPI) for hints on how to create this file'
    nc.Conventions = 'CF-1.4'
    nc.close()

    print "Grid file", nc_outfile, "has been successfully written."

    return nc_outfile


def get_projection_from_file(nc):

    """ This is a copy of the PISM nc2cdo.py code. """

    # First, check if we have a global attribute 'proj4'
    # which contains a Proj4 string:
    try:
        p = Proj(str(nc.proj4))
        print(
            'Found projection information in global attribute proj4, using it')
    except:
        try:
            p = Proj(str(nc.projection))
            print(
                'Found projection information in global attribute projection, using it')
        except:
            try:
                # go through variables and look for 'grid_mapping' attribute
                for var in nc.variables.keys():
                    if hasattr(nc.variables[var], 'grid_mapping'):
                        mappingvarname = nc.variables[var].grid_mapping
                        print(
                            'Found projection information in variable "%s", using it' % mappingvarname)
                        break
                var_mapping = nc.variables[mappingvarname]
                p = Proj(proj="stere",
                         ellps=var_mapping.ellipsoid,
                         datum=var_mapping.ellipsoid,
                         units="m",
                         lat_ts=var_mapping.standard_parallel,
                         lat_0=var_mapping.latitude_of_projection_origin,
                         lon_0=var_mapping.straight_vertical_longitude_from_pole,
                         x_0=var_mapping.false_easting,
                         y_0=var_mapping.false_northing)
            except:
                print('No mapping information found, exiting.')
                sys.exit(1)

    return p


def prepare_ncfile_for_cdo(nc_outfile):

    """ This is a copy of the PISM nc2cdo.py code. """

    # open netCDF file in 'append' mode
    nc = CDF(nc_outfile, 'a')

    # a list of possible x-dimensions names
    xdims = ['x', 'x1']
    # a list of possible y-dimensions names
    ydims = ['y', 'y1']

    # assign x dimension
    for dim in xdims:
        if dim in nc.dimensions.keys():
            xdim = dim
    # assign y dimension
    for dim in ydims:
        if dim in nc.dimensions.keys():
            ydim = dim

    # coordinate variable in x-direction
    x_var = nc.variables[xdim]
    # coordinate variable in y-direction
    y_var = nc.variables[ydim]

    # values of the coordinate variable in x-direction
    easting = x_var[:]
    # values of the coordinate variable in y-direction
    northing = y_var[:]

    # grid spacing in x-direction
    de = np.diff(easting)[0]
    # grid spacing in y-direction
    dn = np.diff(northing)[0]

    # number of grid points in x-direction
    M = easting.shape[0]
    # number of grid points in y-direction
    N = northing.shape[0]

    # number of grid corners
    grid_corners = 4
    # grid corner dimension name
    grid_corner_dim_name = "nv4"

    # array holding x-component of grid corners
    gc_easting = np.zeros((M, grid_corners))
    # array holding y-component of grid corners
    gc_northing = np.zeros((N, grid_corners))
    # array holding the offsets from the cell centers
    # in x-direction (counter-clockwise)
    de_vec = np.array([-de / 2, de / 2, de / 2, -de / 2])
    # array holding the offsets from the cell centers
    # in y-direction (counter-clockwise)
    dn_vec = np.array([-dn / 2, -dn / 2, dn / 2, dn / 2])
    # array holding lat-component of grid corners
    gc_lat = np.zeros((N, M, grid_corners))
    # array holding lon-component of grid corners
    gc_lon = np.zeros((N, M, grid_corners))


    proj = get_projection_from_file(nc)

    # If it does not yet exist, create dimension 'grid_corner_dim_name'
    if grid_corner_dim_name not in nc.dimensions.keys():
        for corner in range(0, grid_corners):
            ## grid_corners in x-direction
            gc_easting[:, corner] = easting + de_vec[corner]
            # grid corners in y-direction
            gc_northing[:, corner] = northing + dn_vec[corner]
            # meshgrid of grid corners in x-y space
            gc_ee, gc_nn = np.meshgrid(
                gc_easting[:, corner], gc_northing[:, corner])
            # project grid corners from x-y to lat-lon space
            gc_lon[:, :, corner], gc_lat[:, :, corner] = proj(
                gc_ee, gc_nn, inverse=True)

        nc.createDimension(grid_corner_dim_name, size=grid_corners)

        var = 'lon_bnds'
        # Create variable 'lon_bnds'
        var_out = nc.createVariable(
            var, 'f', dimensions=(ydim, xdim, grid_corner_dim_name))
        # Assign units to variable 'lon_bnds'
        var_out.units = "degreesE"
        # Assign values to variable 'lon_nds'
        var_out[:] = gc_lon

        var = 'lat_bnds'
        # Create variable 'lat_bnds'
        var_out = nc.createVariable(
            var, 'f', dimensions=(ydim, xdim, grid_corner_dim_name))
        # Assign units to variable 'lat_bnds'
        var_out.units = "degreesN"
        # Assign values to variable 'lat_bnds'
        var_out[:] = gc_lat

    if (not 'lon' in nc.variables.keys()) or (not 'lat' in nc.variables.keys()):
        print("No lat/lon coordinates found, creating them")
        ee, nn = np.meshgrid(easting, northing)
        lon, lat = proj(ee, nn, inverse=True)

    var = 'lon'
    # If it does not yet exist, create variable 'lon'
    if not var in nc.variables.keys():
        var_out = nc.createVariable(var, 'f', dimensions=(ydim, xdim))
        # Assign values to variable 'lon'
        var_out[:] = lon
    else:
        var_out = nc.variables[var]
    # Assign units to variable 'lon'
    var_out.units = "degreesE"
    # Assign long name to variable 'lon'
    var_out.long_name = "Longitude"
    # Assign standard name to variable 'lon'
    var_out.standard_name = "longitude"
    # Assign bounds to variable 'lon'
    var_out.bounds = "lon_bnds"

    var = 'lat'
    # If it does not yet exist, create variable 'lat'
    if not var in nc.variables.keys():
        var_out = nc.createVariable(var, 'f', dimensions=(ydim, xdim))
        var_out[:] = lat
    else:
        var_out = nc.variables[var]
    # Assign units to variable 'lat'
    var_out.units = "degreesN"
    # Assign long name to variable 'lat'
    var_out.long_name = "Latitude"
    # Assign standard name to variable 'lat'
    var_out.standard_name = "latitude"
    # Assign bounds to variable 'lat'
    var_out.bounds = "lat_bnds"

    # Make sure variables have 'coordinates' attribute
    for var in nc.variables.keys():
        if (nc.variables[var].ndim >= 2):
            nc.variables[var].coordinates = "lon lat"

    # lat/lon coordinates must not have mapping and coordinate attributes
    # if they exist, delete them
    for var in ['lat', 'lon', 'lat_bnds', 'lon_bnds']:
        if hasattr(nc.variables[var], 'grid_mapping'):
            delattr(nc.variables[var], 'grid_mapping')
        if hasattr(nc.variables[var], 'coordinates'):
            delattr(nc.variables[var], 'coordinates')

    # If present prepend history history attribute, otherwise create it
    from time import asctime
    histstr = asctime() + \
        ' : grid info for CDO added by nc2cdo.py, a PISM utility\n'
    if 'History' in nc.ncattrs():
        nc.History = histstr + nc.History
    elif 'history' in nc.ncattrs():
        nc.history = histstr + nc.history
    else:
        nc.history = histstr

    for attr in ("projection", "proj4"):
        if hasattr(nc, attr):
            delattr(nc, attr)
    # Write projection attribute
    nc.proj4 = proj.srs
    # Close file
    nc.close()

    print "Prepared file",nc_outfile,"for cdo."
