
"""
Create grids based on the ALBMAP dimensions, that are compatible with PISM.
That means, if they are fed into PISM, PISM will not need to regrid internally.
The here created grids are used as input for the cdo remapping, which is done
for each dataset separately. See for example bedmap2/remap.py

"""

import os, sys
import numpy as np
from pyproj import Proj
import netCDF4

## this hack is needed to import config.py from the project root
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if project_root not in sys.path: sys.path.append(project_root)
import config as cf; reload(cf)
import pism_input.pism_input as pi; reload(pi)

# if true, PISM will not need to regrid
use_PISM_grid=True

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


if __name__ == "__main__":

    for resolution in [5,10,20,30,50]:

        # cdo_targetgrid_file = os.path.join(cf.cdo_remapgridpath,
        #                         'pism_'+str(int(resolution))+'km.nc')
        gridfile = create_grid_for_cdo_remap(cf.cdo_remapgridpath,
                        use_PISM_grid=True, resolution=resolution)
        pi.prepare_ncfile_for_cdo(gridfile)
