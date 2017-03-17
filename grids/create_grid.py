#!/usr/bin/env python
import numpy as np
from pyproj import Proj
from optparse import OptionParser

# try different netCDF modules
try:
    from netCDF4 import Dataset as CDF
except:
    from netCDF3 import Dataset as CDF

# default values
DXY = 15.0  # km

use_PISM_grid = True
PISM_grid_set={}
#PISM_grid_set[15]=[Mx,My,xl,yl,xr,yr]
PISM_grid_set[15]=[400,400,-2795000,-2795000,3190000,3190000]
PISM_grid_set[12]=[500,500,-2796500,-2796500,3191500,3191500]



# set up the option parser
parser = OptionParser()
parser.usage = "usage: %prog [options] FILE"
parser.description = "Create CDO-compliant grid description"
parser.add_option("-g", "--grid_spacing",dest="grid_spacing",type='float',
                  help="use X km grid spacing",
                  metavar="X",default=DXY)

(options, args) = parser.parse_args()
grid_spacing = options.grid_spacing*1e3 # convert km -> m

if len(args) == 0:
    nc_outfile = 'grid_' + str(grid_spacing/1e3) + 'km.nc'
elif len(args) == 1:
    nc_outfile = args[0]
else:
    print('wrong number arguments, 0 or 1 arguments accepted')
    parser.print_help()
    exit(0)


if __name__ == "__main__": 

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

      e0 = -2800000.0
      n0 = -2800000.0
      e1 =  3195000.0
      n1 =  3195000.0
      #x: -2800000 to 3195000
      #y: -2800000 to 3195000

      M = int((e1 - e0)/de) + 1
      N = int((n1 - n0)/dn) + 1

      easting  = np.linspace(e0, e1, M)
      northing = np.linspace(n0, n1, N)


    ee, nn = np.meshgrid(easting,northing)

    print "Grid is created for "+str(int(grid_spacing/1e3))+"km resolution:"
    print M,N,easting[0],northing[0],easting[-1],northing[-1],np.diff(easting)[0],np.diff(northing)[0]


    # Set up SeaRISE Projection
    projection = "+proj=stere +ellps=WGS84 +datum=WGS84 +lon_0=0 +lat_0=-90 +lat_ts=-71 +units=m"

    proj = Proj(projection)

    lon,lat = proj(ee,nn,inverse=True)

    nc = CDF(nc_outfile,'w',format='NETCDF4_CLASSIC')

    nc.createDimension("x", size=easting.shape[0])
    nc.createDimension("y", size=northing.shape[0])

    var = 'x'
    var_out = nc.createVariable(var, 'f', dimensions=("x"))
    var_out.units = "meters";
    var_out[:] = easting

    var = 'y'
    var_out = nc.createVariable(var, 'f', dimensions=("y"))
    var_out.units = "meters";
    var_out[:] = northing

    var = 'lon'
    var_out = nc.createVariable(var, 'f', dimensions=("y","x"))
    var_out.units = "degreesE";
    var_out.bounds = "lon_bounds"
    var_out[:] = lon

    var = 'lat'
    var_out = nc.createVariable(var, 'f', dimensions=("y","x"))
    var_out.units = "degreesN";
    var_out.bounds = "lat_bounds"
    var_out[:] = lat
      
    var = 'dummy'
    var_out = nc.createVariable(var, 'f', dimensions=("y","x"))
    var_out.units = "meters";
#    var_out._FillValue = np.nan
    var_out.long_name = "Just A Dummy"
    var_out.comment = "This is just a dummy variable for CDO."
    var_out.grid_mapping = "mapping"
    var_out.coordinates = "lon lat"    
    var_out[:] = np.nan

    mapping = nc.createVariable("mapping",'c')
    mapping.ellipsoid = "WGS84"
    mapping.false_easting = 0.
    mapping.false_northing = 0.
    mapping.grid_mapping_name = "polar_stereographic"
    mapping.latitude_of_projection_origin = -90.
    mapping.standard_parallel = -71.
    mapping.straight_vertical_longitude_from_pole = 0.

    from time import asctime
    historystr = 'Created ' + asctime() + ' by Maria Martin (PIK), using a script by Andy Aschwanden, University of Alaska Fairbanks, USA \n'
    nc.history = historystr
    nc.comments = 'Thanks to Uwe Schulzweida for hints on how to create this file'
    nc.Conventions = 'CF-1.4'
    nc.close()
