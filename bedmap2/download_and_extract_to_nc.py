"""
matthias.mengel@pik, torsten.albrecht@pik
Download Bedmap 2 data and save to (1km) netcdf file.
"""

import os
import numpy as np
import sys
import netCDF4
import datetime

## this hack is needed to import config.py from the project root
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if project_root not in sys.path: sys.path.append(project_root)
import config as cf; reload(cf)

### Bedmap2   ##########################################################
# Documentation of the data: https://www.bas.ac.uk/project/bedmap-2/
bedmap2_link="https://secure.antarctica.ac.uk/data/bedmap2/bedmap2_bin.zip"
bedmap2_data_path = os.path.join(cf.output_data_path, "bedmap2")
ncout_name = os.path.join(bedmap2_data_path, 'bedmap2_1km.nc')

# if data is not yet extracted in bedmap2_bin
if not os.path.exists(os.path.join(bedmap2_data_path,"bedmap2_bin")):
  print "Downloading bedmap2 binary data."
  os.system("mkdir " + bedmap2_data_path)
  os.system("wget -N " + bedmap2_link + " -P " + bedmap2_data_path)
  os.system("cd "+bedmap2_data_path+" && unzip bedmap2_bin.zip")

# if os.path.isfile(ncout_name):
#   print ncout_name, "already written, do nothing."

data_files = {"topg":"bedmap2_bed.flt",
              "thk":"bedmap2_thickness.flt",
              "mask":'bedmap2_icemask_grounded_and_shelves.flt',
              "bedunc":"bedmap2_grounded_bed_uncertainty.flt",
              "usurf":"bedmap2_surface.flt"}

# taken from bedmap2 readme
N=6667
dx = 1000.0 #m
dy = 1000.0 #m
x = np.linspace(0.0,(N-1)*dx,N)
y = np.linspace(0.0,(N-1)*dy,N)

print "Reading bedmap2 binary files from %s ...\n" % (bedmap2_data_path)

bedm2_vars = {}
for var, file in data_files.iteritems():
  fname = os.path.join(bedmap2_data_path,"bedmap2_bin",file)
  vardata = np.flipud(np.ma.masked_equal(np.reshape(
          np.fromfile(fname,dtype=np.float32),(N,N)),-9999.0))
  print " range of "+str(var)+" = [%.2f, %.2f]" % (vardata.min(),vardata.max())
  bedm2_vars[var] =  vardata


bedm2_attributes = {"bed": {"long_name" : "elevation of bedrock",
                           "valid_range" : (-9000.0, 9000.0),
                           "standard_name" : "bedrock_altitude",
                           "units" : "meters"},
                    "usurf": {"long_name" : "ice upper surface elevation",
                            "valid_range" : (-1000.0, 9000.0),
                            "standard_name" : "surface_altitude",
                            "units" : "meters"},
                    "thk": {"long_name" : "thickness of ice sheet or ice shelf",
                            "valid_range" : (0.0, 9000.0),
                            "standard_name" : "land_ice_thickness",
                            "units" : "meters"},
                    "bedunc": {"long_name" : "uncertainty of bed topography",
                            "valid_range" : (0.0, 9000.0),
                            "standard_name" : "bed_uncertainty",
                            "units" : "meters"},
                    "mask": {"long_name" : "ice-type (ice-free/grounded/floating/ocean) integer mask",
                            "valid_range" : (0.0, 1.0),
                            "standard_name" : "mask",
                            "units" : ""} }

### Write nc-file

print "Writing NetCDF file '%s' ..." % ncout_name
ncout = netCDF4.Dataset(ncout_name,"w",format='NETCDF4_CLASSIC')

# no time dimension needed here
ncout.createDimension('x',size=len(x))
ncout.createDimension('y',size=len(x))
ncx   = ncout.createVariable( 'x','float32',('x',) )
ncy   = ncout.createVariable( 'y','float32',('y',) )
ncx[:] = x
ncy[:] = y

for varname,data in bedm2_vars.iteritems():

  ncvar = ncout.createVariable( varname,'float32',('y','x') )
  ncvar[:] = data
  for att in bedm2_attributes[var]:
    pass

now = datetime.datetime.now().strftime("%B %d, %Y")
ncout.proj4      = "+proj=stere +lon_0=0 +lat_0=-90 +lat_ts=-71 +ellps=WGS84 +datum=WGS84"
ncout.projection = "+proj=stere +lon_0=0 +lat_0=-90 +lat_ts=-71 +ellps=WGS84 +datum=WGS84"
ncout.comment  = cf.authors+" created netcdf bedmap2 file at " + now

ncout.close()
print "Done"

