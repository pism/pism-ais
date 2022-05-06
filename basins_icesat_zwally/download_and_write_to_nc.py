"""
matthias.mengel@pik, torsten.albrecht@pik

This is preprocessing for the ICESat (Zwally) drainage basin data, accessible from
https://earth.gsfc.nasa.gov/cryo/data/polar-altimetry/antarctic-and-greenland-drainage-systems
This script downloads the data and saves it to a netcdf file.

The ICESat drainage basins are available on 1km and updated from the Zwally et al. 2012

"""

import os, glob
import numpy as np
import numpy.ma as ma
import sys, csv, datetime
import netCDF4 as nc
import datetime, math
import pandas as pd

## this hack is needed to import config.py from the project root
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if project_root not in sys.path: sys.path.append(project_root)
import config as cf; reload(cf)
import pism_input.pism_input as pi; reload(pi)

dataset="basins_icesat_zwally"
data_link="https://earth.gsfc.nasa.gov/sites/default/files/lab_cryo/data/polar_ice_altimetry/antarctic_and_greenland_drainage_systems/ant_icesat_modis_mask_1km_ascii.tgz"

# exclude floating ice and icy islands from drainage basins
exclude_shelves_and_islands = True

basins_data_path = os.path.join(cf.output_data_path, dataset)
basins_ascii_file = os.path.join(basins_data_path, 'Ant_ICESat_MODIS_Mask_1km.ascii')
ncout_file = os.path.join(basins_data_path,dataset+"_1km_input.nc")

# if data is not yet extracted in basins_ascii
if not os.path.exists(basins_ascii_file):
  print "Downloading ICESat (Zwally) basin ascii data."
  os.system("mkdir -pv " + basins_data_path)
  os.system("wget -N " + data_link + " -P " + basins_data_path)
  os.system("cd "+basins_data_path+" && tar -xvf ant_icesat_modis_mask_1km_ascii.tgz")


raw_data = pd.read_csv(basins_ascii_file, skiprows=40, header=None, delimiter=r"\s+",
                 names=["x", "y", "Lat", "Lon", "DS ID", "Surface code"],)

# from Ant_ICESat_MODIS_Mask_1km.ascii header
x = np.arange(1950,7311,1)
y = np.arange(2175,6751,1)
X,Y = np.meshgrid(x, y)

latgrid = np.zeros_like(X, dtype=float)
longrid = np.zeros_like(X, dtype=float)
drainage_id_grid = np.zeros_like(X)
surface_code_grid = np.zeros_like(X)

# for performance
xvals = raw_data.loc[:,"x"].values
yvals = raw_data.loc[:,"y"].values
lats = raw_data.loc[:,"Lat"].values
lons = raw_data.loc[:,"Lon"].values
ds_id = raw_data.loc[:,"DS ID"].values
surface_code = raw_data.loc[:,"Surface code"].values

print "looping over all ascii entries, 13.5 million."
for i in raw_data.index[::1]:

    if i%100000 == 0:
        print i,

    xi = np.searchsorted(x,xvals[i])
    yi = np.searchsorted(y,yvals[i])
    latgrid[yi,xi] = lats[i]
    longrid[yi,xi] = lons[i]
    drainage_id_grid[yi,xi] = ds_id[i]
    surface_code_grid[yi,xi] = surface_code[i]

if exclude_shelves_and_islands:

    shelves_or_islands = surface_code_grid > 4
    drainage_id_grid[shelves_or_islands] = 0

ncout = nc.Dataset(ncout_file, 'w')
#ncout.createDimension('time',size=None)
ncout.createDimension('x',size=len(x))
ncout.createDimension('y',size=len(y))
ncx = ncout.createVariable(varname="x",datatype='float_',dimensions=('x'))
ncy = ncout.createVariable(varname="y",datatype='float_',dimensions=('y'))
ncx[:] = x*1.e3
ncy[:] = y*1.e3
ncx.units = "m"
ncy.units = "m"

ncv = ncout.createVariable(varname="lat",datatype='float_',dimensions=('y','x') )
ncv[:] = np.flipud(latgrid)
ncv = ncout.createVariable(varname="lon",datatype='float_',dimensions=('y','x') )
ncv[:] = np.flipud(longrid)
ncv = ncout.createVariable(varname="mask",datatype='int_',dimensions=('y','x') )
ncv[:] = np.flipud(surface_code_grid)
ncv = ncout.createVariable(varname="basins",datatype='int_',dimensions=('y','x') )
ncv[:] = np.flipud(drainage_id_grid)

now = datetime.datetime.now().strftime("%B %d, %Y")

ncout.data_origin = "https://earth.gsfc.nasa.gov/cryo/data/polar-altimetry/antarctic-and-greenland-drainage-systems"
ncout.proj4 = cf.proj4str
ncout.comment  = cf.authors+" created netcdf basins file at " + now
ncout.close()

# prepare the input file for cdo remapping
# this step takes a while for high resolution data (i.e. 1km)
pi.prepare_ncfile_for_cdo(ncout_file)

print "Preprocessed and wrote netcdf file", ncout_file
