"""
matthias.mengel@pik, torsten.albrecht@pik
Download Albmap data and save rename variables to make the PISM compatible.
ALBMAP is documented here:
http://www.earth-syst-sci-data.net/2/247/2010/
"""

import os, sys
import subprocess
import numpy as np
import netCDF4 as nc
## this hack is needed to import config.py from the project root
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if project_root not in sys.path: sys.path.append(project_root)
import config as cf; reload(cf)
import pism_input.pism_input as pi; reload(pi)

dataset = "albmap"
data_path = os.path.join(cf.output_data_path, dataset)

inputfile = os.path.join(data_path, 'albmap_5km_input.nc')

cdo_targetgrid_file, regridded_file = pi.get_filenames_for_cdo(
    cf.cdo_remapgridpath, data_path, dataset, cf.grid_id)

nci = nc.Dataset(inputfile,"r")
ncr = nc.Dataset(regridded_file,"a")

for var in ["x","y","lon","lat","air_temp","bheatflx","precipitation","thk","topg"]:

    ncr.variables[var].units = nci.variables[var].units


for var in ["air_temp","bheatflx","precipitation","thk","topg"]:

    variable = ncr.variables[var][:]
    mask = variable.mask.copy()
    variable = np.array(variable)
    variable[mask] = 0
    ncr.variables[var][:] = variable

nci.close()
ncr.close()

print "Added back units and filled missing values in bound regions with zero in"
print regridded_file
