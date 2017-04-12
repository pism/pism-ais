"""
Compute average schmidtko input for each basin and fill this into the basins values
Missing values are set to averages between the adjacent basins

add basins variable to the output file

ronja.reese@pik-potsdam.de

FIXME: make the resolution chosable...
"""

import os, sys
import numpy as np
import numpy.ma as ma
import netCDF4 as nc
from shutil import copyfile

## this hack is needed to import config.py from the project root
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if project_root not in sys.path: sys.path.append(project_root)
import config as cf; reload(cf)
import pism_input.pism_input as pi; reload(pi)

# save data path
dataset="schmidtko"
# resolution for the output file, set from config.
resolution = cf.resolution # in km

data_path = os.path.join(cf.output_data_path, dataset)
infile = os.path.join(data_path, 'schmidtko_'+str(resolution)+'km.nc')
outfile = os.path.join(data_path, 'schmidtko_'+str(resolution)+'km_means.nc')
print infile

zwally_data_path = os.path.join(cf.output_data_path,"zwally_basins")
basinfile = os.path.join(zwally_data_path, 'zwally_basins_'+str(resolution)+'km.nc')


# Load data
readfile = nc.Dataset(infile, 'r')
temperature = readfile.variables['theta_ocean'][:]
salinity = readfile.variables['salinity_ocean'][:]
x = readfile.variables['x'][:]
y = readfile.variables['y'][:]
readfile.close()

#load basin data
basinfile = nc.Dataset(basinfile, 'r')
basins = basinfile.variables['basins'][:]
basinfile.close()


# compute mean per basin
temperature_means = np.zeros(int(basins.max())+1)
salinity_means = np.zeros(int(basins.max())+1)
tempcount = np.zeros(int(basins.max())+1)
salcount = np.zeros(int(basins.max())+1)

for i in range(len(x)):
  for j in range(len(y)):
    basin_id = int(basins[i,j])
    if np.logical_not(np.isnan(temperature[0,i,j])):
      temperature_means[basin_id] = temperature_means[basin_id] + temperature[0,i,j]
      tempcount[basin_id] += 1
    if np.logical_not(np.isnan(salinity[0,i,j])):
      salinity_means[basin_id] = salinity_means[basin_id] + salinity[0,i,j]
      salcount[basin_id] += 1


for k in range(int(basins.max())+1):
  if tempcount[k]!=0.0:
    temperature_means[k] = temperature_means[k]*1.0/(1.0*tempcount[k])
  if salcount[k]!=0.0:
    salinity_means[k]    = salinity_means[k]*1.0/(1.0*salcount[k])

# Value for basin 1 is missing: setting it te mean between basin 2 and 17:
#temperature_means[1] = 0.5*(temperature_means[2] + temperature_means[17])
#salinity_means[1] = 0.5*(salinity_means[2] + salinity_means[17])
temperature_means[4] = 0.5*(temperature_means[3] + temperature_means[5])
salinity_means[4] = 0.5*(salinity_means[3] + salinity_means[5])

for i in range(len(x)):
  for j in range(len(y)):
    basin_id = int(basins[i,j])
    temperature[0,i,j] = temperature_means[basin_id]
    salinity[0,i,j] = salinity_means[basin_id]



#save in new file

copyfile(infile, outfile)

wrtfile = nc.Dataset(outfile, 'r+')

wrtfile.variables['theta_ocean'][:] = temperature
wrtfile.variables['salinity_ocean'][:] = salinity
wrtfile.variables['theta_ocean'].units = "Celsius"
wrtfile.variables['salinity_ocean'].units = "g/kg"

## add variable basins to the file
basinsvar = wrtfile.createVariable('basins', 'f8', ('time','y', 'x'))
basins_variable = np.zeros_like(temperature)
basins_variable[0,:,:] = basins
basinsvar[:] = basins_variable
basinsvar.units = ""
basinsvar.long_name = "drainage basins"
basinsvar.standard_name = "drainage_basins"

wrtfile.close()
