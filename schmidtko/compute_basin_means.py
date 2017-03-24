#!/usr/bin/env python

import numpy as np
import numpy.ma as ma
import netCDF4 as nc
from shutil import copyfile


infile='schmidtko_data/schmidtko_15km.nc'
basinfile='../basins/basins_data/basins_zwally12_15km.nc'
outfile='schmidtko_data/schmidtko_15km_means.nc'

# Load data
readfile = nc.Dataset(infile, 'r')
temperature = readfile.variables['thetao'][:]
salinity = readfile.variables['salinity'][:]
x = readfile.variables['x'][:]
y = readfile.variables['y'][:]
readfile.close()

#load basin data
basinfile = nc.Dataset(basinfile, 'r')
basins = basinfile.variables['basins'][:]
basinfile.close()


# compute mean per basin
temperature_means = np.zeros(basins.max()+1)
salinity_means = np.zeros(basins.max()+1)
tempcount = np.zeros(basins.max()+1)
salcount = np.zeros(basins.max()+1)

for i in range(len(x)):
	for j in range(len(y)):
		#basin_id = int(basins[0,i,j])
                basin_id = int(basins[0,i,j])
		if np.logical_not(np.isnan(temperature[0,i,j])):
			temperature_means[basin_id] = temperature_means[basin_id] + temperature[0,i,j]
			tempcount[basin_id] = tempcount[basin_id] +1 
		if np.logical_not(np.isnan(salinity[0,i,j])):
			salinity_means[basin_id] = salinity_means[basin_id] + salinity[0,i,j]
			salcount[basin_id] = salcount[basin_id] +1 

for k in range(int(basins.max())+1):
	temperature_means[k] = temperature_means[k]*1.0/(1.0*tempcount[k])
	salinity_means[k]    = salinity_means[k]*1.0/(1.0*salcount[k])

# Value for basin 1 is missing: setting it te mean between basin 2 and 17:
#temperature_means[1] = 0.5*(temperature_means[2] + temperature_means[17])
#salinity_means[1] = 0.5*(salinity_means[2] + salinity_means[17])
temperature_means[4] = 0.5*(temperature_means[3] + temperature_means[5])
salinity_means[4] = 0.5*(salinity_means[3] + salinity_means[5])

print temperature_means
print salinity_means

for i in range(len(x)):
	for j in range(len(y)):
                #basin_id = int(basins[i,j])
		basin_id = int(basins[0,i,j])
		temperature[0,i,j] = temperature_means[basin_id]
		salinity[0,i,j] = salinity_means[basin_id]



#save in new file

copyfile(infile, outfile)

wrtfile = nc.Dataset(outfile, 'r+')

wrtfile.variables['thetao'][:] = temperature
wrtfile.variables['salinity'][:] = salinity

wrtfile.close()
