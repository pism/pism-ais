#!/usr/bin/env python

import gsw
import numpy as np
import numpy.ma as ma
from shutil import copyfile

#Compute potential temperature from conservative temperature in Sunke Data

# try different netCDF modules
try:
    import netCDF4 as nc
except:
    import netCDF3 as nc

datafile ='schmidtko_data/schmidtko_ocean_input.nc'
copyfile('schmidtko_data/schmidtko_ocean.nc', datafile)


print '\nReading data... '

print 'Read ', datafile
infile = nc.Dataset(datafile, 'r+')
temperature = ma.masked_array(infile.variables['thetao'][:]) #time, lat, lon
salinity = ma.masked_array(infile.variables['salinity'][:]) #time,lat,lon



#COMPUTE POTENTIAL TEMPERATURES

compute_potential_temperatures = True
temp_in_kelvin=False


if (compute_potential_temperatures):

	potential_temps 	= np.zeros_like(temperature)	#time,depth,lat,lon
	
	potential_temps[0,:,:] = gsw.pt_from_CT(salinity[0,:,:], temperature[0,:,:]) # needs temperature in degreeCelsius

	if(temp_in_kelvin):
		potential_temps = potential_temps + 273.15
	else:
		potential_temps = potential_temps


#save temperature in file
infile.variables['thetao'][:] = potential_temps

infile.close()

print 'Done'





