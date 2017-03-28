"""
Convert conservative to potential temperatures and absolute to practical salinity in Schmidtko Data
These functions are based on the gsw package, 
for salinity conversion, we use one reference SAAR value (computed in calculate_reference_saar.py)
as this value might not exits close to the continent
ronja.reese@pik-potsdam.de
"""

import os
import numpy as np
import numpy.ma as ma
from shutil import copyfile
import netCDF4 as nc

from gsw_functions import pt_from_CT, SP_from_SA_Antarctica

datafile ='schmidtko_data/schmidtko_ocean_input.nc'
copyfile('schmidtko_data/schmidtko_ocean.nc', datafile)

print 'Reading ', datafile
infile 		= nc.Dataset(datafile, 'r+')
temperature = ma.masked_array(infile.variables['thetao'][:]) #time, lat, lon
salinity 	= ma.masked_array(infile.variables['salinity'][:]) #time,lat,lon
height 	 	= ma.masked_array(infile.variables['height'][:]) #time,lat,lon
lat 	 	= ma.masked_array(infile.variables['lat'][:]) # lat, degrees north 
lon 	 	= ma.masked_array(infile.variables['lon'][:]) # lon, degrees east



compute_potential_temperatures = True
compute_practical_salinity 	   = True
temp_in_kelvin=False


if (compute_potential_temperatures):

	potential_temps 	= np.zeros_like(temperature)	#time,depth,lat,lon
	potential_temps[0,:,:] = pt_from_CT(salinity[0,:,:], temperature[0,:,:]) # needs temperature in degC
	if(temp_in_kelvin):
		potential_temps = potential_temps + 273.15



if (compute_practical_salinity) :
	practical_salinity 	= np.zeros_like(salinity)	#time,depth,lat,lon
	practical_salinity[0,:,:] = SP_from_SA_Antarctica(salinity[0,:,:])


#save converted values in file
infile.variables['thetao'][:] = potential_temps
infile.variables['salinity'][:] = practical_salinity
infile.variables['salinity'].units = '' # now dimensionless
infile.close()

print 'Done'


# # plot the difference between potential and absolute salinity 
# import matplotlib.pyplot as plt

# fig=plt.figure()
# ax = plt.subplot(1,1,1) 

# plt.contourf(practical_salinity[0,:,:]-salinity[0,:,:]) 

# ax.set_xticks([]) 
# ax.set_yticks([]) 

# #plt.show()
# plt.savefig(os.path.join(schmidtko_data_path,"plots/Practical_versus_absolute_salinity.png"))
# plt.clf()

