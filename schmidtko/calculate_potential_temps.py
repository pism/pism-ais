"""
Convert conservative to potential temperatures and absolute to practical salinity in Schmidtko Data
These functions are based on the gsw package,
for salinity conversion, we use one reference SAAR value (computed in calculate_reference_saar.py)
as this value might not exits close to the continent
ronja.reese@pik-potsdam.de
"""

import os, sys
import numpy as np
import numpy.ma as ma
from shutil import copyfile
import netCDF4 as nc

from gsw_functions import pt_from_CT, SP_from_SA_Antarctica

## this hack is needed to import config.py from the project root
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if project_root not in sys.path: sys.path.append(project_root)
import config as cf; reload(cf)

# save data path
dataset="schmidtko"

compute_potential_temperatures = True
compute_practical_salinity     = True
temp_in_kelvin=False


data_path = os.path.join(cf.output_data_path, dataset)
datafile = os.path.join(data_path, 'schmidtko_ocean_input_potentialtemps.nc')
copyfile(os.path.join(data_path, 'schmidtko_ocean_input.nc') , datafile)

print 'Reading ', datafile
infile 		= nc.Dataset(datafile, 'r+')
temperature = ma.masked_array(infile.variables['theta_ocean'][:]) #time, lat, lon
salinity 	= ma.masked_array(infile.variables['salinity_ocean'][:]) #time,lat,lon
height 	 	= ma.masked_array(infile.variables['height'][:]) #time,lat,lon
lat 	 	= ma.masked_array(infile.variables['lat'][:]) # lat, degrees north
lon 	 	= ma.masked_array(infile.variables['lon'][:]) # lon, degrees east




if (compute_potential_temperatures):

	potential_temps 	= np.zeros_like(temperature)	#time,depth,lat,lon
	potential_temps[0,:,:] = pt_from_CT(salinity[0,:,:], temperature[0,:,:]) # needs temperature in degC
	if(temp_in_kelvin):
		potential_temps = potential_temps + 273.15
	potential_temps = potential_temps.filled(np.nan) # set masked values to nan



if (compute_practical_salinity) :
	practical_salinity 	= np.zeros_like(salinity)	#time,depth,lat,lon
	practical_salinity[0,:,:] = SP_from_SA_Antarctica(salinity[0,:,:])
	practical_salinity = practical_salinity.filled(np.nan)


#save converted values in file
infile.variables['theta_ocean'][:] = potential_temps # FIXME this changes fill values from NaNf to -- which makes a problem in cdo...!!!
infile.variables['salinity_ocean'][:] = practical_salinity
infile.variables['salinity_ocean'].units = '' # now dimensionless
infile.close()

print "Temperature and salinity successfully converted and written to", datafile


# # plot the difference between potential and absolute salinity_ocean
# import matplotlib.pyplot as plt

# fig=plt.figure()
# ax = plt.subplot(1,1,1)

# plt.contourf(practical_salinity[0,:,:]-salinity_ocean[0,:,:])

# ax.set_xticks([])
# ax.set_yticks([])

# #plt.show()
# plt.savefig(os.path.join(schmidtko_data_path,"plots/Practical_versus_absolute_salinity_ocean.png"))
# plt.clf()

