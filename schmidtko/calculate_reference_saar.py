"""
Calculate reference saar (salinity anomaly ratio, used to compute practical from absolute) value 
Not necessary for data preprocessing, since value is already inserted in the gsw_functions file
ronja.reese@pik-potsdam.de
"""


import numpy as np
import numpy.ma as ma
from shutil import copyfile
import netCDF4 as nc

import gsw

datafile ='schmidtko_data/schmidtko_ocean.nc'

print 'Reading ', datafile
infile 		= nc.Dataset(datafile, 'r')
temperature = ma.masked_array(infile.variables['thetao'][:]) #time, lat, lon
salinity 	= ma.masked_array(infile.variables['salinity'][:]) #time,lat,lon
height 	 	= ma.masked_array(infile.variables['height'][:]) #time,lat,lon
lat 	 	= ma.masked_array(infile.variables['lat'][:]) # lat, degrees north 
lon 	 	= ma.masked_array(infile.variables['lon'][:]) # lon, degrees east




practical_salinity 	= np.zeros_like(salinity)	#time,depth,lat,lon
pressure		 	= np.zeros_like(height)		#time,depth,lat,lon
saar			 	= np.zeros_like(salinity)		#time,depth,lat,lon
in_ocean			= np.zeros_like(salinity)		#time,depth,lat,lon

for i,lati in enumerate(lat):
	pressure[0,i,:] = gsw.p_from_z(height[0,i,:], lat[i])  # sea pressure in dbar
	for j,lonj in enumerate(lon):

		# using the code below from gsw.SP_from_SA since this was broken (-baltic)
		#practical_salinity[0,i,j] = gsw.SP_from_SA(salinity[0,i,j], pressure[0,i,j] , lon[j], lat[i]) 

	    saar[0,i,j], in_ocean[0,i,j] = gsw.library.SAAR(pressure[0,i,j], lonj, lati)
	    if in_ocean[0,i,j]!=False:
			print saar[0,i,j], in_ocean[0,i,j]

			practical_salinity[0,i,j] = (35.0 / 35.16504) * salinity[0,i,j] / (1.0 + saar[0,i,j])

saar_for_Antarctica = saar[in_ocean==True].mean()

    		

print 'Saar for Antarctica', saar_for_Antarctica





