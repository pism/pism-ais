#!/usr/bin/env python

import numpy as np
import numpy.ma as ma
import netCDF4 as nc

read_data = True
if read_data:
	# Read data from txt file:
	lon = np.zeros(0)
	lat = np.zeros(0)
	depth = np.zeros(0)
	c_t = np.zeros(0)
	c_t_std = np.zeros(0)
	abs_s = np.zeros(0)
	abs_s_std = np.zeros(0)


	f = open('schmidtko_data/Antarctic_shelf_data.txt', 'r')

	for i,line in enumerate(f):
          if i>0:    
		#print(i)
		line = line.strip()
		columns = line.split()
		lon=np.append(lon, float(columns[0]))
		lat=np.append(lat, float(columns[1]))
		depth=np.append(depth, float(columns[2]))
		c_t=np.append(c_t, float(columns[3]))
		c_t_std=np.append(c_t_std, float(columns[4]))
		abs_s=np.append(abs_s, float(columns[5]))
		abs_s_std=np.append(abs_s_std, float(columns[6]))

	f.close()

# Convert to pism grid...
lat_new = np.sort(np.unique(lat)) #are sorted already
lon_new = np.sort(np.unique(lon))

# Move longitude vals to fit 
lon_new = lon_new -360
lon_new[lon_new<-180] = 360 + lon_new[lon_new<-180]
lon_new = np.sort(lon_new)

#create array for these dimensions and fill in values:
fillvalue = np.nan

thetao = np.zeros((1,len(lat_new), len(lon_new))) + fillvalue
salinity = np.zeros((1,len(lat_new), len(lon_new))) + fillvalue

for i in range(len(c_t)):
	# go through all temp and sal vals and fill them into the right place
	ilat = np.in1d(lat_new.ravel(), lat[i]).reshape(lat_new.shape)
	comp_lon = lon[i] -360
	if comp_lon < -180:
		comp_lon = 360 + comp_lon
	ilon = np.in1d(lon_new.ravel(), comp_lon).reshape(lon_new.shape)
	thetao[0,ilat,ilon] = c_t[i]
	salinity[0,ilat,ilon] = abs_s[i]

# Save data as NetCDF file
save_in_file = True
temp_in_kelvin = False


if(save_in_file):
	print 'save data to file....'
	wrtfile = nc.Dataset('schmidtko_data/schmidtko_ocean.nc', 'w', format='NETCDF4_CLASSIC')
   	wrtfile.createDimension('lon', size=len(lon_new))
    	wrtfile.createDimension('lat', size=len(lat_new))
	wrtfile.createDimension('time', size=None)
	wrtfile.createDimension('nv', size=2)

    	nct = wrtfile.createVariable('time', 'float32', ('time',))
	nctb = wrtfile.createVariable('time_bnds', 'float32', ('time','nv'))
    	nclat = wrtfile.createVariable('lat', 'f4', ('lat',))
	nclon = wrtfile.createVariable('lon', 'f4', ('lon',))
	nctemp = wrtfile.createVariable('thetao', 'f4', ('time','lat', 'lon'))
	ncsal = wrtfile.createVariable('salinity', 'f4', ('time','lat', 'lon'))
	
	tm = np.arange(1,1.5,1)
	nct[:] = tm
	nct.units = 'years since 01-01-01 00:00:00'
	nct.bounds = 'time_bnds'

	t_bds = np.zeros([len(tm),2])
	t_bds[:,0] = tm -5
	t_bds[0,:] = tm + 5
	nctb[:] = t_bds
	nctb.units = nct.units

	nclon[:] = lon_new[:]
	nclon.units = 'degrees_east'
	nclat[:] = lat_new[:]
	nclat.units = 'degrees_north'

	ncsal[:] = salinity
	ncsal.units = ''
	nctemp[:] = thetao
	if(temp_in_kelvin):
		nctemp.units='Kelvin'
	else:
		nctemp.units='degree_Celsius'

	wrtfile.close()

print 'Done'
