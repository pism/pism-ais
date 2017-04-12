"""
Download Southern Ocean data from Schmidtko et al, Science 2015
http://science.sciencemag.org/content/346/6214/1227
and create NetCDF file from the txt data.
ronja.reese@pik-potsdam.de, matthias.mengel@pik

FIXME use correct pism names for variables
"""

import os, sys
import numpy as np
import numpy.ma as ma
import netCDF4 as nc

## this hack is needed to import config.py from the project root
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if project_root not in sys.path: sys.path.append(project_root)
import config as cf; reload(cf)

# save data path
dataset="schmidtko"

temp_in_kelvin = False

data_path = os.path.join(cf.output_data_path, dataset)
savefile = os.path.join(data_path, 'schmidtko_ocean_input.nc')

# the original data is provided for download here.
link = "http://www.geomar.de/fileadmin/personal/fb1/po/sschmidtko/Antarctic_shelf_data.txt"

# if data is not yet there, download
downloaded_file = os.path.join(data_path,"Antarctic_shelf_data.txt")
if not os.path.isfile(downloaded_file):
  print "Downloading Schmidtko data."
  os.system("mkdir " + data_path)
  os.system("wget -N " + link + " -P " + data_path)

# Read data from txt file:
lon = np.zeros(0)
lat = np.zeros(0)
depth = np.zeros(0)
c_t = np.zeros(0)
c_t_std = np.zeros(0)
abs_s = np.zeros(0)
abs_s_std = np.zeros(0)


f = open(os.path.join(data_path, 'Antarctic_shelf_data.txt'), 'r')

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

theta_ocean = np.zeros((1,len(lat_new), len(lon_new))) + fillvalue
salinity_ocean = np.zeros((1,len(lat_new), len(lon_new))) + fillvalue
height   = np.zeros((1,len(lat_new), len(lon_new))) + fillvalue

for i in range(len(c_t)):
    # go through all temp and sal vals and fill them into the right place
    ilat = np.in1d(lat_new.ravel(), lat[i]).reshape(lat_new.shape)
    comp_lon = lon[i] -360
    if comp_lon < -180:
        comp_lon = 360 + comp_lon
    ilon = np.in1d(lon_new.ravel(), comp_lon).reshape(lon_new.shape)
    theta_ocean[0,ilat,ilon] = c_t[i]
    salinity_ocean[0,ilat,ilon] = abs_s[i]
    height[0,ilat,ilon] = depth[i]

# Save data as NetCDF file

wrtfile = nc.Dataset(savefile, 'w', format='NETCDF4_CLASSIC')
wrtfile.createDimension('lon', size=len(lon_new))
wrtfile.createDimension('lat', size=len(lat_new))
wrtfile.createDimension('time', size=None)
wrtfile.createDimension('nv', size=2)

nct     = wrtfile.createVariable('time', 'float32', ('time',))
nctb    = wrtfile.createVariable('time_bnds', 'float32', ('time','nv'))
nclat   = wrtfile.createVariable('lat', 'f4', ('lat',))
nclon   = wrtfile.createVariable('lon', 'f4', ('lon',))
nctemp  = wrtfile.createVariable('theta_ocean', 'f4', ('time','lat', 'lon'))
ncsal   = wrtfile.createVariable('salinity_ocean', 'f4', ('time','lat', 'lon'))
nchgt   = wrtfile.createVariable('height', 'f4', ('time','lat', 'lon'))

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

ncsal[:] = salinity_ocean
ncsal.units = 'g/kg' # absolute salinity_ocean
nctemp[:] = theta_ocean
if(temp_in_kelvin):
    nctemp.units='Kelvin'
else:
    nctemp.units='degree_Celsius' # conservative temperature

nchgt[:] = height
nchgt.units = 'm'

wrtfile.close()


print 'Data successfully saved to', savefile
