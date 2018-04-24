"""
Download Southern Ocean data from Schmidtko et al, Science 2015
http://science.sciencemag.org/content/346/6214/1227
and create NetCDF file from the txt data.
ronja.reese@pik-potsdam.de, matthias.mengel@pik

"""

import os, sys
import numpy as np
import numpy.ma as ma
import netCDF4 as nc
import scipy.io
## this hack is needed to import config.py from the project root
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if project_root not in sys.path: sys.path.append(project_root)
import config as cf; reload(cf)

# save data path
dataset="shimada_socean"

temp_in_kelvin = False

data_path = os.path.join(cf.output_data_path, dataset)
savefile = os.path.join(data_path, 'shimada_ocean_input.nc')
if not os.path.exists(data_path): os.makedirs(data_path)

input_file = os.path.join(cf.shimada_socean_data_path,"climatology.mat")
data = scipy.io.loadmat(input_file)


wrtfile = nc.Dataset(savefile, 'w', format='NETCDF4_CLASSIC')
wrtfile.createDimension('time', size=None)
wrtfile.createDimension('lon', size=data["lat"].shape[1])
wrtfile.createDimension('lat', size=data["lat"].shape[0])
wrtfile.createDimension('depth', size=data["depth"].shape[1])
wrtfile.createDimension('nv', size=2)

nct     = wrtfile.createVariable('time', 'float32', ('time',))
nctb    = wrtfile.createVariable('time_bnds', 'float32', ('time','nv'))
ncdep   = wrtfile.createVariable('depth', 'f4', ('depth',))
nclat   = wrtfile.createVariable('lat', 'f4', ('lat',))
nclon   = wrtfile.createVariable('lon', 'f4', ('lon',))
nctemp  = wrtfile.createVariable('theta_ocean', 'f4', ('time','depth','lat','lon'))
ncsal   = wrtfile.createVariable('salinity_ocean', 'f4', ('time','depth','lat', 'lon'))
# nchgt   = wrtfile.createVariable('height', 'f4', ('time','lat', 'lon'))

tm = np.arange(1,1.5,1)
nct[:] = tm
nct.units = 'years since 01-01-01 00:00:00'
nct.bounds = 'time_bnds'

t_bds = np.zeros([len(tm),2])
t_bds[:,0] = tm - 1
t_bds[0,:] = tm + 1
nctb[:] = t_bds
nctb.units = nct.units

nclon[:] = data["lon"][0,:]
nclon.units = 'degrees_east'
nclat[:] = data["lat"][:,0]
nclat.units = 'degrees_north'
ncdep[:] = data["depth"][0,:]
ncdep.units = 'm'

ncsal[:] = np.moveaxis(data["sal"],2,0)[np.newaxis,:]
ncsal.units = 'g/kg' # absolute salinity_ocean
nctemp[:] = np.moveaxis(data["ptm"],2,0)[np.newaxis,:]
if(temp_in_kelvin):
    nctemp.units='Kelvin'
else:
    nctemp.units='degree_Celsius' # conservative temperature

# nchgt[:] = data["w_depth"][np.newaxis,:]
# nchgt.units = 'm'

wrtfile.close()


print 'Data successfully saved to', savefile
