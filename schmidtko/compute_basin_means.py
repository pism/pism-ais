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
import subprocess
import datetime

## this hack is needed to import config.py from the project root
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if project_root not in sys.path: sys.path.append(project_root)
import config as cf; reload(cf)
import pism_input.pism_input as pi; reload(pi)

# save data path
dataset="schmidtko"

data_path = os.path.join(cf.output_data_path, dataset)
infile = os.path.join(data_path, 'schmidtko_'+cf.grid_id+'.nc')
outfile = os.path.join(data_path, 'schmidtko_'+cf.grid_id+'_means.nc')
print infile

zwally_data_path = os.path.join(cf.output_data_path,"basins_icesat_zwally")
basinfile = os.path.join(zwally_data_path, 'basins_icesat_zwally'+cf.grid_id+'.nc')


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


now = datetime.datetime.now().strftime("%B %d, %Y")
wrtfile.Descricption = "Antarctic drainage basins mapped by NASA and modified. Temperature (converted to potential) and salinity averaged over the basins at the depth of the continential shelf." ;
wrtfile.Reference = "Basins from Zwally, H. Jay, Mario B. Giovinetto, Matthew A. Beckley, and Jack L. Saba, 2012, Antarctic and Greenland Drainage  Systems, GSFC Cryospheric Sciences Laboratory, at http://icesat4.gsfc.nasa.gov/cryo_data/ant_grn_drainage_systems.php. Temperature, Salinity from Schmidtko, S., Heywood, K. J., Thompson, A. F., & Aoki, S. (2014). Multidecadal warming of Antarctic waters. Science, 346(6214), 1227-1231. ,at http://www.geomar.de/fileadmin/personal/fb1/po/sschmidtko/Antarctic_shelf_data.txt" ;
wrtfile.proj4 = "+lon_0=0.0 +ellps=WGS84 +datum=WGS84 +lat_ts=-71.0 +proj=stere +x_0=0.0 +units=m +y_0=0.0 +lat_0=-90.0"
wrtfile.comment  = cf.authors+" created netcdf file at " + now


wrtfile.close()

print "Basin mean values successfully written to"
print outfile


## TODO: the following does not fully solve the problem for cdo merge this file with the
##       others. Should we just not put time and depth dimension in the creation of
##       the Schmidtko data at first place?
## Remove the time dimesion
subprocess.check_call("ncwa -O -a time "+outfile+" "+outfile, shell=True)
## delete the time and the height variable
subprocess.check_call("ncks -O -C -x -v time "+outfile+" "+outfile, shell=True)
subprocess.check_call("ncks -O -C -x -v height "+outfile+" "+outfile, shell=True)
