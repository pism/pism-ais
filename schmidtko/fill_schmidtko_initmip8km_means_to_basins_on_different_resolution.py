"""
ronja.reese@pik
Regridding: bring your data to the grid we commonly use for PISM Antarctica
simulations.
This routine uses means_per_basin for 8km initmip and fills them to basins on any 
other grid. Start with low resolution regridding and use this file to move to high
resolution

"""

import os, sys
import jinja2
import shutil
import numpy as np
import netCDF4 as nc
import datetime

## this hack is needed to import config.py from the project root
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if project_root not in sys.path: sys.path.append(project_root)
import config as cf; reload(cf)
import pism_input.pism_input as pi; reload(pi)

dataset="schmidtko"

data_path = os.path.join(cf.output_data_path, dataset)
inputfile = os.path.join(data_path, 'schmidtko_initmip8km_means.nc')

targetgrid_file = os.path.join(cf.cdo_remapgridpath, cf.grid_id+'.nc')
regridded_file = os.path.join(data_path, dataset+"_"+cf.grid_id+"_means.nc")
basins_file = os.path.join(cf.output_data_path,'zwally_basins',"zwally_basins"+"_"+cf.grid_id+".nc") 
basins_input_file = os.path.join(cf.output_data_path,'zwally_basins',"zwally_basins"+"_initmip8km.nc")      


# copy the targetgrid file to the output file
shutil.copyfile(basins_file,regridded_file)

# fill the basins with the corresponding value of schmidtko temps and salinities
# get mean salinity / temperature values

ncf = nc.Dataset(inputfile)
thetao8 = ncf.variables["theta_ocean"][:]
salinityo8 = ncf.variables["salinity_ocean"][:]
ncf.close()

ncf = nc.Dataset(basins_input_file)
basins8 = ncf.variables["basins"][:]
ncf.close()

nbasins = int(np.max(basins8))
temps = np.zeros(nbasins+1)
sals = np.zeros(nbasins+1)

for basin in range(1,nbasins+1):
	#print basin
	temps[basin] = np.mean(thetao8[basins8==basin])
	sals[basin] = np.mean(salinityo8[basins8==basin])
	
print temps
print sals                            
         
# distribute averaged values over basins of new grid
ncf = nc.Dataset(basins_file)
basins = ncf.variables["basins"][:]
ncf.close()

thetao = np.copy(basins)*0.0
salinity = np.copy(basins)*0.0

for basin in range(1,nbasins+1):
	thetao[basins==basin] = temps[basin]
	salinity[basins==basin] = sals[basin]	


# open new file, add variable
ncout = nc.Dataset(regridded_file,'a')
ncv = ncout.createVariable( varname="theta_ocean",datatype='float32',dimensions=('y','x') )
ncv[:] = thetao
ncv = ncout.createVariable( varname="salinity_ocean",datatype='float32',dimensions=('y','x') )
ncv[:] = salinity
ncout.variables['theta_ocean'].units = "Celsius"
ncout.variables['salinity_ocean'].units = "g/kg"

basinsvar = ncout.variables['basins']
basinsvar.long_name = "drainage basins"
basinsvar.standard_name = "drainage_basins"


now = datetime.datetime.now().strftime("%B %d, %Y")
ncout.Descricption = "Antarctic drainage basins mapped by NASA and modified. Temperature (converted to potential) and salinity averaged over the basins at the depth of the continential shelf." ;
ncout.Reference = "Basins from Zwally, H. Jay, Mario B. Giovinetto, Matthew A. Beckley, and Jack L. Saba, 2012, Antarctic and Greenland Drainage  Systems, GSFC Cryospheric Sciences Laboratory, at http://icesat4.gsfc.nasa.gov/cryo_data/ant_grn_drainage_systems.php. Temperature, Salinity from Schmidtko, S., Heywood, K. J., Thompson, A. F., & Aoki, S. (2014). Multidecadal warming of Antarctic waters. Science, 346(6214), 1227-1231. ,at http://www.geomar.de/fileadmin/personal/fb1/po/sschmidtko/Antarctic_shelf_data.txt" ;
ncout.proj4 = "+lon_0=0.0 +ellps=WGS84 +datum=WGS84 +lat_ts=-71.0 +proj=stere +x_0=0.0 +units=m +y_0=0.0 +lat_0=-90.0"
ncout.comment  = cf.authors+" created netcdf file at " + now

ncout.close()

## TODO: the following does not fully solve the problem for cdo merge this file with the
##       others. Should we just not put time and depth dimension in the creation of
##       the Schmidtko data at first place?
## Remove the time dimesion
#subprocess.check_call("ncwa -O -a time "+outfile+" "+outfile, shell=True)
## delete the time and the height variable
#subprocess.check_call("ncks -O -C -x -v time "+outfile+" "+outfile, shell=True)
#subprocess.check_call("ncks -O -C -x -v height "+outfile+" "+outfile, shell=True)

pi.prepare_ncfile_for_cdo(regridded_file)


                                                                                                                      
