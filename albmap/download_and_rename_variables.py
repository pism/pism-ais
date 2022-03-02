"""
matthias.mengel@pik, torsten.albrecht@pik
Download Albmap data and save rename variables to make the PISM compatible.
ALBMAP is documented here:
http://www.earth-syst-sci-data.net/2/247/2010/
"""

import os, sys
import subprocess
## this hack is needed to import config.py from the project root
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if project_root not in sys.path: sys.path.append(project_root)
import config as cf; reload(cf)
import pism_input.pism_input as pi; reload(pi)

dataset = "albmap"
data_path = os.path.join(cf.output_data_path, dataset)
final_filename = os.path.join(data_path,"albmap_5km_input.nc")

# Documentation of the data: https://www.bas.ac.uk/project/bedmap-2/
link = "http://websrv.cs.umt.edu/isis/images/4/4d/Antarctica_5km_dev1.0.nc"
## such file definitions should go to config.py, so that other functions can access them.

# if data is not yet there, download
downloaded_file = os.path.join(data_path,"Antarctica_5km_dev1.0.nc")
if not os.path.isfile(downloaded_file):
  print("Downloading albmap data.")
  os.system("mkdir " + data_path)
  os.system("wget -N " + link + " -P " + data_path)


subprocess.check_call("cp "+downloaded_file+" "+final_filename,shell=True)
# following use NCO (http://nco.sourceforge.net/)
# rename dimensions
subprocess.check_call("ncrename -O -v x1,x -v y1,y -d x1,x -d y1,y "+final_filename,shell=True)
subprocess.check_call("ncrename -O -v time,t -d time,t "+final_filename,shell=True)
# fix polar stereographic parameter
subprocess.check_call("ncatted -O -a standard_parallel,mapping,m,d,-71.0 "+final_filename,shell=True)
# rename usurf for convenience
subprocess.check_call("ncrename -O -v usrf,usurf "+final_filename,shell=True)
# fix surface temperature name and make K
subprocess.check_call('ncap2 -O -s "air_temp=temp+273.15" '+final_filename+' '+final_filename,shell=True)
# choose Van de Berg et al version of accumulation; will treat as
# ice-equivalent snow rate and convert from an ice-equivalent
# thickness rate ("m year-1") to "kg m-2 year-1" by assuming ice
# density of 910 kg m-3
subprocess.check_call('ncap2 -O -s "precipitation=accr*910.0" '+final_filename+" "+
                      final_filename,shell=True)
subprocess.check_call('ncatted -O -a units,precipitation,m,c,"kg m-2 year-1" '+final_filename,shell=True)

# use bheatflx_shapiro as the default bheatflx data
subprocess.check_call('ncrename -O -v bheatflx_shapiro,bheatflx '+final_filename,shell=True)

subprocess.check_call('ncatted -O -a units,bheatflx,m,c,"W m-2" '+final_filename,shell=True)

# delete incorrect standard_name attribute from bheatflx; there is no known standard_name
subprocess.check_call('ncatted -O -a standard_name,bheatflx,d,, '+final_filename,shell=True)

# keep only the fields we actually use at bootstrapping
subprocess.check_call('ncks -O -v x,y,lat,lon,bheatflx,topg,thk,precipitation,air_temp,mapping '+
                      final_filename+' '+final_filename,shell=True)

# remove the time dimension
subprocess.check_call('ncwa -O -a t '+ final_filename+' '+final_filename,shell=True)

# remove the time (t) coordinate variable
subprocess.check_call('ncks -O -x -v t '+ final_filename+' '+final_filename,shell=True)
# create netcdf4 format
subprocess.check_call('ncks -O -4 '+ final_filename+' '+final_filename,shell=True)

pi.prepare_ncfile_for_cdo(final_filename)


# projection metadata
#subprocess.check_call('ncatted -O -a proj4,global,a,c,"+lon_0=0.0 +ellps=WGS84 +datum=WGS84 +lat_ts=-71.0 +proj=stere +x_0=0.0 +units=m +y_0=0.0 +lat_0=-90.0"'+ final_filename,shell=True)


print("  PISM-readable file",final_filename,"successfully created.")
