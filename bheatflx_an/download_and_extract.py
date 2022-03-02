"""
matthias.mengel@pik, torsten.albrecht@pik
Download bheatflux data netcdf file.
"""

import os, sys
import numpy as np
import sys
import netCDF4
import datetime

## this hack is needed to import config.py from the project root
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if project_root not in sys.path: sys.path.append(project_root)
from importlib import reload
import config as cf; reload(cf)
import pism_input.pism_input as pi; reload(pi)

### Bheatflux from An at al. 2015  ##########################################################
# Documentation of the data: http://www.seismolab.org/model/antarctica/lithosphere/#an1-hf
link="http://www.seismolab.org/model/antarctica/lithosphere/AN1-HF.tar.gz"
data_path = os.path.join(cf.output_data_path, "bheatflx_an")
ncout_name = os.path.join(data_path, 'bheatflux_an_input.nc')

# if data is not yet downloaded
if not os.path.exists(ncout_name):
  print("Downloading data on lon-lat grid.")
  os.system("mkdir -p " + data_path)
  os.system("wget " + link + " -P " + data_path)
  print("data_path:" + data_path)
  print("ncout_name:" + ncout_name)
  os.system("cd " + data_path) 
  os.system("tar -zxvf " + os.path.join(data_path, link.split('/')[-1]))
  os.system("mv AN1-HF.grd " + ncout_name)


now = datetime.datetime.now().strftime("%B %d, %Y")

ncout = netCDF4.Dataset(ncout_name,"a",format='NETCDF4_CLASSIC')
#ncout.proj4 = cf.proj4str
ncout.comment  = cf.authors+" downloaded netcdf file at " + now +" from "+link
ncout.close()

# prepare the input file for cdo remapping
# this step takes a while for high resolution data (i.e. 1km)
#pi.prepare_ncfile_for_cdo(ncout_name)

