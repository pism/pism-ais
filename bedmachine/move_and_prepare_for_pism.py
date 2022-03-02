import os, sys
import numpy as np
import sys
import netCDF4 as nc
import datetime
import subprocess


## this hack is needed to import config.py from the project root
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if project_root not in sys.path: sys.path.append(project_root)
from importlib import reload
import config as cf; reload(cf)
import pism_input.pism_input as pi; reload(pi)

filename = "BedMachineAntarctica_2020-07-15_v02"

bedmachine_data_path = os.path.join(cf.output_data_path, "bedmachine")
ncout_name = os.path.join(bedmachine_data_path, filename+"_ReadyForRemapping.nc")

# move file to output path
if not os.path.exists(os.path.join(bedmachine_data_path)):
  print("Moving bedmachine data.")
  os.system("mkdir " + bedmachine_data_path)
subprocess.call("mv %s %s" % (filename+'.nc', bedmachine_data_path), shell=True)
subprocess.call("mv %s %s" % (filename+'.nc.xml', bedmachine_data_path), shell=True)

# copy file 
subprocess.check_call("cp %s %s " % (os.path.join(bedmachine_data_path,filename+'.nc'),ncout_name) ,shell=True)

# variables do not need to be renamed as PISM recongnizes standard names

# raise Lake Vostok
print("Adjusting bed to remove lake Vostok.")

ncf = nc.Dataset(ncout_name, 'r+')
topg_var = ncf.variables["bed"]
mask = ncf.variables["mask"][:]
surface = ncf.variables["surface"][:]
thickness = ncf.variables["thickness"][:]

topg = topg_var[:]
LAKE_VOSTOK=4
print("Before ", np.mean(topg[mask==LAKE_VOSTOK]))
topg[mask==LAKE_VOSTOK] = surface[mask==LAKE_VOSTOK] - thickness[mask==LAKE_VOSTOK] 
print "After ", np.mean(topg[mask==LAKE_VOSTOK])
topg_var[:] = topg

# add projection string
ncf.proj4 = cf.proj4str

ncf.close()


# keep only the fields we actually use at bootstrapping
subprocess.check_call('ncks -O -v x,y,bed,surface,mask,thickness '+
                      os.path.join(bedmachine_data_path,ncout_name)+' '+os.path.join(bedmachine_data_path,ncout_name),shell=True)



