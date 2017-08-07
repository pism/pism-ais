"""
matthias.mengel@pik

A template for time_averaging data.
Applied after remapping. May also be done before it.

"""

import os, sys
import subprocess
import netCDF4 as nc
## this hack is needed to import config.py from the project root
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if project_root not in sys.path: sys.path.append(project_root)
import config as cf; reload(cf)

dataset = "racmo_wessem"
grid = cf.grid_id
tap = [str(t) for t in cf.time_averaging_period]

data_path = os.path.join(cf.output_data_path, dataset)
inputfile = os.path.join(data_path, dataset+"_"+grid+".nc")
timemean_file = os.path.join(data_path, dataset+"_"+grid+"_mean"+tap[0]+"_"+tap[1]+".nc")

cmd = "cdo -O -f nc4c timmean -selyear,"+tap[0]+"/"+tap[1]+" "+inputfile+" "+timemean_file
print cmd
# if you see "cdo: error while loading shared libraries:",
# 'module load cdo'
subprocess.check_call(cmd,shell=True)

## Add the precipition from Racmo HadCM3 simulations over the ocean to the SMB field.
## This is a workaround to have some data over regions that are now ocean. We will update
## once we have the fields from Wessem et al. 2014.

precip_scenario = "c20" # A1B or c20
precip_dataset = "racmo_hadcm3_I2S"
data_path = os.path.join(cf.output_data_path, precip_dataset)

precip_timemean_file = os.path.join(data_path, precip_dataset+"_"+precip_scenario+
    "_"+cf.grid_id+"_mean"+tap[0]+"_"+tap[1]+".nc")

print precip_timemean_file
nchadcm3 = nc.Dataset(precip_timemean_file,"r")
ncw = nc.Dataset(timemean_file,"a")

precip_hadcm3 = nchadcm3.variables["precipitation"][:]
smb_wessem = ncw.variables["smb"][:]
smb_wessem.unshare_mask()
smb_wessem[smb_wessem.mask] = precip_hadcm3[smb_wessem.mask]

print ncw.variables["smb"][:].shape
ncw.variables["smb"][:] = smb_wessem

ncw.close()
nchadcm3.close()

print "added the precip field from"
print precip_timemean_file
print "to", timemean_file

# TODO: remove time dimension, time_bnds, assigned, dir, dtg
