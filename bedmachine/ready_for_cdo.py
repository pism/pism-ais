import os, sys
import numpy as np
import sys
import netCDF4 as nc
import datetime
import subprocess


## this hack is needed to import config.py from the project root
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if project_root not in sys.path: sys.path.append(project_root)
import config as cf; reload(cf)
import pism_input.pism_input as pi; reload(pi)

bedmachine_data_path = os.path.join(cf.output_data_path, "bedmachine")
ncout_name = os.path.join(bedmachine_data_path, filename+"_ReadyForRemapping.nc")

print "Prepare for cdo "+ncout_name
pi.prepare_ncfile_for_cdo(ncout_name)


