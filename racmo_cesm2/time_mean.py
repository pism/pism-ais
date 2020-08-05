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

dataset = "racmo_cesm2"
grid = cf.grid_id
tap = [str(t) for t in cf.time_averaging_period]

data_path = os.path.join(cf.output_data_path, dataset)
inputfile = os.path.join(data_path, dataset+"_"+grid+".nc")
timemean_file = os.path.join(data_path, dataset+"_"+grid+"_mean"+tap[0]+"_"+tap[1]+".nc")

cmd = "cdo -O -f nc4c timmean -selyear,"+tap[0]+"/"+tap[1]+" "+inputfile+" "+timemean_file

# if you see "cdo: error while loading shared libraries:",
# 'module load cdo'
subprocess.check_call(cmd,shell=True)

