"""
matthias.mengel@pik

A template for time_averaging data.
Applied after remapping. May also be done before it.

"""

import os, sys
import subprocess
## this hack is needed to import config.py from the project root
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if project_root not in sys.path: sys.path.append(project_root)
import config as cf; reload(cf)

dataset = "racmo_wessem"
grid = cf.grid_id

data_path = os.path.join(cf.output_data_path, dataset)
inputfile = os.path.join(data_path, dataset+"_"+grid+".nc")
timemean_file = inputfile.rstrip(".nc")+"_timemean.nc"


cmd = "cdo -O -f nc4c timmean "+inputfile+" "+timemean_file
print cmd
# if you see "cdo: error while loading shared libraries:",
# 'module load cdo'
subprocess.check_call(cmd,shell=True)

# TODO: remove time dimension, time_bnds, assigned, dir, dtg
