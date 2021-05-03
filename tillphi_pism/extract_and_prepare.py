"""
matthias.mengel@pik, torsten.albrecht@pik
Extract the till friction angle (tillphi) from a PISM output file
and prepare for remapping with cdo.
"""

# import os, glob
# import numpy as np
# import numpy.ma as ma
# import sys, csv, datetime
# import netCDF4 as nc
# import datetime, math
import os, sys
import subprocess
import netCDF4 as nc

## this hack is needed to import config.py from the project root
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if project_root not in sys.path: sys.path.append(project_root)
import config as cf; reload(cf)
import pism_input.pism_input as pi; reload(pi)

dataset = "tillphi_pism"

data_path = os.path.join(cf.output_data_path, dataset)
if not os.path.exists(data_path): os.makedirs(data_path)

savefile = os.path.join(data_path, dataset)+".nc"

subprocess.check_call("ncks -O -4 -v tillphi "+cf.tillphi_data_path+" "+
                      os.path.join(data_path,dataset+"temp.nc"),shell=True)
# ncks -O -4 -v tillphi cf.tillphi_data_path tillphi_15km.nc
subprocess.check_call('ncap2 -O -s "lon_bnds=double(lon_bnds);lat_bnds=double(lat_bnds)" '+
                      os.path.join(data_path,dataset+"temp.nc")+" "+
                      savefile,shell=True)

subprocess.check_call('ncatted -O -a reference,tillphi,c,c,"optimized till friction angle as described in Albrecht et al., 2020a, https://doi.org/10.5194/tc-14-599-2020" '+ 
                        savefile, shell=True)


pi.prepare_ncfile_for_cdo(os.path.join(data_path,dataset+".nc"))
