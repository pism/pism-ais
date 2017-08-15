"""
    matthias.mengel@pik, torsten.albrecht@pik, ronja.reese@pik

    Merge the two initmip-provided anomaly forcing files of surface mass balance and
    basal melt rate into one file. Initmip provides these on 1km and 16km resolution.
    All other resolutions should to be conservatively regridded from the 1km file.

    # see also for Greenland: https://github.com/pism/pism-gris/blob/master/initMIP/prepare_anomalies.sh
    # by Andy Aschwanden (UAF)

    #Prepare data for initMIP Antarctica
    #Wiki: http://www.climate-cryosphere.org/wiki/index.php?title=InitMIP-Antarctica

    #Downloaded dBasalMelt and dSMB anomaly fields from
    #ftp searise@cryoftp1.gsfc.nasa.gov initMIP directory /ISMIP6/initMIP/AIS

    #Password personal communication with Sophie Nowicki <sophie.nowicki@nasa.gov>

    #official email: ismip6 <ismip6@gmail.com>
    #cc: helene.seroussi@jpl.nasa.gov
"""

import os, sys
import subprocess
import netCDF4 as nc
import numpy as np
import shutil

## this hack is needed to import config.py from the project root
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if project_root not in sys.path: sys.path.append(project_root)
import config as cf; reload(cf)
import pism_input.pism_input as pi; reload(pi)

dataset = "initmip"
data_path = os.path.join(cf.output_data_path, dataset)

if not os.path.exists(data_path): os.makedirs(data_path)

for grid_id in ["1km","16km"]:

    output_file = os.path.join(data_path, dataset+"_"+grid_id+"_input.nc")

    source_file = {"smb": os.path.join(cf.initmip_data_path,'dSMB/smb_anomaly_'+grid_id+'.nc'),
                   "bmr": os.path.join(cf.initmip_data_path,'dBasalMelt/basal_melt_anomaly_'+grid_id+'.nc')
                   }

    shutil.copyfile(source_file["smb"], output_file)
    # merge_these_files = " ".join([source_file[var] for var in source_file.keys()])

    subprocess.check_call('ncks -A -v abmb '+source_file["bmr"]+" "+output_file, shell=True)

    print output_file, "created."

# for one kilometer file, correct the value at the sout pole to lat=-90
grid_id = "1km"
output_file = os.path.join(data_path, dataset+"_"+grid_id+"_input.nc")
ncf = nc.Dataset(output_file,"a")
# south pole point
ncf.variables["lat"][3040,3040] = -90.
ncf.close()

# also prepare for the 1km file remapping to 4km, 8km etc.
pi.prepare_ncfile_for_cdo(output_file)