
"""
matthias.mengel@pik, torsten.albrecht@pik
Regridding: bring your data to the grid we commonly use for PISM Antarctica
simulations. This is equivalent to the ALBMAP grid.
This step will take a while if high resolution data is processed.
Regrid Bedmap2 data to various grid resolution using cdo remapcony.

"""

import os, sys
import jinja2

## this hack is needed to import config.py from the project root
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if project_root not in sys.path: sys.path.append(project_root)
import config as cf; reload(cf)
import pism_input.pism_input as pi; reload(pi)

scenario = "c20" # A1B or c20
dataset = "racmo_hadcm3_I2S"

# regridding_method: choose "bilinear", "integer" or "conservative"
regridding_method = "bilinear"

data_path = os.path.join(cf.output_data_path, dataset)
inputfile = os.path.join(data_path, 'racmo_hadcm3_'+scenario+'_input.nc')

cdo_targetgrid_file, regridded_file = pi.get_filenames_for_cdo(
    cf.cdo_remapgridpath, data_path, dataset, cf.grid_id)

regridded_file = os.path.join(data_path, dataset+"_"+scenario+"_"+cf.grid_id+".nc")

# Create a bash script that handles the regridding.
# Regridding can be a CPU-heavy task. Choose cluster_regridding=True in config.py
# If you want to submit to the cluster using SLURM.
# use 'sbatch cdo_remap.sh' to submit your job.
# Conservative regridding does not work for all datasets yet, use it for bedmap2 or albmap.
# We use cdo, see https://code.zmaw.de/projects/cdo/embedded/index.html
pi.write_regrid_command_file(cf, data_path, dataset, inputfile, cf.grid_id,
                     cdo_targetgrid_file, regridded_file, regridding_method)