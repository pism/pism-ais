
"""
matthias.mengel@pik, torsten.albrecht@pik
Regridding: bring your data to the grid we commonly use for PISM Antarctica
simulations.
This step will take a while if high resolution data is processed.
Regrid Bedmap2 data to various grid resolution using cdo remapcony.

"""

import os, sys
import subprocess

## this hack is needed to import config.py from the project root
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if project_root not in sys.path: sys.path.append(project_root)
from importlib import reload
import config as cf; reload(cf)
import pism_input.pism_input as pi; reload(pi)

dataset="icesat4_basins"
# resolution for the output file, set from config.
# resolution = cf.resolution # in km
# conservative regridding for bedmap2 and albmap data. does
# not yet work for the other datasets.
use_conservative_regridding = False

data_path = os.path.join(cf.output_data_path, dataset)
inputfile = os.path.join(data_path, 'icesat4_basins_1km_input.nc')
# pi.prepare_ncfile_for_cdo(inputfile)

cdo_targetgrid_file, regridded_file = pi.get_filenames_for_cdo(
    cf.cdo_remapgridpath, data_path, dataset, cf.grid_id)

# Create a bash script that handles the regridding.
# Regridding can be a CPU-heavy task. Choose cluster_regridding=True in config.py
# If you want to submit to the cluster using SLURM.
# use 'sbatch cdo_remap.sh' to submit your job.
# Conservative regridding does not work for all datasets yet, use it for bedmap2 or albmap.
# We use cdo, see https://code.zmaw.de/projects/cdo/embedded/index.html

# FIXME: use remapnn instead of remapbil in script.
pi.write_regrid_command_file(cf, data_path, dataset, inputfile, cf.grid_id,
                     cdo_targetgrid_file, regridded_file, cf.regridding_method)


# ## Integer regridding for basin values, can be run interactively.
# subprocess.check_call("cdo remapnn,"+cdo_targetgrid_file+" "+inputfile+" "+
#                       regridded_file, shell=True)
