
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

dataset="zwally_basins"

data_path = os.path.join(cf.output_data_path, dataset)
inputfile = os.path.join(data_path, 'basins_zwally_5km_input.nc')
pi.prepare_ncfile_for_cdo(inputfile)

cdo_targetgrid_file, regridded_file = pi.get_filenames_for_cdo(
    cf.cdo_remapgridpath, data_path, dataset, cf.grid_id)

## Integer regridding for basin values, can be run interactively.
subprocess.check_call("cdo remapnn,"+cdo_targetgrid_file+" "+inputfile+" "+
                      regridded_file, shell=True)
