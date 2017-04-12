
"""
matthias.mengel@pik, torsten.albrecht@pik
Regridding: bring your data to the grid we commonly use for PISM Antarctica
simulations. This is equivalent to the ALBMAP grid.
This step will take a while if high resolution data is processed.
Regrid Bedmap2 data to various grid resolution using cdo remapcony.

"""

import os, sys
import subprocess

## this hack is needed to import config.py from the project root
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if project_root not in sys.path: sys.path.append(project_root)
import config as cf; reload(cf)
import pism_input.pism_input as pi; reload(pi)

dataset="zwally_basins"
# resolution for the output file, set from config.
resolution = cf.resolution # in km
# conservative regridding for bedmap2 and albmap data. does
# not yet work for the other datasets.
use_conservative_regridding = False

data_path = os.path.join(cf.output_data_path, dataset)
inputfile = os.path.join(data_path, 'basins_zwally_5km_input.nc')
regridded_file = os.path.join(data_path, dataset+"_"+str(resolution)+"km.nc")

# the cdo target grids are independent of the specific input dataset.
# they are therefore created beforehand by grids/create_cdo_grid.py
cdo_targetgrid_file = os.path.join(cf.cdo_remapgridpath,'pism_'+str(int(resolution))+'km.nc')

# check if target grid is present.
if not os.path.isfile(cdo_targetgrid_file):
    print "cdo target grid file", cdo_targetgrid_file," does not exist."
    print "run grids/create_cdo_grid.py first."
    sys.exit(0)

## Integer regridding for basin values, can be run interactively.
subprocess.check_call("cdo remapnn,"+cdo_targetgrid_file+" "+inputfile+" "+
                      regridded_file, shell=True)
