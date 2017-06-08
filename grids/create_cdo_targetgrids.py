
"""
Create grids based on the ALBMAP dimensions, that are compatible with PISM.
That means, if they are fed into PISM, PISM will not need to regrid internally.
The here created grids are used as input for the cdo remapping, which is done
for each dataset separately. See for example bedmap2/remap.py

"""

import os, sys
import numpy as np
from pyproj import Proj
import netCDF4

## this hack is needed to import config.py from the project root
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if project_root not in sys.path: sys.path.append(project_root)
import config as cf; reload(cf)
import pism_input.pism_input as pi; reload(pi)

# if true, PISM will not need to regrid
use_PISM_grid=True


if __name__ == "__main__":

    for resolution in [1,2,5,10,12,15,20,30,50]:

        # cdo_targetgrid_file = os.path.join(cf.cdo_remapgridpath,
        #                         'pism_'+str(int(resolution))+'km.nc')
        gridfile = pi.create_grid_for_cdo_remap(cf.cdo_remapgridpath,
                        use_PISM_grid=True, resolution=resolution)
        pi.prepare_ncfile_for_cdo(gridfile)
