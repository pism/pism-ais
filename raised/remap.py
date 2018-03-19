
"""
matthias.mengel@pik-potsdam.de, torsten.albrecht@pik-potsdam.de
Regridding: bring your data to the grid we commonly use for PISM Antarctica
simulations.
Regrid RAISED data to various grid resolution using cdo remapnn.

"""

import os, sys

## this hack is needed to import config.py from the project root
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if project_root not in sys.path: sys.path.append(project_root)
import config as cf; reload(cf)
import pism_input.pism_input as pi; reload(pi)


dataset="raised"

data_path = os.path.join(cf.output_data_path, dataset)

cdo_targetgrid_file, regridded_file = pi.get_filenames_for_cdo(
          cf.cdo_remapgridpath, data_path, dataset, cf.grid_id)

os.system("export REMAP_EXTRAPOLATE=off")

for scenario in ["a","b"]:
  for snaps in [20,15,10,5]:

    infile = os.path.join(data_path,"raised_"+str(snaps)+"ka_scenario"+scenario+"_1kmmask.nc")
    outfile = regridded_file.replace(dataset+"_",dataset+"_"+str(snaps)+"ka_scenario"+scenario+"_")

    pi.prepare_ncfile_for_cdo(infile)
    cmd = "cdo -P 4 remapnn,"+cdo_targetgrid_file+" "+infile+" "+outfile
    print cmd
    os.system(cmd)


#export REMAP_EXTRAPOLATE=off

#for scenario in a b ;
#do
#  for snaps in 20 15 10 5 ;
#  do
#    echo $scenario
#    python nc2cdo.py data/raised_${snaps}ka_scenario${scenario}_1kmmask.nc
#    cdo -P 4 remapnn,${resfile} data/raised_${snaps}ka_scenario${scenario}_1kmmask.nc raised_${res}/raised_${snaps}_${scenario}_${res}.nc
#  done
#done
