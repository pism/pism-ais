
# We here provide examples on how to regrid and merge input data for PISM.
# This is divided into steps. Comment out what you do not need.

import os
import jinja2
import config as cf; reload(cf)
import pism_input.pism_input as pi; reload(pi)


resolution = 5 # in km

# First step: download and prepare your input data with the specific scripts for
# for each datasets. These scripts are located in the data subfolders, i.e. bedmap2/
# albmap/, racmo/ etc.

# Second step: regridding: bring your data to the grid we commonly use for PISM Antarctica
# simulations. This is equivalent to the ALBMAP grid.
# This step will take a while if high resolution data is processed.
# We here provide the example for bedmap2, adapt the paths to your dataset.

# prepare the input file for cdo remapping
# this step takes a while for high resolution data (i.e. 1km)
data_path = os.path.join(cf.output_data_path, "bedmap2")
inputfile = os.path.join(data_path, 'bedmap2_1km_input.nc')
pi.prepare_ncfile_for_cdo(inputfile)

# Create the target grid on which we want to regrid on with cdo. The cdo_targetgrid_file
# is resolution-specific. The create_grid_for_cdo_remap is specific for albmap-like
# grids used in PISM.
cdo_targetgrid_file = pi.create_grid_for_cdo_remap(
    cf.cdo_remapgridpath,resolution=resolution)
pi.prepare_ncfile_for_cdo(cdo_targetgrid_file)

# Regridding is generally a CPU-heavy task. We therefore do not regrid interactively,
# but prepare a script that can be submitted to compute-clusters. The example is specific
# for PIK's cluster using SLURM.
# use 'sbatch cdo_remap.sh' to submit your job.
# Conservative regridding does not work for all datasets yet, use it for bedmap2 or albmap.
# We use cdo, see https://code.zmaw.de/projects/cdo/embedded/index.html

# make jinja aware of templates in the pism_input/tools folder
jinja_env = jinja2.Environment(loader=jinja2.FileSystemLoader(
            searchpath=os.path.join(cf.project_root,"tools")))

scen_template_file = "GENERATED_SCENARIO.SCEN.template"
scen_template = jinja_env.get_template("cdo_remap.sh.template")

regridded_file = os.path.join(data_path, "regridded_"+str(resolution)+"_km.nc")
mapweights = os.path.join(data_path, "mapweights.nc")
use_conservative_regridding = True

out = scen_template.render(user="mengel",
                           use_conservative_regridding = use_conservative_regridding,
                           targetgrid = cdo_targetgrid_file,
                           inputfile = inputfile,
                           mapweights = mapweights,
                           regridded_file = regridded_file,
                          )

with open("tools/cdo_remap.sh", 'w') as f:
    f.write(out)