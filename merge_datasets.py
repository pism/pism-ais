"""

Merge different datasets of the same resolution into one output file.
This output file should be ready for PISM to read.

This script normally does not nead user input, except for custom filenames.
All other user specific setting are in config.py

"""

import os
import subprocess
import config as cf; reload(cf)
import pism_input.pism_input as pi; reload(pi)

# will hold the merged files
data_path = os.path.join(cf.output_data_path, "merged")

# set your custom name here if standard naming is not used
# (for example if you applied a time mean in between)
custom_file_names = {"racmo_hadcm3_I2S":
                     "racmo_hadcm3_I2S_A1B_"+str(cf.resolution)+"km_timemean.nc",
                     "schmidtko":"schmidtko_"+str(cf.resolution)+"km_means.nc"}

merged_filename = ("_").join(cf.datasets_to_merge)+"_"+str(cf.resolution)+"km.nc"
merged_filename = os.path.join(data_path,merged_filename)
# created beforehand by grids/create_cdo_grid.py
cdo_targetgrid_file = os.path.join(cf.cdo_remapgridpath,'pism_'+str(cf.resolution)+'km.nc')


if not os.path.exists(data_path): os.makedirs(data_path)

preselected_datapaths = []

added_lat_lon_mapping = False

for ds in cf.datasets_to_merge:

    custom_file_name = custom_file_names[ds] if ds in custom_file_names else None

    input_datapath = pi.get_path_to_data(cf.output_data_path,ds,cf.resolution,
                                         custom_file_name=custom_file_name)

    # data with selected variables will be written to /merged directory
    preselected_datapath = pi.get_path_to_data(cf.output_data_path,ds,cf.resolution,
                                custom_data_path=data_path)
    selected_variables = ",".join(cf.variables[ds])

    # assume latitude and longitude are present in first dataset.
    # TODO: include mapping variable, need to ensure that it is present
    # in all input files first.
    if not added_lat_lon_mapping:
        selected_variables += ",lon,lat"
        added_lat_lon_mapping = True

    cmd = ("ncks -O -4 -v "+selected_variables+" "+
            input_datapath+" "+preselected_datapath)
    print cmd
    subprocess.check_call(cmd,shell=True)
    preselected_datapaths.append(preselected_datapath)

cmd = "cdo -O -f nc4c merge "+" ".join(preselected_datapaths)+" "+merged_filename
print cmd
# if you see "cdo: error while loading shared libraries:",
# 'module load cdo'
subprocess.check_call(cmd,shell=True)

# add the axes from the gridfile to the merged file
cmd = "ncks -A -v x,y "+cdo_targetgrid_file+" "+merged_filename
subprocess.check_call(cmd,shell=True)