import os
import subprocess
import config as cf; reload(cf)
import pism_input.pism_input as pi; reload(pi)

# will hold the merged files
data_path = os.path.join(cf.output_data_path, "merged")

# set your custom name here if standard naming is not used
custom_file_names = {"racmo_hadcm3_I2S":"racmo_hadcm3_c20_timemean_input.nc"}

merged_filename = ("_").join(cf.datasets_to_merge)+"_"+str(cf.resolution)+"km.nc"
merged_filename = os.path.join(data_path,merged_filename)

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

cmd = "cdo -O merge "+" ".join(preselected_datapaths)+" "+merged_filename
print cmd
# if you see "cdo: error while loading shared libraries:",
# 'module load cdo'
subprocess.check_call(cmd,shell=True)
