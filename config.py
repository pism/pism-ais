"""
This file contains user defined paths and settings.
Normally, this should not be committed unless your changes are major.
"""
import os
import pwd

authors="matthias.mengel@pik-potsdam.de and torsten.albrecht@pik-potsdam.de"

#torsten local and tumble
output_data_path = os.path.expanduser("/p/projects/tumble/pism_input/GitLab/")
#output_data_path = os.path.expanduser("/home/albrecht/Documents/pism/python/pism_input/")

#matthias
output_data_path = os.path.expanduser("~/data/20170316_PismInputData/")
output_data_path = "/p/projects/tumble/mengel/pismInputData/20170316_PismInputData"

# these datasets will be checked for the variables_to_write.
# a variable is taken from the first dataset in the list in which it is found,
# i.e. if datasets = ["bedmap2","albmap"], thk is taken from bedmap2 and air_temp from albmap.
# datasets should be named here same as the subfolder its preprocessing.
datasets = ["bedmap2","albmap"]

# will be written to output_data_path
output_file_name = "pism_antarctica.nc"
# these PISM variables will appear in the file created by create_pism_input.py
variables_to_write = ["thk","topg","precipitation","air_temp"]

# this resolution will be written by appear in the output file.
resolution = 5 # in km


#### No edits needed below that line. ####
output_file_name = os.path.join(output_data_path,output_file_name)
cdo_remapgridpath = os.path.join(output_data_path,"cdo_remapgrids")
project_root = os.path.dirname(os.path.abspath(__file__))
username = pwd.getpwuid(os.getuid()).pw_name