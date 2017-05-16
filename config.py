"""
This file contains user defined paths and settings.
Normally, this should not be committed unless your changes are major.
"""
import os
import pwd

authors="matthias.mengel@pik-potsdam.de and torsten.albrecht@pik-potsdam.de"

# the resolution of the final output files.
resolution = 15 # in km

#torsten local and tumble
#output_data_path = os.path.expanduser("/p/projects/tumble/pism_input/GitLab/")
#output_data_path = os.path.expanduser("/home/albrecht/Documents/pism/python/pism_input/")
output_data_path = os.path.expanduser("/p/projects/tumble/albrecht/pism_input/data/")


# matthias
#output_data_path = os.path.expanduser("~/data/20170316_PismInputData/")
#output_data_path = "/p/projects/tumble/mengel/pismInputData/20170316_PismInputData"

# RACMO data is not freely available and cannot be downloaded,
# so we have to provide an explicit path here
# if racmo data is intended to be used in publications,
# get in contact the racmo authors before.
racmo_data_path = "/p/projects/tumble/mengel/pismSourceData/20170328_RacmoHadCM3_Ice2Sea"

# merge the follwing dataset into one PISM-ready file.
# datasets should be named here as the subfolder of its preprocessing.
# TODO: merging schmidtko data does not work yet due to subtle grid differences.
#       it can be supplied through -ocean cavity -ocean_cavity_file $file directly.
datasets_to_merge = ["bedmap2","albmap","racmo_hadcm3_I2S","schmidtko"]

# choose here which variables should be taken from which dataset
# The basins variable for the PICO model can be passed with the Schmidtko dataset.
variables = {"bedmap2":["thk","topg","usurf"],
             "albmap":["bheatflx"],
             "racmo_hadcm3_I2S":["precipitation","air_temp"],
             "schmidtko":["theta_ocean","salinity_ocean","basins"]}

#### No edits needed below that line. ####
cdo_remapgridpath = os.path.join(output_data_path,"cdo_remapgrids")
project_root = os.path.dirname(os.path.abspath(__file__))
username = pwd.getpwuid(os.getuid()).pw_name
