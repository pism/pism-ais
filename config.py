"""
This file contains user defined paths and settings.
Normally, this should not be committed unless your changes are major.
"""
import os
import pwd

authors="matthias.mengel@pik-potsdam.de and torsten.albrecht@pik-potsdam.de"

# the resolution of the final output files.
resolution = 16 # in km

grid_id = "initmip8km"

# grids, as inferred from PISM output
grids={
    # grids, as inferred from PISM output
    "pism50km":[120,120,-2777500,-2777500,3172500,3172500],
    "pism30km":[200,200,-2787500,-2787500,3182500,3182500],
    "pism20km":[300,300,-2792500,-2792500,3187500,3187500],
    "pism15km":[400,400,-2795000,-2795000,3190000,3190000],
    "pism12km":[500,500,-2796500,-2796500,3191500,3191500],
    "pism10km":[600,600,-2797500,-2797500,3192500,3192500],
    "pism7km":[800,800,-2798750,-2798750,3193750,3193750],
    "pism5km":[1200,1200,-2800000,-2800000,3195000,3195000], # 5km standard Albmap input grid
    "pism3km":[2000,2000,-2801000,-2801000,3196000,3196000],
    "pism2km":[3000,3000,-2801500,-2801500,3196500,3196500],
    "pism1km":[6000,6000,-2802000,-2802000,3197000,3197000],
     # only corners are relevant for initMIP
    "initmip8km":[761,761,-3040000,-3040000,3040000,3040000],
}

#torsten local and tumble
#output_data_path = os.path.expanduser("/p/projects/tumble/pism_input/GitLab/")
#output_data_path = os.path.expanduser("/home/albrecht/Documents/pism/python/pism_input/")

# matthias
#output_data_path = os.path.expanduser("~/data/20170316_PismInputData/")
#output_data_path = "/p/projects/tumble/mengel/pismInputData/20170316_PismInputData"
output_data_path = os.path.expanduser("/p/projects/pism/mengel/pism_input/")

# RACMO data is not freely available and cannot be downloaded,
# so we have to provide an explicit path here
# if racmo data is intended to be used in publications,
# get in contact with the racmo authors before.
# deprecated: RACMO v2.1 i2s data. not suggested to use in PISM runs
racmo_i2s_data_path = "/p/projects/tumble/mengel/pismSourceData/20170328_RacmoHadCM3_Ice2Sea"

# updated and driven by Reanalysis. Preferred to use:
racmo_wessem_data_path = "/p/projects/tumble/mengel/pismSourceData/20170626_RacmoData_Wessem_etal"


# the till friction angle (tillphi) can be infered from PISM inversion.
# for now, we rely on an inversion run from Torsten on 15km.
# The code here allows to remap the 15km to other resolutions.
tillphi_data_path = "/p/tmp/albrecht/pism17/pismOut/forcing/forcing2300_TPSO/results/result_fit_15km_50000yrs.nc"

# merge the follwing dataset into one PISM-ready file.
# datasets should be named here as the subfolder of its preprocessing.
# TODO: merging schmidtko data does not work yet due to subtle grid differences.
#       it can be supplied through -ocean cavity -ocean_cavity_file $file directly.
datasets_to_merge = ["bedmap2","albmap","racmo_wessem",
                     # "schmidtko",
                     "tillphi_pism"]

# choose here which variables should be taken from which dataset
# The basins variable for the PICO model can be passed with the Schmidtko dataset.
variables = {"bedmap2":["thk","topg","usurf"],
             "albmap":["bheatflx"],
             "racmo_wessem":["smb","air_temp"],
             # "schmidtko":["theta_ocean","salinity_ocean","basins"],
             "tillphi_pism":["tillphi"]}

#### No edits needed below that line. ####
cdo_remapgridpath = os.path.join(output_data_path,"cdo_remapgrids")
project_root = os.path.dirname(os.path.abspath(__file__))
username = pwd.getpwuid(os.getuid()).pw_name
