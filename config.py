"""
This file contains user defined paths and settings.
Normally, this should not be committed unless your changes are major.
"""
import os
import pwd

authors="matthias.mengel@pik-potsdam.de and torsten.albrecht@pik-potsdam.de"

proj4str="+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

# grid_id determines also the resolution
#grid_id = "initmip32km"
grid_id = "initmip16km"
#grid_id = "initmip8km"
#grid_id = "initmip4km"
#grid_id = "initmip1km"
#grid_id = "pism7km"

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
    "initmip1km":[6081,6081,-3040000,-3040000,3040000,3040000],
    "initmip2km":[3041,3041,-3040000,-3040000,3040000,3040000],
    "initmip4km":[1521,1521,-3040000,-3040000,3040000,3040000],
    "initmip8km":[761,761,-3040000,-3040000,3040000,3040000],
    "initmip16km":[381,381,-3040000,-3040000,3040000,3040000],
    "initmip32km":[191,191,-3040000,-3040000,3040000,3040000],
    "initmip64km":[95,95,-3040000,-3040000,3040000,3040000]}

# if true, prepare regridding bash grid to be submitted and not run interactively.
cluster_regridding = True
# cdo has a feature to extrapolate data into data-empty regions
use_cdo_extrapolation=True

# regridding_method: choose "bilinear", "integer" or "conservative"
# use conservative remapping to remap fine bed topography to coarser grids,
# works only for bedmap2 and albmap data.
regridding_method = "bilinear"

time_averaging_period = [1986,2005]

#torsten local and tumble
#output_data_path = os.path.expanduser("/p/projects/tumble/pism_input/GitLab/")
#output_data_path = os.path.expanduser("/home/albrecht/Documents/pism/python/pism_input/")
#output_data_path = os.path.expanduser("/p/projects/tumble/albrecht/pism_input/data/")
#output_data_path = os.path.expanduser("/p/tmp/garbe/projects/LARMIP/forcingData/")
#output_data_path = os.path.expanduser("/home/albrecht/Documents/pism/data/pism-ais/")

# matthias
#output_data_path = os.path.expanduser("/p/projects/pism/mengel/pism_input/")

# julius
output_data_path = os.path.expanduser("/p/projects/pism/pism_input_data/Velocity/")

# RACMO data is not freely available and cannot be downloaded,
# so we have to provide an explicit path here
# if racmo data is intended to be used in publications,
# get in contact with the racmo authors before.
# deprecated: RACMO v2.1 i2s data. not suggested to use in PISM runs
racmo_i2s_data_path = "/p/projects/pism/pism_input_data/Racmo_Ice2Sea"

# updated and driven by Reanalysis. Preferred to use:
racmo_wessem_data_path = "/p/projects/pism/pism_input_data/Racmo_ERAInterim"

# driven by HadCM3.
racmo_hadcm3_data_path = "/p/projects/pism/pism_input_data/Racmo_HadCM3"

# driven by CESM2
racmo_cesm2_data_path = "/p/projects/pism/pism_input_data/Racmo_CESM2"


# data is from Keisha Shimada, personal communication, based on
# DOI: 10.1175/JTECH-D-16-0075.1
# please du not use in publications (yet)
shimada_socean_data_path = "/p/projects/tumble/mengel/pismSourceData/20180424_Shimada_SouthernOceanClimatology"

# Bedmachine Antarctica from Mathieu Morlighem, personal communication
# please do not use in publications (yet)
bedmachine_path = "/p/projects/pism/pism_input_data/BedMachine"


# the till friction angle (tillphi) can be infered from PISM inversion.
# for now, we rely on an inversion run from Torsten on 16km.
tillphi_data_path = "/p/projects/tumble/mengel/pismSourceData/20170807_PismTillPhiFromTorsten/result_fit_16km_50000yrs_2411_TPSO.nc"

#paleo_time_input = "/home/albrecht/Documents/pism/data/paleo_timeseries/"
paleo_time_input = "/p/projects/tumble/pism_input/Paleo/"

raised_data_path = "/p/projects/tumble/pism_input/Raised/data/"

# anomaly data used for initMIP experiments
# Downloaded dBasalMelt and dSMB anomaly fields from
# ftp searise@cryoftp1.gsfc.nasa.gov initMIP directory /ISMIP6/initMIP/AIS
initmip_data_path = "/p/projects/tumble/pism_input/ISMIP6/initMIP/AIS/"
# PISM initial (background) state for initMIP experiments.
# Note: this depends on the PISM ice sheet state.
#initmip_pism_out = "/p/tmp/albrecht/pism17/pismOut/forcing/forcing2308_TPSO/results/result_forcing_16km_205000yrs.nc"
#initmip_pism_out = "/p/tmp/albrecht/pism17/pismOut/forcing/forcing2294f_LGM/results/result_forcing_16km_205000yrs.nc"
#initmip_pism_out = "/p/tmp/albrecht/pism18/pismOut/pism_paleo/pism1.0_paleo05_5073/paleo.nc"
#initmip_pism_out = "/p/tmp/albrecht/pism18/pismOut/larmip/larmip2302/results/snap_constant_8km_100000yrs.nc_-50000.000.nc"
#initmip_pism_out = "/p/tmp/garbe/projects/AIS_Equilibrium/pismpik_ant_equi_larmip8km_1.0_10/snapshots_29000.000.nc"
#initmip_pism_out = "/p/tmp/garbe/projects/AIS_Equilibrium/pismpik_ant_equi_larmip8km_1.0_10/snapshots_58000.000.nc"
#initmip_pism_out = "/p/tmp/albrecht/pism18/pismOut/larmip/larmip2304/results/result_constant_8km_100000yrs.nc"
initmip_pism_out = "/p/tmp/albrecht/pism18/pismOut/forcing/forcing2542_TPSO/results/snapshots_2000.000.nc"
#initmip_pism_out = "/p/projects/pism/albrecht/abumip/results/forcing2543e6_TPSO/results/result_ctrl_4km_1000yrs.nc"
initmip_pism_out = "/p/tmp/mengel/pism_out/opttphi_056_initmip4km_nomassoceanconst_amundsenm0p37/snapshots_2450.000.nc"

# merge the follwing dataset into one PISM-ready file.
# datasets should be named here as the subfolder of its preprocessing.
# TODO: merging schmidtko data does not work yet due to subtle grid differences.
#       it can be supplied through -ocean cavity -ocean_cavity_file $file directly.
datasets_to_merge = ["bedmap2",
                     "albmap","racmo_wessem",
                     #"schmidtko",
                     "tillphi_pism"]

# choose here which variables should be taken from which dataset
# The basins variable for the PICO model can be passed with the Schmidtko dataset.
variables = {"bedmap2":["thk","topg","usurf"],
             "albmap":["bheatflx"],
             "racmo_wessem":["climatic_mass_balance","ice_surface_temp"],
             #"schmidtko":["theta_ocean","salinity_ocean","basins"],
             "tillphi_pism":["tillphi"]}

#### No edits needed below that line. ####
cdo_remapgridpath = os.path.join(output_data_path,"cdo_remapgrids")
project_root = os.path.dirname(os.path.abspath(__file__))
username = pwd.getpwuid(os.getuid()).pw_name
