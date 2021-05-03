#!/usr/bin/env python

"""
In case of gaps in RACMO data (missing months), we can apply a workaround and replace the missing month by interpolating between the preceding and subsequent months.

julius.garbe@pik, 2021
"""

import shutil as su
import netCDF4 as nc

#############################

file_path = '/p/projects/pism/garbe/2018_PISM_Input_Data/racmo_cesm2/initmip8km/'
file_in = 'racmo_cesm2_hist_1950_2014_monthly_initmip8km.nc'
file_out = file_in.replace('.nc','_patched.nc')   # this will be the patched file

i_missing = [65,]   # indices of missing months

varnames = ['albedo','air_temp','ice_surface_temp','snowmelt']

#############################

# Duplicate original file first to avoid overwriting it (depending on file size, this step might take a while)
su.copyfile(file_path+file_in, file_path+file_out)

# Read in the data
data = nc.Dataset(file_path+file_out, mode="r+", format="NETCDF4_CLASSIC")

vardata = {}

for var in varnames:

    vardata[var] = data.variables[var]

    for i in i_missing:
        vardata[var][i] = (vardata[var][i-1] + vardata[var][i+1]) / 2.

data.close()

