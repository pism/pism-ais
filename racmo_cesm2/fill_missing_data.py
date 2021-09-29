#!/usr/bin/env python

"""
In case of gaps in the RACMO data (missing timeslices in some variables), we can apply a workaround 
and fill the missing data by interpolating between the preceding and subsequent timeslices.

julius.garbe@pik
"""

import os, sys
import shutil
import netCDF4 as nc
## this hack is needed to import config.py from the project root
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if project_root not in sys.path: sys.path.append(project_root)
import config as cf; reload(cf)

dataset = "racmo_cesm2"
period = "hist" # historical: "hist" / SSP5-8.5: "ssp585"
period_years = "1950_2014" # hist: "1950_2014" / ssp585: "2015_2099"
timestep = "monthly" # time step of input data ("monthly" or "yearly")
grid = cf.grid_id

data_path = os.path.join(cf.output_data_path, dataset)
inputfile = os.path.join(data_path, grid, dataset+"_"+period+"_"+period_years+"_"+timestep+"_"+grid+".nc")
outputfile = inputfile.replace('.nc','_patched.nc')   # this will be the patched file

i_missing = [65,]   # indices of missing timeslices

varnames = ['albedo','air_temp','snowmelt']    # variables showing the data gaps

#############################

# Duplicate original file first to avoid overwriting it
if not os.path.exists(outputfile):
    print "Copy original file... (depending on file size, this might take a while)"
    shutil.copyfile(inputfile, outputfile)

# Read in the data
data = nc.Dataset(outputfile, mode="r+", format="NETCDF4_CLASSIC")

vardata = {}

print "Patch variables..."

for var in varnames:
    print var

    vardata[var] = data.variables[var]

    for i in i_missing:
        vardata[var][i] = (vardata[var][i-1] + vardata[var][i+1]) / 2.

data.close()

print "Done."
