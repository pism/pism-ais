#!/usr/bin/env python 

########################################################################
## This file uses the EDC temperature reconstructions from 
## Jouzel et al., 2007, to create a moving average ocean forcing
## by ricarda.winkelmann@pik-potsdam.de, julius.garbe@pik-potsdam.de
## and torsten.albrecht@pik-potsdam.de, 
########################################################################

from numpy import zeros
from netCDF4 import Dataset as NC
import os, sys
import numpy as np

project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if project_root not in sys.path: sys.path.append(project_root)
import config as cf; reload(cf)


""" GET FORCING TIMESERIES """
surface_temp_file="timeseries_jouzel07_temp.nc"

ncin      = NC(surface_temp_file, 'r')
DeltaT_in = ncin.variables['delta_T'][:]
years_in  = ncin.variables['time'][:]#-50.0
ncin.close()

print "EDC timeseries has "+str(len(years_in))+" points"


ratio_surface = 1.8 # ratio of Antarctic surface warming to global mean temperature increase
ratio_ocean   = 0.75 # ratio of Antarctic ocean warming at shelf depth to global mean temperature increase

tbin = 500.0
print "Running mean with time bins of "+str(tbin)+" years"
time_bin = np.arange(tbin*np.int(years_in[0]/tbin),2.0*tbin,tbin)
DeltaT_bin=np.zeros(len(time_bin))
count_bin = np.zeros_like(DeltaT_bin)

for i,year in enumerate(years_in):
  indxav, = np.where( time_bin == tbin*np.int(year/tbin))[0]
  #print i,year,indxav,time_bin[indxav]
  DeltaT_bin[indxav]+=DeltaT_in[i]*ratio_ocean
  count_bin[indxav]+=1.0

DeltaT_bin/=count_bin
DeltaT_bin[-1]=-DeltaT_bin[-2]
#print DeltaT_bin,count_bin

""" VARIABLES FOR NC-FILE """
tvar = time_bin # time vector for output file
dTvar = DeltaT_bin # warming vector for output file

"""
define dimensions in NetCDF file
"""
outputFile="timeseries_jouzel07_tempocean_mean.nc"
ncout = NC(outputFile, 'w',format='NETCDF3_CLASSIC')
timedim = ncout.createDimension('time', None)

    
"""
define variables, set attributes, write data
"""
    
nct  = ncout.createVariable( 'time','f4', dimensions=('time',) )
setattr(nct, 'units', 'years since 1-1-1')
ncdT   = ncout.createVariable( "delta_T",'f4', dimensions=('time',) )
    
nct[:] = tvar
ncdT[:] = dTvar
    
"""
attributes in NetCDF file
"""
setattr(ncout, 'Conventions', 'CF-1.3') # only global attribute
setattr(nct, 'units', 'years since 1-1-1')
setattr(ncdT, 'long_name', 'Ocean Temperature (variation from present)')
setattr(ncdT, 'standard_name', 'ocean_temperature')
setattr(ncdT, 'units', 'Kelvin')
setattr(ncdT, 'interpolation', 'linear')   
ncout.close()


ncattedcommand='ncatted -O -a calendar,time,c,c,"365_day" '+outputFile
# set calendar convention as PISM
os.system(ncattedcommand)

print "NetCDF file ",outputFile," created. \n"
print "Done. \n"