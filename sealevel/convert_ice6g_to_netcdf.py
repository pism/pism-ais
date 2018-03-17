#!/usr/bin/env python 

########################################################################
## This file downloads the eustatic sea-level reconstructions from 
## Stuhne & Peltier, 2015
## by torsten.albrecht@pik-potsdam.de
########################################################################

from numpy import zeros
from pylab import figure, plot, axis, xlabel, ylabel, show, legend
from netCDF4 import Dataset as NC
import os


datfile = "esl_incl_daocean.txt"


############################################################################
f = open(datfile)

datalength=124
datastart=2

time = zeros([datalength-datastart])
slev = zeros([datalength-datastart])


for linecount,line in enumerate(f.readlines()):
  for entrycount,entry in enumerate(line.split(' ')):
    if linecount>=datastart and linecount<datalength:
      print entrycount,entry
      if entrycount==0:
        time[linecount-datastart] = float(entry)*1.0e3
      elif entrycount==2:
        slev[linecount-datastart] = float(entry)
f.close()


#####################################################################################

# temperature
fig2=figure(2, figsize=(10,5));
ax2 = fig2.add_subplot(111)
ax2.plot(time[:],slev[:], linewidth=2, alpha=0.9, color='b',label="ICE6G:ESL15")
ax2.axis([-35000,0,-140,40])
ax2.legend()
ylabel("dSL (C)")
xlabel("age (y BP)")
show()


#########################################################################

timeseries = 'timeseries_peltier15_sl.nc'
ncf = NC(timeseries,'w', format='NETCDF3_CLASSIC')

# define time dimension, then time variable, then attributes
timedim = ncf.createDimension('time', None)
yearvar = ncf.createVariable('time', 'f4', dimensions=('time',))
setattr(yearvar, 'units', 'years since 1950-01-01')
#os.system("chmod 775 " + timeseries + " 2> /dev/null")

sealevel=ncf.createVariable('delta_SL', 'f4', dimensions=('time',))
setattr(sealevel, 'units', 'meters')
setattr(sealevel, 'interpolation', 'linear')
sealevel.long_name = 'Relative Sea Level (variation from present)'
setattr(sealevel, 'standard_name', 'global_average_sea_level_change')

yearvar[:] = time[:] #+50.0
sealevel[:] = slev[:]
# close
ncf.close()


ncattedcommand='ncatted -O -a calendar,time,c,c,"365_day" '+timeseries
# set calendar convention as PISM
os.system(ncattedcommand)

print " time series written into NetCDF file ",timeseries
