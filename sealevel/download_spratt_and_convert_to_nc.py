#!/usr/bin/env python 

########################################################################
## This file downloads the global sea-level reconstructions 
## from Spratt & Lisiecki, 2016
## by torsten.albrecht@pik-potsdam.de
########################################################################

from numpy import zeros
from netCDF4 import Dataset as NC
import os, sys

project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if project_root not in sys.path: sys.path.append(project_root)
import config as cf; reload(cf)

datain = cf.paleo_time_input
datfile = "spratt2016.txt"
link = "https://www1.ncdc.noaa.gov/pub/data/paleo/contributions_by_author/spratt2016/"+datfile
cmd = 'wget '+link+' -O '+datfile
if not os.path.isfile(datain+datfile): 
  print "Downloading "+datfile
  os.system(cmd)
  os.system("mv "+datfile+" "+datain)

############################################################################
f = open(datain+datfile)
print datain+datfile

datalength=527
datastart=97

time = zeros([datalength-datastart])
slev = zeros([datalength-datastart])


for linecount,line in enumerate(f.readlines()):
  for entrycount,entry in enumerate(line.split("\t")):
    if linecount>=datastart and linecount<datalength:
      if entrycount==0:
        time[linecount-datastart] = float(entry)*(-1.0e3)
      elif entrycount==1:
        slev[linecount-datastart] = float(entry)
f.close()
#
#########################################################################

timeseries = 'timeseries_spratt16_sl.nc'
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
sealevel[:] = slev[:]-slev[0]
# close
ncf.close()


ncattedcommand='ncatted -O -a calendar,time,c,c,"365_day" '+timeseries
# set calendar convention as PISM
os.system(ncattedcommand)

print " time series written into NetCDF file ",timeseries
