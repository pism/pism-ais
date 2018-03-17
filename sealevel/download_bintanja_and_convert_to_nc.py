#!/usr/bin/env python 

########################################################################
## This file downloads the global sea-level reconstructions 
## from Bintanja et al., 2008
## by torsten.albrecht@pik-potsdam.de
########################################################################

from numpy import zeros
from netCDF4 import Dataset as NC
import os, sys

project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if project_root not in sys.path: sys.path.append(project_root)
import config as cf; reload(cf)

datain = cf.paleo_time_input
datfile = "bintanja2008.txt"
link = "ftp://ftp.ncdc.noaa.gov/pub/data/paleo/contributions_by_author/bintanja2008/"+datfile
cmd = 'wget -r '+link+' -O '+datfile
if not os.path.isfile(datain+datfile): 
  print "Downloading "+datfile
  os.system(cmd)
  os.system("mv "+datfile+" "+datain)

############################################################################
f = open(datain+datfile)

datalength=30110
datastart=109

time = zeros([datalength-datastart])
slev = zeros([datalength-datastart])


for linecount,line in enumerate(f.readlines()):
  for entrycount,entry in enumerate(line.split('    ')):
    if linecount>=datastart and linecount<datalength:
      #print entrycount,entry
      if entrycount==0:
        time[linecount-datastart] = float(entry)*(-1.0e3)
      elif entrycount==8:
        slev[linecount-datastart] = float(entry)*(-1.0)
f.close()

#########################################################################

timeseries = 'timeseries_bintanja08_sl.nc'
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
