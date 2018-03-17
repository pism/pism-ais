#!/usr/bin/env python 

########################################################################
## This file downloads the WDC temperature reconstructions from 
## Cuffey et al., 2016, see http://www.usap-dc.org/view/dataset/600377
## by torsten.albrecht@pik-potsdam.de
########################################################################

from numpy import zeros
from netCDF4 import Dataset as NC
import os, sys

project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if project_root not in sys.path: sys.path.append(project_root)
import config as cf; reload(cf)

datain = cf.paleo_time_input
datfile = "WDT_cuffey2016Eq2.txt"
link = "http://www.usap-dc.org/dataset/ldeo/NSF-ANT05-39232/2017-01-12_10-26-36/"+datfile
cmd = 'wget -r '+link+' -O '+datfile

if not os.path.isfile(datain+datfile): 
  print "Please downloading "+datfile+" from "+link+" and save to "+datain
  #os.system(cmd)
  #os.system("mv "+datfile+" "+datain)


# jouzel_07 T ###########################################################################
f = open(datain+datfile)

datalength=39141
datastart=24

time = zeros([datalength-datastart])
temp = zeros([datalength-datastart])


for linecount,line in enumerate(f.readlines()):
  for entrycount,entry in enumerate(line.split('  ')):
    if linecount>=datastart and linecount<datalength:
      if entrycount==1:
        time[linecount-datastart] = float(entry)*1.0e3
      elif entrycount==2:
        temp[linecount-datastart] = float(entry)
f.close()

timeBP = 32
#print temp[timeBP],time[timeBP]
temp -= temp[timeBP]


#########################################################################

timeseries = 'timeseries_cuffey16_temp.nc'
ncf = NC(timeseries,'w', format='NETCDF3_CLASSIC')

# define time dimension, then time variable, then attributes
timedim = ncf.createDimension('time', None)
yearvar = ncf.createVariable('time', 'f4', dimensions=('time',))
setattr(yearvar, 'units', 'years since 1950-01-01')
#os.system("chmod 775 " + timeseries + " 2> /dev/null")

# define variable
temperature=ncf.createVariable('delta_T', 'f4', dimensions=('time',))
setattr(temperature, 'units', 'Kelvin')
setattr(temperature, 'interpolation', 'linear')
temperature.long_name = 'Antarctic Temperature Anomaly'

yearvar[:] = time[:] #+50.0
temperature[:] = temp[:]
# close
ncf.close()



ncpdqcommand="ncpdq -O --rdr=-time "+timeseries+" "+timeseries
# reverse time dimension so that
ncapcommand='ncap2 -O -s "time=-time" '+timeseries+' '+timeseries
# times follow same convention as PISM
ncattedcommand='ncatted -O -a calendar,time,c,c,"365_day" '+timeseries
# set calendar convention as PISM
os.system(ncpdqcommand)
os.system(ncapcommand)
os.system(ncattedcommand)

print " time series written into NetCDF file ",timeseries
