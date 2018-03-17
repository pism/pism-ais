#!/usr/bin/env python 

########################################################################
## This file downloads the WDC temperature reconstructions from 
## Cuffey et al., 2016, see http://www.usap-dc.org/view/dataset/600377
## by torsten.albrecht@pik-potsdam.de
########################################################################

from numpy import zeros
from pylab import figure, plot, axis, xlabel, ylabel, show, legend
from netCDF4 import Dataset as NC
import os


datfile = "WDT_cuffey2016Eq2.txt"
link = "http://www.usap-dc.org/dataset/ldeo/NSF-ANT05-39232/2017-01-12_10-26-36/"+datfile
cmd = 'wget -r '+link+' -O '+datfile
#os.system(cmd)


# jouzel_07 T ###########################################################################
f = open(datfile)

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
print temp[timeBP],time[timeBP]
temp -= temp[timeBP]


#####################################################################################

# temperature
fig2=figure(2, figsize=(10,5));
ax2 = fig2.add_subplot(111)
ax2.plot(time[3:-3],temp[3:-3], linewidth=2, alpha=0.9, color='b',label="WDC:cuffey16")
ax2.axis([0,35000,-12,8])
ax2.legend()
ylabel("dT_AA (C)")
xlabel("age (y BP, EDC3)")
show()


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

yearvar[:] = time[0:-3] #+50.0
temperature[:] = temp[0:-3]
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
