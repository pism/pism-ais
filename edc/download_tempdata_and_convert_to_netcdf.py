#!/usr/bin/env python 

########################################################################
## This file downloads the EDC temperature reconstructions from 
## Jouzel et al., 2007, by torsten.albrecht@pik-potsdam.de
########################################################################

from numpy import zeros #, squeeze
from pylab import figure, plot, axis, xlabel, ylabel, show, legend
from netCDF4 import Dataset as NC
import os


datfile = "edc3deuttemp2007.txt"
link = "ftp://ftp.ncdc.noaa.gov/pub/data/paleo/icecore/antarctica/epica_domec/"+datfile
cmd = 'wget -r '+link+' -O '+datfile
os.system(cmd)


# jouzel_07 T ###########################################################################
f = open(datfile)

datalength=5892
datastart=104
dataexcept=[170,219,536]

time = zeros([datalength-datastart])
temp = zeros([datalength-datastart])

count=0

for linecount,line in enumerate(f.readlines()):
  for entrycount,entry in enumerate(line.split('      ')):
    if linecount>=datastart and linecount<datalength:
      #print entrycount,entry
      if entrycount==0:
        if float(entry) in dataexcept:
          count+=1
          break
      elif entrycount==2:
        time[linecount-datastart-count] = float(entry)
      elif entrycount==4:
        temp[linecount-datastart-count] = float(entry)
f.close()


#####################################################################################

# temperature
fig2=figure(2, figsize=(10,5));
ax2 = fig2.add_subplot(111)
ax2.plot(time[3:-3],temp[3:-3], linewidth=2, alpha=0.9, color='b',label="EDC:jouzel07")
ax2.axis([0,35000,-12,8])
ax2.legend()
ylabel("dT_AA (C)")
xlabel("age (y BP, EDC3)")
show()


#########################################################################

timeseries = 'timeseries_jouzel07_temp.nc'
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
