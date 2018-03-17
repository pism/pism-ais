#!/usr/bin/env python 

########################################################################
## This file downloads the WDC accumulation reconstructions from 
## Fudge et al., 2016, see http://www.usap-dc.org/view/dataset/601004
## by torsten.albrecht@pik-potsdam.de
########################################################################

from numpy import zeros, concatenate
from pylab import figure, plot, axis, xlabel, ylabel, show, legend, axhline
from netCDF4 import Dataset as NC
import os


datfile = "WDC_accumulation_combined_movavg.csv"
#link = "http://www.usap-dc.org/dataset/usap-dc/601004/2018-01-16T22:27:48.0Z/"+datfile
#cmd = 'wget -r '+link+' -O '+datfile
#os.system(cmd)


# fudge 16 P ###########################################################################
f = open(datfile)

datalength=1358
datastart=3

time = zeros([datalength-datastart])
accum = zeros([datalength-datastart])

for linecount,line in enumerate(f.readlines()):
  for entrycount,entry in enumerate(line.split(' ')):
    if linecount>=datastart and linecount<datalength:
      #print entrycount,entry
      entry = entry.replace(',','.')
      if entrycount==0:
        time[linecount-datastart] = float(entry)
      elif entrycount==1:
        accum[linecount-datastart] = float(entry)
f.close()

#print time
timeBP = 0
print accum[timeBP],time[timeBP]

frac = accum/accum[timeBP]
anomaly = accum-accum[timeBP]



#####################################################################################

# accum
fig2=figure(2, figsize=(10,5));
ax2 = fig2.add_subplot(111)
ax2.axhline(1,linestyle='dashed')
#ax2.plot(time,accum, linewidth=2, alpha=0.9, color='b',label="WDC:fudge16")
#ax2.plot(time,anomaly, linewidth=2, alpha=0.9, color='b',label="WDC:fudge16")
#ax2.plot(time2,accum2, linewidth=2, alpha=0.9, color='g',label="WDC:fudge16")
ax2.plot(time,frac, linewidth=2, alpha=0.9, color='g',label="WDC:fudge16")
#ax2.axis([0,35000,0.,0.5])
#ax2.axis([0,35000,-1,1])
ax2.legend()
ylabel("P (m ice a-1)")
xlabel("age (y BP)")
show()


#########################################################################

timeseries = 'timeseries_fudge16_accum_movavg50.nc'
ncf = NC(timeseries,'w', format='NETCDF3_CLASSIC')

# define time dimension, then time variable, then attributes
timedim = ncf.createDimension('time', None)
yearvar = ncf.createVariable('time', 'f4', dimensions=('time',))
setattr(yearvar, 'units', 'years since 1950-01-01')
#os.system("chmod 775 " + timeseries + " 2> /dev/null")

# define variable
accumulation=ncf.createVariable('delta_P', 'f4', dimensions=('time',))
#setattr(accumulation, 'units', 'm year-1')
setattr(accumulation, 'units', 'kg m-2 year-1')
setattr(accumulation, 'interpolation', 'linear')
accumulation.long_name = 'WDC Precipitation Anomaly'


accumfrac=ncf.createVariable('frac_P', 'f4', dimensions=('time',))
setattr(accumfrac, 'units', '1')
setattr(accumfrac, 'interpolation', 'linear')
accumfrac.long_name = 'WDC Precipitation Multiplier'


yearvar[:] = time[:] #+50.0
accumulation[:] = anomaly[:]*910.0 #kg/m2
accumfrac[:] = frac[:]


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
