#!/usr/bin/env python

import sys, os


albmap_nc="Antarctica_5km_dev1.0.nc"
albmap_link="http://websrv.cs.umt.edu/isis/images/4/4d/"+albmap_nc
if not os.path.isfile(albmap_nc): 
  print "Downloading "+albmap_nc
  os.system("wget -nc "+albmap_link)


outname="imbrie_mcintyre06_sl.nc"
if not os.path.isfile(outname): 
  print "Extracting "+outname
os.system("ncks -A -v sealevel_time_series,sealeveltimes "+albmap_nc+' '+outname)
os.system("ncpdq -O --rdr=-sealeveltimes "+outname+' '+outname)
os.system("ncap2 -O -s 'sealeveltimes=-sealeveltimes' "+outname+' '+outname)
os.system("ncatted -O -a calendar,sealeveltimes,c,c,365_day "+outname)
os.system("ncrename -v sealevel_time_series,delta_SL "+outname)
os.system("ncrename -d sealeveltimes,time -v sealeveltimes,time "+outname)

