#!/usr/bin/env python

import sys, os

project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if project_root not in sys.path: sys.path.append(project_root)
import config as cf; reload(cf)

albmap_nc="Antarctica_5km_dev1.0.nc"
albmap_link="http://websrv.cs.umt.edu/isis/images/4/4d/"+albmap_nc

dataset = "albmap"
data_path = os.path.join(cf.output_data_path, dataset)
albmap_orig = os.path.join(data_path,albmap_nc)

if not os.path.isfile(albmap_orig):
  print "Downloading "+albmap_nc
  os.system("wget -nc "+albmap_link)
  os.system("mv "+albmap_nc+' '+albmap_orig)


outname="timeseries_specmap_sl.nc"
if not os.path.isfile(outname): 
  print "Extracting "+outname
  os.system("ncks -A -v sealevel_time_series,sealeveltimes "+albmap_orig+' '+outname)
  os.system("ncpdq -O --rdr=-sealeveltimes "+outname+' '+outname)
  os.system("ncap2 -O -s 'sealeveltimes=-sealeveltimes' "+outname+' '+outname)
  os.system("ncatted -O -a calendar,sealeveltimes,c,c,365_day "+outname)
  os.system("ncrename -v sealevel_time_series,delta_SL "+outname)
  os.system("ncrename -d sealeveltimes,time -v sealeveltimes,time "+outname)

