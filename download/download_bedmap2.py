#!/usr/bin/env python

# This file is part of bedmap2pism.
# Copyright (C) 2015-2016
authors="matthias.mengel@pik-potsdam.de and torsten.albrecht@pik-potsdam.de"

# downloads BEDMAP2
# creates 1km dataset NetCDF files and adjusts metadata,
# depends on wget and NCO (ncrename, ncap2, ncatted, ncpdq, ncks)



# bedmap2pism is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# bedmap2pism is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with bedmap2pism.  If not, see <http://www.gnu.org/licenses/>.

import sys, os, imp
import datetime
import netCDF4 as nc
import numpy as np
from PISMNC import PISMDataset as PNC
from mpl_toolkits.basemap import interp

#import download_data; reload(download_data)

### Bedmap2   ##########################################################
# Documentation of the data: https://www.bas.ac.uk/project/bedmap-2/
bedmap2_bin="bedmap2_bin"
bedmap2_link="https://secure.antarctica.ac.uk/data/bedmap2/"+bedmap2_bin+".zip"
bedmap2_data="bedmap2_data"
bedmap2_nc= 'bedmap2_1km.nc'
ncbm2_name=bedmap2_data+"/"+bedmap2_nc
N=6667
x0_bm2= 3333500.
y0_bm2=x0_bm2
if not os.path.exists(bedmap2_bin):
  print "Downloading "+bedmap2_bin+"\n"
  os.system("wget " + bedmap2_link)
  os.system("unzip " + bedmap2_bin+".zip")
  os.system("rm " + bedmap2_bin+".zip")


if not os.path.isfile(bedmap2_data+"/"+bedmap2_nc):
  os.system("mkdir "+bedmap2_data)
  print "Reading bedmap2 binary files from %s ...\n" % (bedmap2_bin)

  fname = bedmap2_bin + '/bedmap2_bed.flt'
  bed = np.flipud(np.ma.masked_equal(np.reshape(np.fromfile(fname,dtype=np.float32),(N,N)),-9999.0))
  fname = bedmap2_bin + '/bedmap2_thickness.flt'
  thk = np.flipud(np.ma.masked_equal(np.reshape(np.fromfile(fname,dtype=np.float32),(N,N)),-9999.0))
  fname = bedmap2_bin + '/bedmap2_icemask_grounded_and_shelves.flt'
  mask = np.flipud(np.ma.masked_equal(np.reshape(np.fromfile(fname,dtype=np.float32),(N,N)),-9999.0))
  fname = bedmap2_bin + '/bedmap2_grounded_bed_uncertainty.flt'
  bedunc = np.flipud(np.ma.masked_equal(np.reshape(np.fromfile(fname,dtype=np.float32),(N,N)),-9999.0))
  fname = bedmap2_bin + '/bedmap2_surface.flt'
  usurf = np.flipud(np.ma.masked_equal(np.reshape(np.fromfile(fname,dtype=np.float32),(N,N)),-9999.0))

  print " range of bed = [%.2f, %.2f]" % (bed.min(),bed.max())
  print " range of thk = [%.2f, %.2f]" % (thk.min(),thk.max())
  print " range of mask = [%.2f, %.2f]" % (mask.min(),mask.max())
  print " range of bedunc = [%.2f, %.2f]" % (bedunc.min(),bedunc.max())
  print " range of usurf = [%.2f, %.2f]" % (usurf.min(),usurf.max())

  #get rid off NaN
  thk[thk.mask]=0.0
  bed[bed.mask]=-5000.0
  usurf[usurf.mask]=0.0
  bedunc[bedunc.mask]=0.0
  mask[mask.mask]=2.0

  #########################################################################

  print "writing NetCDF file '%s' ..." % ncbm2_name
  try:
      nc = PNC(ncbm2_name, 'w', format='NETCDF4_CLASSIC')
  except:
      print("can't open file %s for writing" % ncbm2_name)
      exit(1)

  print " writing x,y ..."
  dx = 1000.0
  dy = 1000.0
  x = np.linspace(-(N-1)*dx/2.0,(N-1)*dx/2.0,N)
  y = np.linspace(-(N-1)*dy/2.0,(N-1)*dy/2.0,N)
  nc.create_dimensions(x, y, time_dependent = False)



  print " writing topg ..."
  nc.define_2d_field("topg", time_dependent = False,
                     attrs = {"long_name" : "elevation of bedrock",
                              "valid_range" : (-9000.0, 9000.0),
                              "standard_name" : "bedrock_altitude",
                              "units" : "meters"})
  nc.write_2d_field("topg", bed)

  print " writing usurf ..."
  nc.define_2d_field("usurf", time_dependent = False,
                     attrs = {"long_name" : "ice upper surface elevation",
                              "valid_range" : (-1000.0, 9000.0),
                              "standard_name" : "surface_altitude",
                              "units" : "meters"})
  nc.write_2d_field("usurf", usurf)

  print " writing thk ..."
  nc.define_2d_field("thk", time_dependent = False,
                     attrs = {"long_name" : "thickness of ice sheet or ice shelf",
                              "valid_range" : (0.0, 9000.0),
                              "standard_name" : "land_ice_thickness",
                              "units" : "meters"})
  nc.write_2d_field("thk", thk)

  print " writing bedunc ..."
  nc.define_2d_field("bedunc", time_dependent = False,
                     attrs = {"long_name" : "uncertainty of bed topography",
                              "valid_range" : (0.0, 9000.0),
                              "standard_name" : "bed_uncertainty",
                              "units" : "meters"})
  nc.write_2d_field("bedunc", bedunc)

  print " writing mask ..."
  nc.define_2d_field("mask", time_dependent = False,
                     attrs = {"long_name" : "ice-type (ice-free/grounded/floating/ocean) integer mask",
                              "valid_range" : (0.0, 2.0),
                              "standard_name" : "mask",
                              "units" : ""})
  nc.write_2d_field("mask", mask)


  now = datetime.datetime.now().strftime("%B %d, %Y")
  #nc.projection = "+proj=stere +lon_0=0 +lat_0=-90 +lat_ts=-71 +ellps=WGS84 +datum=WGS84"
  nc.proj4 = "+proj=stere +lat_0=90 +lat_ts=71 +lon_0=-39 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
  nc.Comment  = "albrecht@pik-potsdam.de created netcdf bedmap2 file at " + now


  nc.close()
  print "done"
