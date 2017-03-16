#!/usr/bin/env python

# This file is part of bedmap2pism.
# Copyright (C) 2015-2016 
authors="matthias.mengel@pik-potsdam.de and torsten.albrecht@pik-potsdam.de"

# downloads ALBMAP, BEDMAP2, Rignot velocity and Arthern accumulation
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

import download_data; reload(download_data)


compare_bm2_alb = True

### Albmap   ##########################################################
albmap_nc="Antarctica_5km_dev1.0.nc"
albmap_link="http://websrv.cs.umt.edu/isis/images/4/4d/"+albmap_nc
albmap_data="albmap_data"
ncalb_name=albmap_data+"/albmap_pism_5km.nc"
if not os.path.isfile(albmap_data+"/"+albmap_nc): 
  print "Downloading "+albmap_nc
  os.system("wget -nc "+albmap_link)
  os.system("mkdir "+albmap_data)
  os.system("mv "+albmap_nc+" "+albmap_data+"/")

### Rignot Velocity   ##########################################################
rignot_bin="rignot_bin"
rignot_nc="rignot_velocity_1km.nc"
rignot_link="https://secure.antarctica.ac.uk/data/bedmap2/resources/Rignot_velocity/bin.zip"
rignot_data="rignot_data"
N2 = 5602
x0_rig=2800500.
y0_rig=2801500.
if not os.path.exists(rignot_bin): 
  print "Downloading "+rignot_data
  os.system("wget "+rignot_link)
  os.system("mv bin.zip "+rignot_bin+".zip")
  os.system("unzip " + rignot_bin+".zip -d "+rignot_bin)
  os.system("rm " + rignot_bin+".zip")


### Arthern accumulation   ##########################################################
arthern_bin="arthern_bin"
arthern_nc="arthern_accumulation_1km.nc"
arthern_link="https://secure.antarctica.ac.uk/data/bedmap2/resources/Arthern_accumulation/Arthern_accumulation_bin.zip"
arthern_data="arthern_data"
N3a=7899
N3b=8300
x0_art=3949500.
y0_art=x0_art
if not os.path.exists(arthern_bin): 
  print "Downloading "+arthern_link.split("/")[-1]
  os.system("wget "+arthern_link)
  os.system("mv Arthern_accumulation_bin.zip "+arthern_bin+".zip")
  os.system("unzip " + arthern_bin+".zip -d "+arthern_bin)
  os.system("rm " + arthern_bin+".zip")

### Bedmap2   ##########################################################
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


  ### add Rignot velocity data   ##########################################
  print "\nReading Rignot binary file from %s ...\n" % (rignot_bin)
  x0=(x0_bm2-x0_rig)/1000.
  y0=(y0_bm2-y0_rig)/1000.
  fname = rignot_bin + '/rignot_velocity_bedmap2_grid.flt'
  vel = np.flipud(np.ma.masked_equal(np.reshape(np.fromfile(fname,dtype=np.float32),(N2,N2)),-9999.0))
  print " range of vel = [%.2f, %.2f]" % (vel.min(),vel.max())
  vel_bm2=np.ones([N,N])*(-9999.0)
  vel_bm2[x0:x0+N2,y0:y0+N2]=vel
  vel_bm2=np.ma.masked_array(vel_bm2, mask=(vel_bm2==-9999.))


  ### add Arthern accumulation data   ##########################################
  print "\nReading Arthern binary file from %s ...\n" % (arthern_bin)
  x0=(x0_art-x0_bm2)/1000.
  y0=(y0_art-y0_bm2)/1000.
  fname = arthern_bin + '/arthern_accumulation_bedmap2_grid.flt'
  accum = np.flipud(np.ma.masked_equal(np.reshape(np.fromfile(fname,dtype=np.float32),(N3b,N3a)),-9999.0))
  #fname = arthern_bin + '/arthern_accumulation_rms_bedmap2_grid.flt'
  print " range of accum = [%.2f, %.2f]" % (accum.min(),accum.max())
  accum_bm2=np.zeros([N,N])
  accum_bm2=accum[x0:x0+N,y0:y0+N]
  accum_bm2=np.ma.masked_array(accum_bm2, mask=(accum_bm2==-9999.))



  ### Write nc-file
  print "\nWriting NetCDF file '%s' ...\n" % ncbm2_name
  try:
    #nc = PNC(ncbm2_name, 'w', format='NETCDF3_CLASSIC')
    nc = PNC(ncbm2_name, 'w', format='NETCDF4_CLASSIC')
  except:
    print("can't open file %s for writing" % ncbm2_name)
    exit(1)


  dx = 1000.0 #m
  dy = 1000.0 #m
  x = np.linspace(0.0,(N-1)*dx,N)
  y = np.linspace(0.0,(N-1)*dy,N)
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
                            "valid_range" : (0.0, 1.0),
                            "standard_name" : "mask",
                            "units" : ""})
  nc.write_2d_field("mask", mask)

  print " writing velocity ..."
  nc.define_2d_field("velocity", time_dependent = False,
                   attrs = {"long_name" : "observed surface velocity after Rignot et al.",
                            "valid_range" : (0.0, 9999.),
                            "standard_name" : "velocity",
                            "units" : "m/year"})
  nc.write_2d_field("velocity", vel_bm2)

  print " writing accumulation ..."
  nc.define_2d_field("accum", time_dependent = False,
                   attrs = {"long_name" : "inverted accumulation after Arthern et al.",
                            "valid_range" : (-9999., 9999.),
                            "standard_name" : "accumulation",
                            "units" : "kg/m2/year"})
  nc.write_2d_field("accum", accum_bm2)


  now = datetime.datetime.now().strftime("%B %d, %Y")
  nc.proj4      = "+proj=stere +lon_0=0 +lat_0=-90 +lat_ts=-71 +ellps=WGS84 +datum=WGS84"
  nc.projection = "+proj=stere +lon_0=0 +lat_0=-90 +lat_ts=-71 +ellps=WGS84 +datum=WGS84"
  nc.Comment  = authors+" created netcdf bedmap2 file at " + now

  nc.close()
  print "\nDone"



### albmap preprocessing analog to pism searise experiment

if not os.path.isfile(ncalb_name):
  print "\nWriting NetCDF file '%s' ...\n" % ncalb_name

  os.system("cp "+albmap_data+"/"+albmap_nc+" "+ncalb_name)
  # following use NCO (http://nco.sourceforge.net/)
  # rename dimensions
  os.system("ncrename -O -v x1,x -v y1,y -d x1,x -d y1,y "+ncalb_name)
  os.system("ncrename -O -v time,t -d time,t "+ncalb_name)
  # fix polar stereographic parameter
  os.system("ncatted -O -a standard_parallel,mapping,m,d,-71.0 "+ncalb_name)
  # rename usurf for convenience
  #os.system("ncrename -O -v usrf,usurf "+ncalb_name)
  # fix surface temperature name and make K
  os.system("ncap2 -O -s 'air_temp=temp+273.15' "+ncalb_name+" "+ncalb_name)
  os.system("ncatted -O -a units,air_temp,m,c,'K' "+ncalb_name)
  # choose Van de Berg et al version of accumulation; will treat as ice-equivalent snow rate
  os.system("ncrename -O -v accr,precipitation "+ncalb_name)
  os.system("ncatted -O -a units,precipitation,m,c,'m/year' "+ncalb_name)
  # use bheatflx_shapiro as the default bheatflx data and 
  os.system("ncrename -O -v bheatflx_shapiro,bheatflx "+ncalb_name)
  os.system("ncatted -O -a units,bheatflx,m,c,'W m-2' "+ncalb_name)
  # delete incorrect standard_name attribute from bheatflx; there is no known standard_name
  os.system("ncatted -O -a standard_name,bheatflx,d,, "+ncalb_name)
  # keep only the fields we actually use at bootstrapping
  os.system("ncks -O -v x,y,lat,lon,bheatflx,topg,thk,precipitation,air_temp,mapping "+ncalb_name+" "+ncalb_name)

  print "\nDone"

#################################################################################