#!/usr/bin/env python

# This file is part of bedmap2pism.
# Copyright (C) 2015-2016 
authors="matthias.mengel@pik-potsdam.de and torsten.albrecht@pik-potsdam.de"

# creates combined 5km dataset NetCDF file from Bedmap2 and Albmap data
# adds metadata, ready for PISM

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

#import sys, os, imp
import os
import datetime
import netCDF4 as nc
import numpy as np
from PISMNC import PISMDataset as PNC
from mpl_toolkits.basemap import interp

import download_data; reload(download_data)

ncalb_bm2_name = "albmap_bedmap2_5km.nc"
ncalb_name="albmap_data/albmap_pism_5km.nc"
ncbm2_name="bedmap2_data/bedmap2_1km.nc"
bedmap2_link="https://secure.antarctica.ac.uk/data/bedmap2/bedmap2_bin.zip"
rignot_link="https://secure.antarctica.ac.uk/data/bedmap2/resources/Rignot_velocity/bin.zip"
arthern_link="https://secure.antarctica.ac.uk/data/bedmap2/resources/Arthern_accumulation/Arthern_accumulation_bin.zip"
x0_bm2= 3333500.

compare_bm2_alb = True

#############################################################################################################
### Combine Albmap and Bedmap2
print "Making PISM-readable file combining Albmap and Bedmap2 on 5km resolution"

os.system("cp " + ncalb_name +" "+ ncalb_bm2_name)

ncbm2 = nc.Dataset(ncbm2_name, 'r')
ncalb = nc.Dataset(ncalb_bm2_name, 'a')

# get mask, topo, thk from bedmap 2
xbm2 = ncbm2.variables['x'][:]
ybm2 = ncbm2.variables['y'][:]
topg = ncbm2.variables['topg'][:]
thk  = ncbm2.variables['thk'][:]
mask  = ncbm2.variables['mask'][:]
usurf  = ncbm2.variables['usurf'][:]

# get surface rignot velocity from bedmap 2
vel = ncbm2.variables['velocity'][:]

# get arthern accumulation from bedmap 2
accumulation = ncbm2.variables['accum'][:]

# get others from albmap
xalb = ncalb.variables['x'][:]
yalb = ncalb.variables['y'][:]
#lat  = ncalb.variables['topg'][:]
precip = ncalb.variables['precipitation'][:]
artm   = ncalb.variables['air_temp'][:]

xgrid, ygrid = np.meshgrid(xalb,yalb)
# adjust bedm2 to centered x,y, see bedmap2 readme file
xbm2 -= x0_bm2
ybm2 -= x0_bm2

thkbm2 = np.asarray((interp(thk, xbm2, ybm2, xgrid, ygrid )))
topgbm2 = np.asarray((interp(topg, xbm2, ybm2, xgrid, ygrid )))
maskbm2 = np.asarray((interp(mask, xbm2, ybm2, xgrid, ygrid )))
usurfbm2 = np.asarray((interp(usurf, xbm2, ybm2, xgrid, ygrid )))
thkbm2[thkbm2 > 10000.] = 0.
topgbm2[topgbm2 > 10000.] = -9999
velbm2 = np.asarray((interp(vel, xbm2, ybm2, xgrid, ygrid )))
accumbm2 = np.asarray((interp(accumulation, xbm2, ybm2, xgrid, ygrid )))

### add some difference field to compare albmap and bedma2 topg and thk
if compare_bm2_alb:

  ncmsk  = ncalb.createVariable( 'mask','float32',('t','y','x') )
  ncusurf  = ncalb.createVariable( 'usurf','float32',('t','y','x') )
  ncthkold = ncalb.createVariable( 'thk_alb','float32',('t','y','x') )
  nctopgold = ncalb.createVariable( 'topg_alb','float32',('t','y','x') )
  nctopgdiff = ncalb.createVariable( 'topg_diff','float32',('t','y','x') )
  ncthkdiff = ncalb.createVariable( 'thk_diff','float32',('t','y','x') )

  ncthkold[:] = ncalb.variables['thk'][:]
  ncthkold.units  = ncalb.variables['thk'].units
  ncthkold.standard_name = ncalb.variables['thk'].standard_name
  ncthkold.long_name = ncalb.variables['thk'].long_name

  ncthkdiff[:] = thkbm2- ncalb.variables['thk'][:]
  ncthkdiff.units  = ncalb.variables['thk'].units
  ncthkdiff.standard_name = "bedm2_alb_thk"
  ncthkdiff.long_name = "bedmap2 albmap thickness difference"

  nctopgold[:] = ncalb.variables['topg'][:]
  nctopgold.units  = ncalb.variables['topg'].units
  nctopgold.standard_name = ncalb.variables['topg'].standard_name
  nctopgold.long_name = ncalb.variables['topg'].long_name

  nctopgdiff[:] = topgbm2 - ncalb.variables['topg'][:]
  nctopgdiff.units  = ncalb.variables['topg'].units
  nctopgdiff.standard_name = "bedm2_alb_topg"
  nctopgdiff.long_name = "bedmap2 albmap topograpy difference"

  ncmsk[:] = maskbm2
  ncmsk.units         =  ncbm2.variables['mask'].units
  ncmsk.standard_name =  ncbm2.variables['mask'].standard_name
  ncmsk.long_name     =  ncbm2.variables['mask'].long_name
  ncmsk.valid_range   =  ncbm2.variables['mask'].valid_range

ncalb.variables['thk'][:] = thkbm2
ncalb.variables['topg'][:] = topgbm2
ncalb.variables['usurf'][:] = usurfbm2
ncalb.variables['thk'].valid_range = [0.,9999.]
ncalb.variables['topg'].valid_range = [-9999.,9999.]
ncalb.variables['topg'].source=bedmap2_link
ncalb.variables['topg'].reference="Fretwell et al. (2013), Bedmap2: improved ice bed, surface and thickness datasets for Antarctica. see publication http://www.the-cryosphere.net/7/375/2013/tc-7-375-2013.pdf"

ncvel = ncalb.createVariable( 'velocity','float32',('y','x') )
ncvel.units=ncbm2.variables['velocity'].units
ncvel.long_name=ncbm2.variables['velocity'].long_name
ncvel.standard_name=ncbm2.variables['velocity'].standard_name
ncvel.source = rignot_link
ncvel.reference = "Rignot, E., J. Mouginot, and B. Scheuchl (2011), Ice Flow of the Antarctic Ice Sheet, Science, doi 10.1126/science.1208336."
#velbm2=np.ma.masked_array(velbm2, mask=(velbm2==-9999.))
ncalb.variables['velocity'][:] = velbm2

ncalb.variables['precipitation'][:] = accumbm2
ncalb.variables['precipitation'].long_name = "accumulation after Arthern et al."
ncalb.variables['precipitation'].units = "kg/m2/year"
ncalb.variables['precipitation'].reference = "Arthern, R. J., D. P. Winebrenner, and D. G. Vaughan (2006), Antarctic snow accumulation mapped using polarization of 4.3-cm wavelength microwave emission, J. Geophys. Res., 111, D06107, doi:10.1029/2004JD005667."
ncalb.variables['precipitation'].source = arthern_link

now = datetime.datetime.now().strftime("%B %d, %Y")
ncalb.mergeComment = authors+" merged Albmap and Bedmap2 and Rignot velocities and Arthern accumulation at " + now

ncalb.close()
ncbm2.close()



