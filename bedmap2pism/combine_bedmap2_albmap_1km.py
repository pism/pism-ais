#!/usr/bin/env python

# This file is part of bedmap2pism.
# Copyright (C) 2015-2016 
authors="matthias.mengel@pik-potsdam.de and torsten.albrecht@pik-potsdam.de"

# creates combined 1km dataset NetCDF file from Bedmap2 and Albmap data
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

ncalb_bm2_name = "albmap_bedmap2_1km.nc"
ncalb_name="albmap_data/albmap_pism_5km.nc"
ncbm2_name="bedmap2_data/bedmap2_1km.nc"
bedmap2_link="https://secure.antarctica.ac.uk/data/bedmap2/bedmap2_bin.zip"
rignot_link="https://secure.antarctica.ac.uk/data/bedmap2/resources/Rignot_velocity/bin.zip"
arthern_link="https://secure.antarctica.ac.uk/data/bedmap2/resources/Arthern_accumulation/Arthern_accumulation_bin.zip"
x0_bm2= 3333500.

compare_bm2_alb = True
convert_to_64bit = False
chunklen = 550

#############################################################################################################
### Combine Albmap and Bedmap2
print "Making PISM-readable file combining Albmap and Bedmap2 on 1km resolution"

os.system("cp " + ncbm2_name +" "+ ncalb_bm2_name)

ncalb = nc.Dataset(ncalb_name, 'r')
ncbm2 = nc.Dataset(ncalb_bm2_name, 'a')

# get mask, topo, thk from bedmap 2
xbm2    = ncbm2.variables['x'][:]
ybm2    = ncbm2.variables['y'][:]
topg    = ncbm2.variables['topg'][:]
thk     = ncbm2.variables['thk'][:]
mask    = ncbm2.variables['mask'][:]
usurf   = ncbm2.variables['usurf'][:]
bedunc  = ncbm2.variables['bedunc'][:]

# get surface rignot velocity from bedmap 2
vel = ncbm2.variables['velocity'][:]

# get arthern accumulation from bedmap 2
#accumulation = ncbm2.variables['accum'][:]

# adjust bedm2 to centered x,y, see bedmap2 readme file
xbm2 -= x0_bm2
ybm2 -= x0_bm2 #- 395000.
xgrid, ygrid = np.meshgrid(xbm2,ybm2)


# get others from albmap
xalb      = ncalb.variables['x'][:]
yalb      = ncalb.variables['y'][:]
precip    = ncalb.variables['precipitation'][:][0]
artm      = ncalb.variables['air_temp'][:][0]
bheatflx  = ncalb.variables['bheatflx'][:][0]
lon       = ncalb.variables['lon'][0,:,:]
lat       = ncalb.variables['lat'][0,:,:]
#xgrid, ygrid = np.meshgrid(xalb,yalb)

ncairtemp       = ncbm2.createVariable( 'air_temp','float32',('y','x') )
ncprecipitation = ncbm2.createVariable( 'precipitation','float32',('y','x') )
ncbheatfl       = ncbm2.createVariable( 'bheatflx','float32',('y','x') )
ncmap           = ncbm2.createVariable( 'mapping','c' )
nclon           = ncbm2.createVariable( 'lon','float32',('y','x') )
nclat           = ncbm2.createVariable( 'lat','float32',('y','x') )


if compare_bm2_alb:

  thk_alb        = ncalb.variables['thk'][0,:,:]
  topg_alb       = ncalb.variables['topg'][0,:,:]

  ncthkold    = ncbm2.createVariable( 'thk_alb','float32',('y','x') )
  nctopgold   = ncbm2.createVariable( 'topg_alb','float32',('y','x') )
  ncthkdiff   = ncbm2.createVariable( 'thk_diff','float32',('y','x') )
  nctopgdiff  = ncbm2.createVariable( 'topg_diff','float32',('y','x') )

  ncthkold.units    = ncalb.variables['thk'].units
  nctopgold.units   = ncalb.variables['thk'].units
  ncthkdiff.units   = ncalb.variables['thk'].units
  nctopgdiff.units  = ncalb.variables['thk'].units

  ncthkdiff.standard_name = "bedm2_alb_thk"
  ncthkdiff.long_name = "bedmap2 albmap thickness difference"
  nctopgdiff.standard_name = "bedm2_alb_topg"
  nctopgdiff.long_name = "bedmap2 albmap topograpy difference"



ncmap.proj4      = "+proj=stere +lon_0=0 +lat_0=-90 +lat_ts=-71 +ellps=WGS84 +datum=WGS84"
ncbm2.projection = "+proj=stere +lon_0=0 +lat_0=-90 +lat_ts=-71 +ellps=WGS84 +datum=WGS84"


print "Write Albmap attributes"
for varname in ['thk','topg','air_temp','precipitation','bheatflx','lon','lat','mapping']:
  #print varname
  for att in ncalb.variables[varname].ncattrs():
    exec "ncbm2.variables[varname]." + att + "= ncalb.variables[varname]." + att


#nc.reDef(ncbm2)
# save memory by writing in chunks.
chunklen = 550
for xi in range(0,len(xbm2),chunklen):
  for yi in range(0,len(ybm2),chunklen):

    if yi==0:
      print " %.0f %%" % (100.0*xi/len(xbm2))
    #print xi, yi

    xgrid_red = xgrid[xi:xi+chunklen,yi:yi+chunklen]
    ygrid_red = ygrid[xi:xi+chunklen,yi:yi+chunklen]

    # regrid albmap stuff
    ncairtemp[xi:xi+chunklen,yi:yi+chunklen] = interp(artm, xalb, yalb, xgrid_red, ygrid_red)
    ncprecipitation[xi:xi+chunklen,yi:yi+chunklen] = interp(precip, xalb, yalb, xgrid_red, ygrid_red)
    ncbheatfl[xi:xi+chunklen,yi:yi+chunklen] = interp(bheatflx, xalb, yalb, xgrid_red, ygrid_red)
    nclon[xi:xi+chunklen,yi:yi+chunklen] = interp(lon, xalb, yalb, xgrid_red, ygrid_red)
    nclat[xi:xi+chunklen,yi:yi+chunklen] = interp(lat, xalb, yalb, xgrid_red, ygrid_red)

    # write flipped bedm2 data
    #ncbm2.variables['topg'][xi:xi+chunklen,yi:yi+chunklen]  = topg[xi:xi+chunklen,yi:yi+chunklen]
    #ncbm2.variables['thk'][xi:xi+chunklen,yi:yi+chunklen]   = thk[xi:xi+chunklen,yi:yi+chunklen]
    #ncbm2.variables['mask'][xi:xi+chunklen,yi:yi+chunklen]  = mask[xi:xi+chunklen,yi:yi+chunklen]

    #ncbm2.variables['usurf'][xi:xi+chunklen,yi:yi+chunklen]  = usurf[xi:xi+chunklen,yi:yi+chunklen]
    #ncbm2.variables['bedunc'][xi:xi+chunklen,yi:yi+chunklen]  = bedunc[xi:xi+chunklen,yi:yi+chunklen]
    #ncbm2.variables['velocity'][xi:xi+chunklen,yi:yi+chunklen]  = vel[xi:xi+chunklen,yi:yi+chunklen]

    if compare_bm2_alb:

      ncthkold[xi:xi+chunklen,yi:yi+chunklen] = interp(thk_alb, xalb, yalb, xgrid_red, ygrid_red)
      nctopgold[xi:xi+chunklen,yi:yi+chunklen] = interp(topg_alb, xalb, yalb, xgrid_red, ygrid_red)


ncbm2.variables['x'][:] = xbm2
ncbm2.variables['y'][:] = ybm2


### clean up (to avoid errors in the PISM boot)
thk[thk > 10000.] = 0.
topg[topg > 9000.] = -9990.
topg[topg < -9000.] = -9999.

nclon[nclon < -180.] = -180.0
nclon[nclon > 180.] = 180.0
nclat[nclat < -90.] = -90.0
nclat[nclat > 90.] = 90.0

ncbheatfl[ncbheatfl < 0.] = 0.
ncbheatfl[ncbheatfl > 0.2] = 0.2
ncprecipitation[ncprecipitation < 0.] = 0.0
ncprecipitation[ncprecipitation > 10.] = 10.0
ncairtemp[ncairtemp < 200.] = 200.0
ncairtemp[ncairtemp > 280.] = 280.0



### add some difference field to compare albmap and bedma2 topg and thk
if compare_bm2_alb:

  #ncthkold[thk == -9999.] = -9999.
  #nctopgold[topg == -9999.] = -9999.
  ncthkdiff[:]   = thk[:] - ncthkold[:]
  nctopgdiff[:]  = topg[:] - nctopgold[:]


  #mask?

#ncbm2.variables['thk'].valid_range = [0.,9999.]
ncbm2.variables['thk'].source=bedmap2_link
ncbm2.variables['thk'].reference="Fretwell et al. (2013), Bedmap2: improved ice bed, surface and thickness datasets for Antarctica. see publication http://www.the-cryosphere.net/7/375/2013/tc-7-375-2013.pdf"

ncbm2.variables['topg'].valid_range = [-9999.,9999.]
ncbm2.variables['topg'].source=bedmap2_link
ncbm2.variables['topg'].reference="Fretwell et al. (2013), Bedmap2: improved ice bed, surface and thickness datasets for Antarctica. see publication http://www.the-cryosphere.net/7/375/2013/tc-7-375-2013.pdf"

ncbm2.variables['velocity'].source = rignot_link
ncbm2.variables['velocity'].reference = "Rignot, E., J. Mouginot, and B. Scheuchl (2011), Ice Flow of the Antarctic Ice Sheet, Science, doi 10.1126/science.1208336."
#velbm2=np.ma.masked_array(velbm2, mask=(velbm2==-9999.))

ncbm2.variables['precipitation'].reference = "Arthern, R. J., D. P. Winebrenner, and D. G. Vaughan (2006), Antarctic snow accumulation mapped using polarization of 4.3-cm wavelength microwave emission, J. Geophys. Res., 111, D06107, doi:10.1029/2004JD005667."
ncbm2.variables['precipitation'].source = arthern_link

now = datetime.datetime.now().strftime("%B %d, %Y")
ncbm2.mergeComment = authors+" merged Albmap and Bedmap2 and Rignot velocities and Arthern accumulation at " + now


ncbm2.close()
ncalb.close()

if convert_to_64bit:
  os.system("ncpdq -O --fl_fmt=64bit "+ncalb_bm2_name +" "+ ncalb_bm2_name)

