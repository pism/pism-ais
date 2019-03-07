"""
torsten.albrecht@pik-potsdam.de
Download basal heatflux data from Purucker, 2013 for given lon-lat points and triangulate to PISM compatible grid.

Purucker, M. E.: Geothermal heat flux data set based on low res- olution observations collected by the CHAMP satellite be- tween
2000 and 2010, and produced from the MF-6 model following the technique described in Fox Maule et al. (2005), 2013

retrieved from:
https://core2.gsfc.nasa.gov/research/purucker/heatflux_updates.html
"""


import os, sys
import numpy as np
import numpy.ma as ma
import netCDF4 as nc
import subprocess

## this hack is needed to import config.py from the project root
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if project_root not in sys.path: sys.path.append(project_root)

import config as cf; reload(cf)
import pism_input.pism_input as pi; reload(pi)

plot_figure=False

dataset = "bheatflx_purucker"
data_path = os.path.join(cf.output_data_path, dataset)
final_filename = os.path.join(data_path,"bheatflx_purucker_"+cf.grid_id+".nc")

data_file = "heatflux_mf7_foxmaule05.txt"
link = "https://core2.gsfc.nasa.gov/research/purucker/"+data_file


## such file definitions should go to config.py, so that other functions can access them.

# if data is not yet there, download
downloaded_file = os.path.join(data_path,data_file)
if not os.path.isfile(downloaded_file):
  print "Downloading basal heatflux data."
  os.system("mkdir " + data_path)
  os.system("wget -N " + link + " -P " + data_path)

gridfile = os.path.join(cf.cdo_remapgridpath,cf.grid_id+'.nc')
griddat = nc.Dataset(gridfile, 'r')
x_pism = griddat.variables["x"][:]
y_pism = griddat.variables["y"][:]
#m_pism = griddat.variables["mask"][:]
#print len(x_pism),len(y_pism),np.shape(m_pism)
griddat.close()
xi, yi = np.meshgrid(x_pism,y_pism)


# Read data from txt file:
lon = np.zeros(0)
lat = np.zeros(0)
hfl = np.zeros(0)

f = open(os.path.join(data_path, data_file), 'r')

for i,line in enumerate(f):
        line = line.strip()
        columns = line.split()
        lon=np.append(lon, float(columns[0]))
        lat=np.append(lat, float(columns[1]))
        hfl=np.append(hfl, float(columns[2]))

f.close()

# Define Projection ######################################################################

import pyproj

# Define a projection with Proj4 notation from Bedmachine
#http://spatialreference.org/ref/epsg/wgs-84-antarctic-polar-stereographic/
wgs84=pyproj.Proj(cf.proj4str)
x0,y0 = wgs84(180.0,-90.0)
x, y = wgs84(lon[:],lat[:])


# Triangulation #####################################

#https://matplotlib.org/gallery/images_contours_and_fields/triinterp_demo.html#sphx-glr-gallery-images-contours-and-fields-triinterp-demo-py
import matplotlib.tri as mtri
triang = mtri.Triangulation(x, y)
interp_lin = mtri.LinearTriInterpolator(triang, hfl)
heatflux = interp_lin(xi, yi)

#create array for these dimensions and fill in values:
fillvalue = 70.0 #np.nan
heatflux = np.ma.masked_array(heatflux,fill_value=fillvalue)
heatflux[heatflux.mask==True] = fillvalue

# Plot ############################################

if plot_figure:

  import matplotlib.pyplot as plt

  fig = plt.figure(figsize=(8,12))

  plt.subplot(211)
  plt.tricontourf(triang, hfl,np.arange(0,185,5),extent="both")
  plt.triplot(triang, 'ko-',lw=0.1,alpha=0.3,markersize=2)
  plt.title('Triangular grid')
  plt.colorbar()
  plt.axis([x_pism[0],x_pism[-1],y_pism[0],y_pism[-1]])
  plt.xticks([])
  plt.yticks([])

  # Plot linear interpolation to quad grid.
  plt.subplot(212)
  plt.contourf(xi, yi, heatflux,np.arange(0,185,5),extent="both")
  #plt.plot(xi, yi, 'k-', lw=0.1, alpha=0.2)
  #plt.plot(xi.T, yi.T, 'k-', lw=0.1, alpha=0.2)
  plt.title("Linear interpolation")
  plt.colorbar()
  #plt.contour(xi, yi, m_pism[:,:],[0.5,1.5],colors='k', lws=0.5, alpha=0.9,zorder=10)
  plt.axis([x_pism[0],x_pism[-1],y_pism[0],y_pism[-1]])
  plt.xticks([])
  plt.yticks([])

  plt.tight_layout()
  #plt.savefig("interpolation.pdf",format="pdf")
  plt.show()


# Save data as NetCDF file ####################################

wrtfile = nc.Dataset(final_filename, 'w', format='NETCDF4_CLASSIC')
wrtfile.createDimension('x', size=len(x_pism))
wrtfile.createDimension('y', size=len(y_pism))
wrtfile.createDimension('time', size=None)
wrtfile.createDimension('nv', size=2)

nct     = wrtfile.createVariable('time', 'f8', ('time',))
nctb    = wrtfile.createVariable('time_bnds', 'f8', ('time','nv'))
ncx   = wrtfile.createVariable('x', 'f8', ('x',))
ncy   = wrtfile.createVariable('y', 'f8', ('y',))
nchfl   = wrtfile.createVariable('bheatflx', 'f8', ('time','y', 'x'))

tm = np.arange(1,1.5,1)
nct[:] = tm
nct.units = 'years since 01-01-01 00:00:00'
nct.bounds = 'time_bnds'
nct.calendar = "365_day"

ncy.units = 'meters'
ncx.units = 'meters'

t_bds = np.zeros([len(tm),2])
t_bds[:,0] = tm -5
t_bds[0,:] = tm + 5
nctb[:] = t_bds
nctb.units = nct.units
nctb.calendar = nct.calendar

ncy[:] = y_pism
ncy.long_name = 'Y-coordinate in Cartesian system'
ncy.standard_name = "projection_y_coordinate"
ncy.spacing_meters = 1000.
ncy.units = 'm'
ncx[:] = x_pism
ncx.long_name = 'X-coordinate in Cartesian system'
ncx.units = 'm'
ncx.standard_name = "projection_x_coordinate"
ncx.spacing_meters = 1000.


nchfl[0,:] = heatflux[:]/1000.0
nchfl.units = 'W m-2'
nchfl.long_name = "geothermal heat flux - Purucker, 2013" ;

import datetime
now = datetime.datetime.now().strftime("%B %d, %Y")
wrtfile.inputdata = link
wrtfile.proj4 = cf.proj4str
wrtfile.projection = cf.proj4str
wrtfile.comment  = cf.authors+" created netcdf geothermal heatflux file at " + now
wrtfile.institution = 'Potsdam Institute for Climate Impact Research (PIK), Germany'
wrtfile.citation = "Purucker, M. E.: Geothermal heat flux data set based on low res- olution observations collected by the CHAMP satellite between 2000 and 2010, and produced from the MF-6 model following the technique described in Fox Maule et al. (2005), 2013"

wrtfile.close()

#subprocess.check_call('python ../tools/nc2cdo.py '+final_filename,shell=True)
subprocess.check_call('ncks -A -v lon,lat '+gridfile+' '+final_filename,shell=True)


print 'Data successfully saved to', final_filename
