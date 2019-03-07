"""
matthias.mengel@pik, torsten.albrecht@pik
Download basal heatflux data from Martos et al., 2017 and rename variables to make the PISM compatible.
Data are documented here:
https://doi.pangaea.de/10.1594/PANGAEA.882503

Paper is found here:
http://onlinelibrary.wiley.com/doi/10.1002/2017GL075609/abstract

DATA DESCRIPTION:
Citation: Martos, Yasmina M (2017): Antarctic geothermal heat flux distribution and estimated Curie Depths, links to gridded files. PANGAEA, https://doi.pangaea.de/10.1594/PANGAEA.882503 (DOI registration in progress), 
  Supplement to: Martos, Yasmina M; Catalan, Manuel; Jordan, Tom A; Golynsky, Alexander V; Golynsky, Dmitry A; Eagles, Graeme; Vaughan, David G (accepted): Heat flux distribution of Antarctica unveiled. Geophysical Research Letters, https://doi.org/10.1002/2017GL075609
Abstract: Antarctica is the largest reservoir of ice on Earth. Understanding its ice sheet dynamics is crucial to unraveling past global climate change and making robust climatic and sea level predictions. Of the basic parameters that shape and control ice flow, the most poorly known is geothermal heat flux. Direct observations of heat flux are difficult to obtain in Antarctica, and until now continent-wide heat flux maps have only been derived from low-resolution satellite magnetic and seismological data. We present a high resolution heat flux map and associated uncertainty derived from spectral analysis of the most advanced continental compilation of airborne magnetic data. Small-scale spatial variability and features consistent with known geology are better reproduced than in previous models, between 36% and 50%. Our high-resolution heat-flux map and its uncertainty distribution provide an important new boundary condition to be used in studies on future subglacial hydrology, ice-sheet dynamics and sea-level change.
Coverage: LATITUDE: -90.000000 * LONGITUDE: 0.000000
Event(s): pan-Antarctica * LATITUDE: -90.000000 * LONGITUDE: 0.000000
Comment:  Gridded data at 15 km. X and Y coordinates are in polar stereographic projection.
Parameter(s): File content (Content) * PI: Martos, Yasmina M (yasmartos@gmail.com)
  File name (File name) * PI: Martos, Yasmina M (yasmartos@gmail.com)
  File format (File format) * PI: Martos, Yasmina M (yasmartos@gmail.com)
  File size [kByte] (File size) * PI: Martos, Yasmina M (yasmartos@gmail.com)
  Uniform resource locator/link to file (URL file) * PI: Martos, Yasmina M (yasmartos@gmail.com)
License:  Creative Commons Attribution 3.0 Unported (CC-BY)
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

dataset = "bheatflx_martos"
data_path = os.path.join(cf.output_data_path, dataset)
final_filename = os.path.join(data_path,"bheatflx_martos_15km_input.nc")

data_file = "Antarctic_GHF.xyz"
data_u_file = data_file.replace(".xyz","_uncertainty.xyz")

link = "http://store.pangaea.de/Publications/Martos-etal_2017/"+data_file
linku = "http://store.pangaea.de/Publications/Martos-etal_2017/Antarctic_GHF_uncertainty.xyz"


## such file definitions should go to config.py, so that other functions can access them.

# if data is not yet there, download
downloaded_file = os.path.join(data_path,data_file)
if not os.path.isfile(downloaded_file):
  print "Downloading albmap data."
  os.system("mkdir " + data_path)
  os.system("wget -N " + link + " -P " + data_path)
  os.system("wget -N " + linku + " -P " + data_path)


# Read data from txt file:
xval = np.zeros(0)
yval = np.zeros(0)
bhflx = np.zeros(0)
bhflxu = np.zeros(0)


f = open(os.path.join(data_path, data_file), 'r')
for i,line in enumerate(f):
        line = line.strip()
        columns = line.split()
        xval=np.append(xval, float(columns[0]))
        yval=np.append(yval, float(columns[1]))
        bhflx=np.append(bhflx, float(columns[2]))
f.close()

f = open(os.path.join(data_path, data_u_file), 'r')
for i,line in enumerate(f):
        line = line.strip()
        columns = line.split()
        #xvalu=np.append(xvalu, float(columns[0]))
        #yvalu=np.append(yvalu, float(columns[1]))
        bhflxu=np.append(bhflxu, float(columns[2]))
f.close()

# Convert to pism grid...
x_new = np.sort(np.unique(xval)) #are sorted already
y_new = np.sort(np.unique(yval))

print "Create a "+ str(len(x_new)) +" x " + str(len(y_new)) +" grid on 15km resolution"

#create array for these dimensions and fill in values:
#fillvalue = 70.0
#fillvalue = np.nan
fillvalue = -9.e+33
bheatflx = np.zeros((len(y_new), len(x_new))) + fillvalue
bheatflxu = np.zeros((len(y_new), len(x_new))) + fillvalue

for i in range(len(bhflx)):
    # go through all vals and fill them into the right place
    ix = np.in1d(x_new.ravel(), xval[i]).reshape(x_new.shape)
    iy = np.in1d(y_new.ravel(), yval[i]).reshape(y_new.shape)
    bheatflx[iy,ix] = bhflx[i]*1.0e-3
    bheatflxu[iy,ix] = bhflxu[i]*1.0e-3

# Save data as NetCDF file
wrtfile = nc.Dataset(final_filename, 'w', format='NETCDF4_CLASSIC')
wrtfile.createDimension('x', size=len(x_new))
wrtfile.createDimension('y', size=len(y_new))

ncx   = wrtfile.createVariable('x', 'f4', ('x',))
ncy   = wrtfile.createVariable('y', 'f4', ('y',))
ncbhflx  = wrtfile.createVariable('bheatflx', 'f4', ('y', 'x'),fill_value=fillvalue)
ncbhflxu  = wrtfile.createVariable('bheatflx_uncertainty', 'f4', ('y', 'x'),fill_value=fillvalue)

ncbhflx[:] = bheatflx
ncbhflx.long_name = "geothermal heat flux - Martos et al., 2017" ;
ncbhflx.units = "W m-2" 
ncbhflx.source = data_file+" downloaded from https://doi.pangaea.de/10.1594/PANGAEA.882503" ;
ncbhflx.reference = "Martos, Yasmina M; Catalan, Manuel; Jordan, Tom A; Golynsky, Alexander V; Golynsky, Dmitry A; Eagles, Graeme; Vaughan, David G (accepted): Heat flux distribution of Antarctica unveiled. Geophysical Research Letters, https://doi.org/10.1002/2017GL075609"
ncbhflx.missing_value = fillvalue

ncbhflxu[:] = bheatflxu
ncbhflxu.units = ncbhflx.units
ncbhflx.long_name = "uncertainty of " + ncbhflx.long_name

ncx[:] = x_new
ncy[:] = y_new
ncx.units = "meters"
ncy.units = "meters"

wrtfile.proj4 = cf.proj4str

wrtfile.close()

# prepare the input file for cdo remapping
# this step takes a while for high resolution data (i.e. 1km)
pi.prepare_ncfile_for_cdo(final_filename)
#subprocess.check_call('python ../tools/nc2cdo.py '+final_filename,shell=True)

print 'Data successfully saved to', final_filename
