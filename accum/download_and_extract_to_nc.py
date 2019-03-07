"""
matthias.mengel@pik, torsten.albrecht@pik
Download Arthern accumulation data and save to (1km) netcdf file.
"""

import os, sys
import numpy as np
import sys
import netCDF4
import datetime

## this hack is needed to import config.py from the project root
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if project_root not in sys.path: sys.path.append(project_root)
import config as cf; reload(cf)

### Arthern accumulation   ##########################################################
# Documentation of the data: https://www.bas.ac.uk/project/bedmap-2/
## such file definitions should go to config.py, so that other functions can access them.
arthern_link="https://secure.antarctica.ac.uk/data/bedmap2/resources/Arthern_accumulation/Arthern_accumulation_bin.zip"
arthern_data_path = os.path.join(cf.output_data_path, "accum")
ncout_name = os.path.join(arthern_data_path, 'accum_1km_input.nc')

# if data is not yet extracted in bedmap2_bin
if not os.path.exists(os.path.join(arthern_data_path,"accum_bin")):
  print "Downloading arthern acumulation binary data."
  os.system("mkdir -p " + os.path.join(arthern_data_path, 'accum_bin'))
  os.system("wget -N " + arthern_link + " -P " + arthern_data_path+'/')
  os.system("unzip "+os.path.join(arthern_data_path,"Arthern_accumulation_bin.zip")+" -d "+os.path.join(arthern_data_path, 'accum_bin/'))
  os.system("rm "+os.path.join(arthern_data_path,"Arthern_accumulation_bin.zip"))

if os.path.isfile(ncout_name):
  print "Accum file", ncout_name
  print "was already written, do nothing."
  sys.exit(0)

data_files = {"accum":"arthern_accumulation_bedmap2_grid.flt",
              "rms":"arthern_accumulation_rms_bedmap2_grid.flt"}

data_fills = {"accum":0.0,
              "rms":0.0}

# taken from arthern readme
Nx=7899
Ny=8300
x0_art=-3949500.
y0_art=x0_art
dx = 1000.0 #m
dy = dx
#x = np.linspace(-(N-1)*dx/2.0,(N-1)*dx/2.0,N)
#y = np.linspace(-(N-1)*dy/2.0,(N-1)*dy/2.0,N)
x = np.arange(x0_art,x0_art+(Nx)*dx,dx)
y = np.arange(y0_art,y0_art+(Ny)*dy,dy)


print "Reading arthern binary files from %s ...\n" % (arthern_data_path)

accum_vars = {}
for var, file in data_files.iteritems():
  fname = os.path.join(arthern_data_path,"accum_bin",file)
  vardata = np.flipud(np.ma.masked_equal(np.reshape(
          np.fromfile(fname,dtype=np.float32),(Ny,Nx)),-9999.0))
  print np.shape(vardata)

  print " range of "+str(var)+" = [%.2f, %.2f]" % (vardata.min(),vardata.max())
  #get rid off NaN
  vardata[vardata.mask]=data_fills[var]

  accum_vars[var] =  vardata


accum_attributes = {"accum": {"long_name" : "snow accumulation",
                           "units" : "m year-1"},
                    "rms": {"long_name" : "RMSE snow accumulation",
                            "units" : "m year-1"}}

print "Writing NetCDF file '%s' ..." % ncout_name
ncout = netCDF4.Dataset(ncout_name,"w",format='NETCDF4_CLASSIC')

# no time dimension needed here
ncout.createDimension('x',size=len(x))
ncout.createDimension('y',size=len(y))
ncx   = ncout.createVariable( 'x','float64',('x',) )
ncy   = ncout.createVariable( 'y','float64',('y',) )
ncx[:] = x
ncy[:] = y

for varname,data in accum_vars.iteritems():

  ncvar = ncout.createVariable( varname,'float64',('y','x') ) #double precision
  ncvar[:] = data
  for att in accum_attributes[var]:
    pass

now = datetime.datetime.now().strftime("%B %d, %Y")
#antarctica
ncout.proj4 = cf.proj4str
ncout.comment  = cf.authors+" created netcdf accumulation file at " + now
ncout.citation = "Arthern, R. J., D. P. Winebrenner, and D. G. Vaughan (2006), Antarctic snow accumulation mapped using polarization of 4.3-cm wavelength microwave emission, J. Geophys. Res., 111, D06107, doi:10.1029/2004JD005667."

ncout.close()
print "Done"

