"""
matthias.mengel@pik, torsten.albrecht@pik
Download Bedmap 2 data and save to (1km) netcdf file.
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
import pism_input.pism_input as pi; reload(pi)

### Bedmap2   ##########################################################
# Documentation of the data: https://www.bas.ac.uk/project/bedmap-2/
bedmap2_link="https://secure.antarctica.ac.uk/data/bedmap2/bedmap2_bin.zip"
reference_short = "Fretwell, P. et al: Bedmap2: improved ice bed, surface and thickness datasets for Antarctica, The Cryosphere, 7, 375-393, doi:10.5194/tc-7-375-2013, 2013."
## such file definitions should go to config.py, so that other functions can access them.
bedmap2_data_path = os.path.join(cf.output_data_path, "bedmap2")
ncout_name = os.path.join(bedmap2_data_path, 'bedmap2_1km_input.nc')

# if data is not yet extracted in bedmap2_bin
if not os.path.exists(os.path.join(bedmap2_data_path,"bedmap2_bin")):
  print "Downloading bedmap2 binary data."
  os.system("mkdir " + bedmap2_data_path)
  os.system("wget -N --no-check-certificate " + bedmap2_link + " -P " + bedmap2_data_path)
  os.system("cd "+bedmap2_data_path+" && unzip bedmap2_bin.zip")

data_files = {"topg":"bedmap2_bed.flt",
              "thk":"bedmap2_thickness.flt",
              "mask":"bedmap2_icemask_grounded_and_shelves.flt",
              "bedunc":"bedmap2_grounded_bed_uncertainty.flt",
              "usurf":"bedmap2_surface.flt"}

data_fills = {"thk":0.0,
              "topg":-5000.0,
              "usurf":0.0,
              "bedunc":0.0,
              "mask":2.0}

# taken from bedmap2 readme
N=6667
dx = 1000.0 #m
dy = 1000.0 #m
x = np.linspace(-(N-1)*dx/2.0,(N-1)*dx/2.0,N)
y = np.linspace(-(N-1)*dy/2.0,(N-1)*dy/2.0,N)


print "Reading bedmap2 binary files from %s ...\n" % (bedmap2_data_path)

bedm2_vars = {}
for var, file in data_files.iteritems():
  fname = os.path.join(bedmap2_data_path,"bedmap2_bin",file)
  vardata = np.flipud(np.ma.masked_equal(np.reshape(
          np.fromfile(fname,dtype=np.float32),(N,N)),-9999.0))

  print " range of "+str(var)+" = [%.2f, %.2f]" % (vardata.min(),vardata.max())
  #get rid off NaN
  vardata[vardata.mask]=data_fills[var]

  bedm2_vars[var] =  vardata


bedm2_attributes = {"topg": {"long_name" : "elevation of bedrock",
                           "standard_name" : "bedrock_altitude",
                           "units" : "meters", "reference" : reference_short },
                    "usurf": {"long_name" : "ice upper surface elevation",
                            "standard_name" : "surface_altitude",
                            "units" : "meters", "reference" : reference_short },
                    "thk": {"long_name" : "thickness of ice sheet or ice shelf",
                            "standard_name" : "land_ice_thickness",
                            "units" : "meters", "reference" : reference_short },
                    "bedunc": {"long_name" : "uncertainty of bed topography",
                            "standard_name" : "bed_uncertainty",
                            "units" : "meters", "reference" : reference_short},
                    "mask": {"long_name" : "ice-type (ice-free/grounded/floating/ocean) integer mask",
                            "standard_name" : "mask",
                            "units" : "", "reference" : reference_short} }

print "Writing NetCDF file '%s' ..." % ncout_name
ncout = netCDF4.Dataset(ncout_name,"w",format='NETCDF4_CLASSIC')

# no time dimension needed here
ncout.createDimension('x',size=len(x))
ncout.createDimension('y',size=len(x))
ncx   = ncout.createVariable( 'x','float64',('x',) )
ncy   = ncout.createVariable( 'y','float64',('y',) )
ncx[:] = x
ncy[:] = y

for varname,data in bedm2_vars.iteritems():

  ncvar = ncout.createVariable( varname,'float64',('y','x') ) #double precision
  ncvar[:] = data
  for att in bedm2_attributes[varname]:
    setattr(ncvar,att,bedm2_attributes[varname][att])

now = datetime.datetime.now().strftime("%B %d, %Y")
#antarctica
ncout.proj4 = cf.proj4str
ncout.comment  = cf.authors+" created netcdf bedmap2 file at " + now
ncout.reference = "Fretwell, P., Pritchard, H. D., Vaughan, D. G., Bamber, J. L., Barrand, N. E., Bell, R., Bianchi, C., Bingham, R. G., Blankenship, D. D., Casassa, G., Catania, G., Callens, D., Conway, H., Cook, A. J., Corr, H. F. J., Damaske, D., Damm, V., Ferraccioli, F., Forsberg, R., Fujita, S., Gim, Y., Gogineni, P., Griggs, J. A., Hindmarsh, R. C. A., Holmlund, P., Holt, J. W., Jacobel, R. W., Jenkins, A., Jokat, W., Jordan, T., King, E. C., Kohler, J., Krabill, W., Riger-Kusk, M., Langley, K. A., Leitchenkov, G., Leuschen, C., Luyendyk, B. P., Matsuoka, K., Mouginot, J., Nitsche, F. O., Nogi, Y., Nost, O. A., Popov, S. V., Rignot, E., Rippin, D. M., Rivera, A., Roberts, J., Ross, N., Siegert, M. J., Smith, A. M., Steinhage, D., Studinger, M., Sun, B., Tinto, B. K., Welch, B. C., Wilson, D., Young, D. A., Xiangbin, C., and Zirizzotti, A.: Bedmap2: improved ice bed, surface and thickness datasets for Antarctica, The Cryosphere, 7, 375-393, doi:10.5194/tc-7-375-2013, 2013."
ncout.source = bedmap2_link
ncout.documentation = "https://www.bas.ac.uk/project/bedmap-2/"
ncout.close()

#remove inconsistency of ice thickness at Lake Wostok
os.system("ncap2 -O -s 'where(mask==0) thk=usurf-topg' "+ncout_name+" "+ncout_name)

# prepare the input file for cdo remapping
# this step takes a while for high resolution data (i.e. 1km)
pi.prepare_ncfile_for_cdo(ncout_name)

