#!/usr/bin/env python

"""
matthias.mengel@pik, torsten.albrecht@pik
Download Mouginot velocity data and save to 1km netcdf file.

Preprocess the ~12x1.3GB Antarctic ice velocity dataset from NASA MEASURES project
https://nsidc.org/data/NSIDC-0720/

See infos http://www.mdpi.com/2072-4292/9/4/364
Currently not available without Earthdata login
"""


import numpy as np
import netCDF4 as NC
#import datetime
import os, sys, glob
import subprocess

## this hack is needed to import config.py from the project root
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if project_root not in sys.path: sys.path.append(project_root)
import config as cf; reload(cf)
import pism_input.pism_input as pi; reload(pi)

dataset = "annualvel"
data_path = os.path.join(cf.output_data_path, dataset)

if not os.path.exists(data_path): os.makedirs(data_path)

output_filename = os.path.join(data_path,dataset+"_1km_input.nc")

input_filenames = "/p/projects/tumble/pism_input/Velocity/mouginot_rignot17_annual/Antarctica_ice_velocity_*_1km_v1.nc"
url = "https://daacdata.apps.nsidc.org/pub/DATASETS/nsidc0720_MEASURES_yearly_antarc_vel_v01/"


def run(commands):
    """Run a list of commands (or one command given as a string)."""
    if isinstance(commands, (list, tuple)):
        for cmd in commands:
            print "Running '%s'..." % cmd
            subprocess.call(cmd.split(' '))
    else:
        run([commands])


def preprocess_ice_velocity(fn,fnout,cnt):


    commands = [#"wget -nc %s%s" % (url, input_filename), # NSIDC supports compression on demand!
                "cp %s %s" % (fn,fnout),
                "ncks -O -3 -v VX,VY,lon,lat,coord_system %s %s" % (fnout, fnout),
                "ncrename -v VX,vx -v VY,vy " + fnout
                ]


    if not os.path.exists(fnout):
        run(commands)

    nc = NC.Dataset(fnout, 'a')

    if 'time' not in nc.variables:
        nc.createDimension('time', size=None)
        nct = nc.createVariable('time', 'f8', ('time',))
        nct[:] = cnt
        nct.units = 'years since 01-01-01 00:00:00'


    # Metadata provided with the dataset describes the *full* grid, so it is a
    # lot easier to modify this file instead of adding grid information to the
    # "cutout" file.
    if 'x' not in nc.variables and 'y' not in nc.variables:
        nx = nc.nx
        ny = nc.ny
        x_min = float(nc.xmin.strip().split(' ')[0])
        y_max = float(nc.ymax.strip().split(' ')[0])
        x_max = y_max
        y_min = x_min

        x = np.linspace(x_min, x_max, nx)
        y = np.linspace(y_max, y_min, ny)

        nc.projection = "+proj=stere +ellps=WGS84 +datum=WGS84 +lon_0=0 +lat_0=-90 +lat_ts=-71 +units=m"

        try:
            x_var = nc.createVariable('x', 'f8', ('x',))
            y_var = nc.createVariable('y', 'f8', ('y',))
        except:
            x_var = nc.variables['x']
            y_var = nc.variables['y']

        x_var[:] = x
        y_var[:] = y

        x_var.units = "meters"
        x_var.standard_name = "projection_x_coordinate"

        y_var.units = "meters"
        y_var.standard_name = "projection_y_coordinate"

        # fix units of 'vx' and 'vy'
        nc.variables['vx'].units = "m / year"
        nc.variables['vy'].units = "m / year"


    # Compute and save the velocity magnitude
    if 'magnitude' not in nc.variables:

            vx = np.copy(nc.variables['vx'][:])
            vy = np.copy(nc.variables['vy'][:])

            v_magnitude = np.zeros_like(vx)

            v_magnitude = np.sqrt(vx**2 + vy**2)

            magnitude = nc.createVariable('v_magnitude', 'f8', ('time','y', 'x'),fill_value=0.0)
            #magnitude = nc.createVariable('v_magnitude', 'f8', ('x', 'y'))
            magnitude.units = "m / year"
            magnitude[0,:] = v_magnitude

    nc.close()



    subprocess.check_call('ncatted -O -a proj4,global,a,c,"+lon_0=0.0 +ellps=WGS84 +datum=WGS84 +lat_ts=-71.0 +proj=stere +x_0=0.0 +units=m +y_0=0.0 +lat_0=-90.0 +init=epsg:3031" '+fnout,shell=True)

    # make all variables double (some already are).
    subprocess.check_call('ncap2 -O -s "vx=double(vx);vy=double(vy);v_magnitude=double(v_magnitude);x=double(x);y=double(y)" '+
                          fnout+" "+fnout,shell=True)


if __name__ == "__main__":

  all_file_out = ""

  for filename in sorted(glob.glob(input_filenames)):

    filenameout = filename.split("/")[-1]
    filenameout = os.path.join(data_path,dataset+"_"+filenameout)
    filenameout = filenameout.replace(".nc","_input.nc")
    all_file_out+=filenameout+" "

    year = np.float(filename.split("/")[-1].split("_")[4])-1.5

    preprocess_ice_velocity(filename,filenameout,year)


  subprocess.check_call('ncrcat -v v_magnitude,lon,lat '+all_file_out+" "+output_filename,shell=True)

  # prepare the input file for cdo remapping
  # this step takes a while for high resolution data (i.e. 1km)
  pi.prepare_ncfile_for_cdo(output_filename)

  print " Rignot file",output_filename,"successfully preprocessed."


