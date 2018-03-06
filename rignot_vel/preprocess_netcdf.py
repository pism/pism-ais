#!/usr/bin/env python

"""
matthias.mengel@pik, torsten.albrecht@pik
Download Rignot velocity data and save to (450m) netcdf file.

Preprocess the ~4GB Antarctic ice velocity dataset from NASA MEASURES project
http://nsidc.org/data/nsidc-0484.html

See infos https://nsidc.org/data/docs/measures/nsidc0484_rignot/
Currently not available without Earthdata login
"""


import numpy as np
import netCDF4 as NC
#import datetime
import os, sys
import subprocess

## this hack is needed to import config.py from the project root
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if project_root not in sys.path: sys.path.append(project_root)
import config as cf; reload(cf)
import pism_input.pism_input as pi; reload(pi)

dataset = "rignotvel"
data_path = os.path.join(cf.output_data_path, dataset)

if not os.path.exists(data_path): os.makedirs(data_path)

output_filename = os.path.join(data_path,dataset+"_450m_input.nc")

input_filename = "/p/projects/tumble/pism_input/Velocity/mouginot_rignot17/antarctica_ice_velocity_450m_v2.nc"
url = "https://daacdata.apps.nsidc.org/DATASETS/nsidc0484_MEASURES_antarc_vel_v02/"

#url = "ftp://anonymous@sidads.colorado.edu/pub/DATASETS/nsidc0484_MEASURES_antarc_vel_V01/900m/"
#url = " ftp://n5eil01u.ecs.nsidc.org/SAN/MEASURES/NSIDC-0484.001/1996.01.01/"
#input_filename = "antarctica_ice_velocity_900m.nc"


def run(commands):
    """Run a list of commands (or one command given as a string)."""
    if isinstance(commands, (list, tuple)):
        for cmd in commands:
            print "Running '%s'..." % cmd
            subprocess.call(cmd.split(' '))
    else:
        run([commands])


def preprocess_ice_velocity():

    commands = [#"wget -nc %s%s" % (url, input_filename), # NSIDC supports compression on demand!
                #"ncrename -d nx,x -d ny,y -O %s %s" % (input_filename, input_filename),
                "cp %s %s" % (input_filename,output_filename),
                "ncrename -v VX,vx -v VY,vy " + output_filename,
                #"ncks -O -v vx,vy,x,y,mapping %s %s" % (output_filename, output_filename)
                "ncks -O -v vx,vy,x,y,mapping,STDX,STDY %s %s" % (output_filename, output_filename)
                ]


    if not os.path.exists(output_filename):
        run(commands)

    nc = NC.Dataset(output_filename, 'a')


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


    y_var = nc.variables['y']
    #print y_var[:],np.flipud(y_var[:])
    y_var[:] = np.flipud(y_var[:])


    # Compute and save the velocity magnitude
    if 'magnitude' not in nc.variables:
            vx = nc.variables['vx'][:]
            vy = nc.variables['vy'][:]

            v_magnitude = np.zeros_like(vx)

            v_magnitude = np.sqrt(vx**2 + vy**2)

            magnitude = nc.createVariable('v_magnitude', 'f8', ('y', 'x'))
            #magnitude = nc.createVariable('v_magnitude', 'f8', ('x', 'y'))
            magnitude.units = "meter/year"

            magnitude[:] = v_magnitude


            sx = nc.variables['STDX'][:]
            sy = nc.variables['STDY'][:]
            v_deviation = np.zeros_like(vx)
            v_deviation = np.sqrt(sx**2 + sy**2)
            #v_deviation = 0.5*(sx+sy)/np.sqrt(vx**2 + vy**2)
            uncertainty = nc.createVariable('v_std', 'f8', ('y', 'x'))
            uncertainty.units = "meter/year"
            uncertainty[:] = v_deviation


    nc.close()

    subprocess.check_call('ncatted -O -a proj4,global,a,c,"+lon_0=0.0 +ellps=WGS84 +datum=WGS84 +lat_ts=-71.0 +proj=stere +x_0=0.0 +units=m +y_0=0.0 +lat_0=-90.0 +init=epsg:3031" '+output_filename,shell=True)

    # make all variables double (some already are).
    subprocess.check_call('ncap2 -O -s "vx=double(vx);vy=double(vy);v_magnitude=double(v_magnitude);x=double(x);y=double(y);v_std=double(v_std)" '+
                          output_filename+" "+output_filename,shell=True)

if __name__ == "__main__":

    preprocess_ice_velocity()

    # prepare the input file for cdo remapping
    # this step takes a while for high resolution data (i.e. 1km)
    pi.prepare_ncfile_for_cdo(output_filename)

    print " Rignot file",output_filename,"successfully preprocessed."


