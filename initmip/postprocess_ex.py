#!/usr/bin/env python
# Copyright (C) 2016 Andy Aschwanden (https://github.com/pism/pism-gris/blob/master/initMIP)
# Modified by Torsten Albrecht for initMIP Antarctica 2017

import os
import glob
import numpy as np
import json
import csv
import cf_units
try:
    import subprocess32 as sub
except:
    import subprocess as sub
    
from argparse import ArgumentParser
from netCDF4 import Dataset as CDF

#from initMIP.resources_ismip6 import *
import sys
sys.path.append('resources/')
from resources_ismip6 import *

## this hack is needed to import config.py from the project root
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if project_root not in sys.path: sys.path.append(project_root)
import config as cf; reload(cf)

dataset = "initmip"
data_path = os.path.join(cf.output_data_path, dataset)
resolution = cf.resolution

#python postprocessing_ex.py /p/tmp/albrecht/pism17/pismOut/initmip/initmip2200/initdata/pism_bedmap2_racmo_uplift_velrignot_lgmokill_15km.nc /p/tmp/albrecht/pism17/pismOut/initmip/initmip2200/results/extra_mbcotrl_15km_105yrs_latlon.nc -n 1 -e ctrl -r bil -t 16000 --id 1
#-w searise_grid_16000m_bil_weights.nc

# Set up the option parser
parser = ArgumentParser()
parser.description = "Script to make ISMIP6-conforming 2D time series."
parser.add_argument("INIT_FILE", nargs=1)
parser.add_argument("EXP_FILE", nargs=1)
parser.add_argument("-n", '--n_procs', dest="n_procs", type=int,
                    help='''number of cores/processors. default=4.''', default=4)
parser.add_argument("-e", "--experiment", dest="experiment",
                    choices=['ctrl', 'asmb', 'abmb'],
                    help="Experiment type", default='ctrl')
parser.add_argument('--id', dest="id", type=str,
                    help='''Experiemnt ID''', default='1')
parser.add_argument("-r", "--remap_method", dest="remap_method",
                    choices=['ycon', 'bil'],
                    help="Remapping method", default='con')
parser.add_argument("-t", "--target_resolution", dest="target_resolution", type=int,
                    choices=[1000, 8000, 16000],
                    help="Horizontal grid resolution", default=1000)
parser.add_argument("-w", "--override_weights_file",
                    dest="override_weights_file", action="store_true",
                    help="Override weights file", default=False)

options = parser.parse_args()
experiment = options.experiment
infile = options.EXP_FILE[0]
n_procs = options.n_procs
override_weights_file = options.override_weights_file
remap_method = options.remap_method
target_resolution = options.target_resolution
target_res_km = target_resolution / 1000

target_grid_filename = 'searise_grid_{}m.nc'.format(target_resolution)

# Need to get grid resolution from file
nc = CDF(infile, 'r')
pism_grid_dx = int(round(nc.variables['run_stats'].grid_dx_meters))
pism_vars_av = nc.variables.keys()
nc.close()
PISM_GRID_RES_ID = str(pism_grid_dx / 1000)
TARGET_GRID_RES_ID = str(target_resolution / 1000)
ID = options.id

IS = 'AIS'
GROUP = 'PIK'
MODEL = 'PISM' + ID + 'PAL'
EXP = experiment
TYPE = '_'.join([EXP, '0' + TARGET_GRID_RES_ID])
INIT = '_'.join(['init', '0' + TARGET_GRID_RES_ID])
project = '{IS}_{GROUP}_{MODEL}'.format(IS=IS, GROUP=GROUP, MODEL=MODEL)

pism_stats_vars = ['pism_config',
                   'run_stats']
pism_proj_vars = ['cell_area',
                  'mapping',
                  'lat',
                  'lat_bnds',
                  'lon',
                  'lon_bnds']

#from Thomas Kleiner
#ismip6_vars_file=open('INITMIP.json')
#initmip_vars_dict=json.loads(ismip6_vars_file.read())
#inimip_vars_dict=initmip_vars_dict["variables"]
#print inimip_vars_dict["lithk"]

ismip6_vars_dict = get_ismip6_vars_dict('resources/ismip6vars.csv', 2)
ismip6_to_pism_dict = dict((k, v.pism_name) for k, v in ismip6_vars_dict.iteritems())
pism_to_ismip6_dict = dict((v.pism_name, k) for k, v in ismip6_vars_dict.iteritems())

pism_copy_vars = [x for x in (ismip6_to_pism_dict.values() + pism_proj_vars)] #+ pism_stats_vars

mask_var = 'sftgif' 


def prepare_inputfile(infile,outfile):

  resfile = infile.replace("extra_","result_")

  print "Prepare PISM output for ismip6 compatibility and save as '{}'".format(outfile)
  ncks_cmd = ['ncks', '-O', '-4', '-L', '3',
             infile,
             outfile]
  sub.call(ncks_cmd)

  #lonlat_list = ['lat','lon', 'lat_bnds', 'lon_bnds', 'cell_area']
  #if all([x in pism_vars_av for x in lonlat_list]):
  print "  Add 'lon' and 'lat' coordinates and 'cell_area' variable"
  ncks_cmd = ['ncks', '-A', '-4', '-L', '3',
            '-v', ','.join(['lat','lon', 'lat_bnds', 'lon_bnds', 'cell_area']),
            resfile,
            outfile]
  sub.call(ncks_cmd)

  if 'mapping' not in pism_vars_av:
    print "  Add 'mapping' information"
    nc = CDF(outfile, 'a')
    mapping = nc.createVariable("mapping", 'c')
    mapping.ellipsoid = "WGS84"
    mapping.false_easting = 0.
    mapping.false_northing = 0.
    mapping.grid_mapping_name = "polar_stereographic"
    mapping.latitude_of_projection_origin = -90.
    mapping.standard_parallel = -71.
    mapping.straight_vertical_longitude_from_pole = 0.
    #except:
    #print "  Mapping seems to exist!"
    nc.close()

  print "  Add 'projection' information"
  ncatted_cmd = ["ncatted", '-O',
                "-a", '''proj4,global,o,c,+lon_0=0.0 +ellps=WGS84 +datum=WGS84 +lat_ts=-71.0 +proj=stere +x_0=0.0 +units=m +y_0=0.0 +lat_0=-90.0''',
                outfile]
  sub.call(ncatted_cmd)


  print "  Convert added variables to single precision" 
  ncap2_cmd = ['ncap2', '-O', '-s',
            '''"x=float(x);y=float(y);lon=float(lon);lat=float(lat);lat_bnds=float(lat_bnds);lon_bnds=float(lon_bnds);"''',
            outfile,
            outfile]
  sub.call(ncap2_cmd)

  print "  Add variable 'base'"
  ncap2_cmd = ['ncap2', '-O', '-s',
            'base = usurf - thk',
            outfile,
            outfile]
  sub.call(ncap2_cmd)
  ncatted_cmd = ["ncatted", '-O',
                "-a", '''long_name,base,o,c,ice lower surface elevation''',
                "-a", '''standard_name,base,o,c,base_altitude''',
                outfile]
  sub.call(ncatted_cmd)

  print "  Add variable 'ligroundf'"
  ncap2_cmd = 'ncap2 -O -s "ligroundf = -2e9 * discharge_flux / discharge_flux;" '
  ncap2_cmd += outfile+' '+outfile
  sub.check_call(ncap2_cmd,shell=True)
  ncatted_cmd = ["ncatted","-O",
                  "-a", '''standard_name,ligroundf,a,c,land_ice_specific_mass_flux_due_at_grounding_line''',
                  "-a", '''long_name,ligroundf,o,c,average ice flux across grounding line over reporting interval''',
                  "-a", '''comment,ligroundf,o,c,not available in PISM''',
                  "-a", '''units,ligroundf,o,c,kg m-2 s-1''',
                  outfile]
  sub.call(ncatted_cmd)

  print "Finished preprocessing and save to "+outfile


if __name__ == "__main__":


    project_dir = os.path.join(data_path, GROUP, MODEL, TYPE)
    if not os.path.exists(project_dir):
        os.makedirs(project_dir)

    init_dir = os.path.join(data_path, GROUP, MODEL, INIT)
    if not os.path.exists(init_dir):
        os.makedirs(init_dir)

    tmp_dir = os.path.join(data_path,'_'.join(['tmp', MODEL]))
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)

    tmp_filename = 'tmp_{}.nc'.format(EXP)
    tmp_file = os.path.join(data_path, tmp_dir, tmp_filename)
    try:
        os.remove(tmp_file)
    except OSError:
        pass

    infile_ismip6=infile.replace('.nc','_ismip6.nc')
    prepare_inputfile(infile,infile_ismip6)

    # Check if request variables are present
    nc = CDF(infile_ismip6, 'r')
    for m_var in pism_copy_vars:
        if m_var not in nc.variables:
            print("Requested variable '{}' missing".format(m_var))
    nc.close()
    print('Copy {} to {}'.format(infile_ismip6, tmp_file))
    cmd = ['ncks', '-O', '-d', 'time,1,',
           '-v', '{}'.format(','.join(pism_copy_vars)),
           infile_ismip6, tmp_file]
    sub.call(cmd)
    
    
    # Make the file ISMIP6 conforming
    make_spatial_vars_ismip6_conforming(tmp_file, ismip6_vars_dict)
    # Should be temporary until new runs
    ncatted_cmd = ["ncatted",
                   "-a", '''bounds,lat,o,c,lat_bnds''',
                   "-a", '''bounds,lon,o,c,lon_bnds''',
                   "-a", '''coordinates,lat_bnds,d,,''',
                   "-a", '''coordinates,lon_bnds,d,,''',
                   "-a", '''history,global,d,,''',
                   "-a", '''history_of_appended_files,global,d,,''',
                   tmp_file]
    sub.call(ncatted_cmd)


    out_filename = '{project}_{exp}.nc'.format(project=project, exp=EXP)
    out_file = os.path.join(tmp_dir, out_filename)
    try:
        os.remove(out_file)
    except OSError:
        pass

    #cdo remap
    if np.int(PISM_GRID_RES_ID) != np.int(target_res_km):

      print('Prepare the remap from {}km to {}km'.format(str(PISM_GRID_RES_ID),str(target_res_km)))

      # Create source grid definition file
      source_grid_filename = 'source_grid.nc'
      source_grid_file = os.path.join(tmp_dir, source_grid_filename)
      print('create source grid file {}'.format(source_grid_file))
      ncks_cmd = ['ncks', '-O', '-v', 'thk,mapping', infile_ismip6, source_grid_file]
      sub.call(ncks_cmd)

      #nc2cdo_cmd = ['nc2cdo.py', source_grid_file]
      nc2cdo_cmd = [os.path.join(os.path.dirname(os.path.dirname(cf.output_data_path)), 
                    'tools', 'nc2cdo.py'), source_grid_file]
      sub.call(nc2cdo_cmd)

    
      # If exist, remove target grid description file
      target_grid_file = os.path.join(tmp_dir, target_grid_filename)
      try:
        os.remove(target_grid_file)
      except OSError:
        pass

      # Create target grid description file
      print('create target grid file {}'.format(target_grid_file))
      create_searise_grid(IS, target_grid_file, target_resolution)
    
      # Generate weights if weights file does not exist yet
      cdo_weights_filename = 'searise_grid_{resolution}m_{method}_weights.nc'.format(resolution=target_resolution, method=remap_method)
      cdo_weights_file = os.path.join(tmp_dir, cdo_weights_filename)
      if (not os.path.isfile(cdo_weights_file)) or (override_weights_file is True):
        print('Generating CDO weights file {}'.format(cdo_weights_file))
        if n_procs > 1:
            cdo_cmd = ['cdo', '-P', '{}'.format(n_procs),
                       'gen{method},{grid}'.format(method=remap_method, grid=target_grid_file),
                       source_grid_file,
                       cdo_weights_file]
        else:
            cdo_cmd = ['cdo',
                       'gen{method},{grid}'.format(method=remap_method, grid=target_grid_file),
                       source_grid_file,
                       cdo_weights_file]            
        sub.call(cdo_cmd)
    
      # Remap to SeaRISE grid    
      print('Remapping to SeaRISE grid')
      if n_procs > 1:
        cdo_cmd = ['cdo', '-P', '{}'.format(n_procs),
                   'remap,{},{}'.format(target_grid_file, cdo_weights_file),
                   tmp_file,
                   out_file]
      else:
        cdo_cmd = ['cdo',
                   'remap,{},{}'.format(target_grid_file, cdo_weights_file),
                   tmp_file,
                   out_file]
      sub.call(cdo_cmd)
    else: #already on seaRISE grid
      cp_cmd = ['cp', tmp_file, out_file]
      sub.call(cp_cmd)

    #Convert to netcdf4 format
    ncks_cmd = ['ncks', '-O', '-4', '-L', '3',
                out_file,
                out_file]
    sub.call(ncks_cmd)

    
    # Adjust the time axis
    print('Adjusting time axis')
    adjust_time_axis(IS,out_file)

    for m_var in ismip6_vars_dict.keys():
        final_file = '{}/{}_{}_{}.nc'.format(project_dir, m_var, project, EXP)
        print('Finalizing variable {}'.format(m_var))
        # Generate file
        print('  Copying to file {}'.format(final_file))
        ncks_cmd = ['ncks', '-O', '-4', '-L', '3',
                    '-v', m_var ,
                    #'-v', ','.join([m_var,'lat','lon', 'lat_bnds', 'lon_bnds']),
                    out_file,
                    final_file]
        sub.call(ncks_cmd)
        # Add stats vars
        #print('  Adding config/stats variables')
        #ncks_cmd = ['ncks', '-A',
        #            '-v', ','.join(pism_stats_vars),
        #            tmp_file,
        #            final_file]
        #sub.call(ncks_cmd)
        # Add coordinate vars and mapping
        #print('  Adding coordinate and mapping variables')
        #ncks_cmd = ['ncks', '-A', '-v', 'x,y,mapping',
        #            target_grid_file,
        #            final_file]
        #sub.call(ncks_cmd)



        # flip signs for some fluxes to comply with arbitrary sign convention
        #if m_var in ('libmassbf', 'licalvf'):
        #    cmd = ['ncap2', '-O', '-s', '''"{var}={var}*-1;"'''.format(var=m_var),
        #                final_file,
        #                final_file]
        #    sub.call(cmd)
       
        if ismip6_vars_dict[m_var].do_mask == 1:

            print('  Mask ice free areas')
            # add mask variable
            cmd = ['ncks', '-A', '-v', '{var}'.format(var=mask_var),
                        out_file,
                        final_file]

            sub.call(cmd)
            # mask where mask==0
            cmd = 'ncap2 -O -s "where({maskvar}==0) {var}=-2e9;" '.format(maskvar=mask_var, var=m_var)
            cmd += final_file+' '+final_file
            sub.check_call(cmd,shell=True)
 

            cmd = ["ncatted", '-O',
                   "-a", '''_FillValue,{var},o,f,-2e9'''.format(var=m_var),
                   "-a", '''missing_value,{var},o,f,-2e9'''.format(var=m_var),
                   final_file]
            sub.call(cmd)

            # remove mask variable
            cmd = ['ncks', '-O', '-x', '-v', '{var}'.format(var=mask_var),
                    final_file,
                    final_file]
            sub.call(cmd)

        # Update attributes
        print('  Adjusting attributes')

        # remove lon lat variable
        var_list='lat_bnds,lon_bnds,lat,lon' #,time_bnds'
        rm_cmd = ['ncks', '-C','-O', '-x', '-v', '{var}'.format(var=var_list),
                  final_file,
                  final_file]
        sub.call(rm_cmd)

        ncap2_cmd = 'ncap2 -O -s "x=float(x);y=float(y);time=float(time);time_bounds=float(time_bounds);" '
        ncap2_cmd += final_file+' '+final_file
        sub.check_call(ncap2_cmd,shell=True)

        nc = CDF(final_file, 'a')
        try:
            nc_var = nc.variables[m_var]
            #nc_var.coordinates = 'lat lon'
            #nc_var.mapping = 'mapping'
            nc_var.standard_name = ismip6_vars_dict[m_var].standard_name
            nc_var.units = ismip6_vars_dict[m_var].units
            nc_var.pism_intent = None
            if ismip6_vars_dict[m_var].state == 1:
                nc_var.cell_methods = 'time: mean (interval: 5 year)' 
        except:
            pass
        nc.Conventions = 'CF-1.6'
        nc.institution = 'Potsdam Institute for Climate Impact Research (PIK), Germany'
        nc.contact = 'torsten.albrecht@pik-potsdam.de and matthias.mengel@pik-potsdam.de'
        nc.source = 'PISM (https://github.com/talbrecht/pism_pik; branch: pik_newdev_paleo_07; commit: 5d9d88e)'
        nc.title = 'ISMIP6 AIS InitMIP'
        nc.close()

        ncatted_cmd = ["ncatted","-hO",
                   "-a", '''nco_openmp_thread_number,global,d,,''',
                   "-a", '''command,global,d,,''',
                   "-a", '''history,global,d,,''',
                   "-a", '''history_of_appended_files,global,d,,''',
                   "-a", '''NCO,global,d,,''',
                   "-a", '''CDI,global,d,,''',
                   "-a", '''_NCProperties,global,d,,''', 
                   final_file]
        sub.call(ncatted_cmd)

        print('  Done finalizing variable {}'.format(m_var))

        if EXP in ('ctrl'):
            init_file = '{}/{}_{}_{}.nc'.format(init_dir, m_var, project, 'init')
            print('  Copying time 0 to file {}'.format(init_file))
            ncks_cmd = ['ncks', '-O', '-4', '-L', '3',
                        '-d', 'time,0',
                        '-v', m_var,
                        final_file,
                        init_file]
            sub.call(ncks_cmd)
