#!/usr/bin/env python
# Copyright (C) 2016 Andy Aschwanden (https://github.com/pism/pism-gris/blob/master/initMIP)
# Modified by Torsten Albrecht for initMIP Antarctica 2017

import os
import numpy as np
import csv
import cf_units
try:
    import subprocess32 as sub
except:
    import subprocess as sub
    
from argparse import ArgumentParser
from netCDF4 import Dataset as CDF

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


# Set up the option parser
parser = ArgumentParser()
parser.description = "Script to make ISMIP6-conforming scalar time series."
parser.add_argument("EXP_FILE", nargs=1)
parser.add_argument("-e", "--experiment", dest="experiment",
                    choices=['ctrl', 'asmb', 'abmb'],
                    help="Output size type", default='ctrl')
parser.add_argument('--id', dest="id", type=str,
                    help='''Experiemnt ID''', default='1')
parser.add_argument("-t", "--target_resolution", dest="target_resolution", type=int,
                    choices=[1000, 8000, 16000],
                    help="Horizontal grid resolution", default=8000)

options = parser.parse_args()
experiment = options.experiment
infile = options.EXP_FILE[0]
target_resolution = options.target_resolution

# Need to get grid resolution from file
nc = CDF(infile, 'r')
pism_grid_dx = int(round(nc.variables['run_stats'].grid_dx_meters))
nc.close()
PISM_GRID_RES_ID = str(pism_grid_dx / 1000)
TARGET_GRID_RES_ID = str(target_resolution / 1000)
ID = options.id

IS = 'AIS'
GROUP = 'PIK'
#MODEL = 'PISM' + '_' + PISM_GRID_RES_ID + 'km_' + ID
MODEL = 'PISM' + ID + 'PAL'
EXP = experiment
TYPE = '_'.join([EXP, str(TARGET_GRID_RES_ID).zfill(2)])
INIT = '_'.join(['init', str(TARGET_GRID_RES_ID).zfill(2)])
project = '{IS}_{GROUP}_{MODEL}'.format(IS=IS, GROUP=GROUP, MODEL=MODEL)
pism_stats_vars = ['pism_config',
                   'run_stats']


ismip6_vars_dict = get_ismip6_vars_dict('resources/ismip6vars.csv', 1)
ismip6_to_pism_dict = dict((k, v.pism_name) for k, v in ismip6_vars_dict.iteritems())
pism_to_ismip6_dict = dict((v.pism_name, k) for k, v in ismip6_vars_dict.iteritems())

pism_copy_vars = [x for x in (ismip6_to_pism_dict.values())] #+ pism_stats_vars)]

if __name__ == "__main__":

    project_dir = os.path.join(data_path, GROUP, MODEL, TYPE)
    if not os.path.exists(project_dir):
        os.makedirs(project_dir)

    init_dir = os.path.join(data_path, GROUP, MODEL, INIT)
    if not os.path.exists(init_dir):
        os.makedirs(init_dir)
    
    out_filename = 'scalar_{project}_{exp}.nc'.format(project=project, exp=EXP)
    out_file = os.path.join(project_dir, out_filename)
    try:
        os.remove(out_file)
    except OSError:
        pass


    infile_ismip6=infile.replace('.nc','_ismip6.nc')

    print "Prepare PISM output for ismip6 compatibility and save as '{}'".format(infile_ismip6)
    ncks_cmd = ['ncks', '-O', '-4', '-L', '3',
             infile,
             infile_ismip6]
    sub.call(ncks_cmd)

    #"total over ice domain of top surface ice mass flux = sub_shelf_ice_flux + grounded_basal_ice_flux" ;
    print "  Add variable 'tendlibmassbf' "
    ncap2_cmd = ['ncap2', '-O', '-s',
            'tendlibmassbf = sub_shelf_ice_flux + grounded_basal_ice_flux;',
            infile_ismip6,
            infile_ismip6]
    sub.call(ncap2_cmd)

    # Check if request variables are present
    nc = CDF(infile_ismip6, 'r')
    for m_var in pism_copy_vars:
        if m_var not in nc.variables:
            print("Requested variable '{}' missing".format(m_var))

    nc.close()
    print('  Removing times < 0 in file {}'.format(out_file))
    cmd = ['ncks', '-O',
           '-d', 'time,4,-1',
           '-v', '{}'.format(','.join(pism_copy_vars)),
           infile_ismip6, out_file]
    sub.call(cmd)

    print "  Convert added variables to single precision" 
    for m_var in pism_copy_vars:

       ncap2_cmd = 'ncap2 -O -s "{}=float({});time=float(time);time_bounds=float(time_bounds)" '.format(m_var,m_var)
       ncap2_cmd += out_file+' '+out_file
       sub.check_call(ncap2_cmd,shell=True)

       ncatted_cmd = ["ncatted", '-O',
                   "-a", '''_FillValue,{var},o,f,-2e9'''.format(var=m_var),
                   "-a", '''missing_value,{var},o,f,-2e9'''.format(var=m_var),
                   out_file]
       sub.call(ncatted_cmd)


    # Adjust the time axis
    print('Adjusting time axis')
    adjust_time_axis(IS,out_file)
    make_scalar_vars_ismip6_conforming(out_file, ismip6_vars_dict)


    # Update attributes
    print('Adjusting attributes')
    nc = CDF(out_file, 'a')
    nc.Conventions = 'CF-1.6'
    nc.institution = 'Potsdam Institute for Climate Impact Research (PIK), Germany'
    nc.contact = 'torsten.albrecht@pik-potsdam.de and matthias.mengel@pik-potsdam.de'
    nc.source = 'PISM (https://github.com/talbrecht/pism_pik; branch: pik_newdev_paleo_07; commit: 9ae1674'
    #del nc.variables["pism_config"]
    #del nc.variables["run_stats"]
    #nc.variables["pism_config"] = None
    #nc.variables["run_stats"] = None
    nc.close()

    # remove mask variable
    #cmd = ['ncks', '-O', '-x', '-v', 'pism_config,run_stats',
    #        out_file,
    #        out_file]
    #sub.call(cmd)

    ncatted_cmd = ["ncatted","-hO",
                   "-a", '''nco_openmp_thread_number,global,d,,''',
                   "-a", '''command,global,d,,''',
                   "-a", '''history,global,d,,''',
                   "-a", '''history_of_appended_files,global,d,,''',
                   "-a", '''NCO,global,d,,''',
                   "-a", '''CDI,global,d,,''',
                   "-a", '''_NCProperties,global,d,,''',
                   out_file]
    sub.call(ncatted_cmd)


    print('Finished processing scalar file {}'.format(out_file))

    if EXP in ('ctrl'):
        init_file = '{}/scalar_{}_{}.nc'.format(init_dir, project, 'init')
        print('  Copying time 0 to file {}'.format(init_file))
        ncks_cmd = ['ncks', '-O', '-4', '-L', '3',
                    '-d', 'time,1',
                    out_file,
                    init_file]
        sub.call(ncks_cmd)
    
