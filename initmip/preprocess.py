#!/usr/bin/env python

# see also for Greenland: https://github.com/pism/pism-gris/blob/master/initMIP/prepare_anomalies.sh

#Prepare data for initMIP Antarctica
#Wiki: http://www.climate-cryosphere.org/wiki/index.php?title=InitMIP-Antarctica

#Downloaded dBasalMelt and dSMB anomaly fields from 
#ftp searise@cryoftp1.gsfc.nasa.gov initMIP directory /ISMIP6/initMIP/AIS

#Password personal communication with Sophie Nowicki <sophie.nowicki@nasa.gov>

#official email: ismip6 <ismip6@gmail.com>
#cc: helene.seroussi@jpl.nasa.gov

import os, sys
try:
    import subprocess32 as sub
except:
    import subprocess as sub
from netCDF4 import Dataset as CDF
#import create_anomalies

## this hack is needed to import config.py from the project root
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if project_root not in sys.path: sys.path.append(project_root)
import config as cf; reload(cf)

dataset = "initmip"
data_path = os.path.join(cf.output_data_path, dataset)
resolution = cf.resolution

data_resolution = 16 # or: 8,1
IM_filename = os.path.join(data_path,'initmip_'+str(data_resolution)+'km_input.nc')
IM_data_path = cf.initmip_data_path

PD_pism_out = cf.initmip_pism_out
try:
  pism_experiment=initmip_pism_out.split("/")[-3].split("_")[0] #specific naming 
except:
  pism_experiment='climate'
PD_pism_climate = os.path.join(data_path,'pism_'+pism_experiment+'_'+str(resolution)+'km.nc')
#PD_pism_climate = os.path.join(data_path,'pism_f2308_'+str(resolution)+'km.nc')

final_filename = os.path.join(data_path,'initmip_'+str(data_resolution)+'km_forcing.nc')
final_filename_ctrl = os.path.join(data_path,'initmip_'+str(data_resolution)+'km_control.nc')

# if anomaly fields data is not yet there in one file, create it
if not os.path.isfile(IM_filename):
  if not os.path.exists(data_path):
    os.makedirs(data_path)
  dsmb_file = os.path.join(IM_data_path,'dSMB/smb_anomaly_'+str(data_resolution)+'km.nc')
  dbmb_file = os.path.join(IM_data_path,'dBasalMelt/basal_melt_anomaly_'+str(data_resolution)+'km.nc')

  if os.path.isfile(dsmb_file) and os.path.isfile(dbmb_file):

    cmd_cp = ['cp', dsmb_file, IM_filename]
    sub.call(cmd_cp)

    cmd_ncks = ['ncks', '-A', '-v' , 'abmb', 
             dbmb_file, IM_filename]

    sub.call(cmd_ncks)
  else:
    print 'Anomaly data for '+str(data_resolution)+'km are not available, try 1km and remap'
    #FIXME: if needed at all


# get PISM initial climate forcing
var_list = 'x,y,effective_climatic_mass_balance,effective_ice_surface_temp,effective_shelf_base_temperature,effective_shelf_base_mass_flux'
cmd_ncks = ['ncks', '-A', '-v' , var_list,
             PD_pism_out, PD_pism_climate]
sub.call(cmd_ncks)

cmd_ncre = ['ncrename', '-O', 
             '-v' , 'effective_climatic_mass_balance,climatic_mass_balance',
             '-v' , 'effective_ice_surface_temp,ice_surface_temp',
             '-v' , 'effective_shelf_base_mass_flux,shelfbmassflux',
             '-v' , 'effective_shelf_base_temperature,shelfbtemp',
             PD_pism_climate]
#sub.call(cmd_ncre)

#add step forcing
ca_cmd = ['python','create_anomalies.py', '--force_file', IM_filename,
         '--background_file', PD_pism_climate, 
         '{}'.format(final_filename)]
#print ca_cmd
sub.call(ca_cmd)

#add x and y and mapping
ncks_cmd = ['ncks', '-A', '-v' , 'x,y',
           IM_filename, final_filename]
sub.call(ncks_cmd)

nc = CDF(final_filename, 'a')
try:
    mapping = nc.createVariable("mapping", 'c')
    mapping.ellipsoid = "WGS84"
    mapping.false_easting = 0.
    mapping.false_northing = 0.
    mapping.grid_mapping_name = "polar_stereographic"
    mapping.latitude_of_projection_origin = -90.
    mapping.standard_parallel = -71.
    mapping.straight_vertical_longitude_from_pole = 0.
except:
    print "  Mapping seems to exist!"
nc.close()

#change units
ncatted_cmd = ['ncatted', 
               '-a', '''units,climatic_mass_balance,o,c,kg m-2 year-1''',
               '-a', '''standard_name,climatic_mass_balance,o,c,"land_ice_surface_specific_mass_balance"''', 
               '-a', '''grid_mapping,climatic_mass_balance,o,c,"mapping"''', 
               '-a', '''grid_mapping,ice_surface_temp,o,c,"mapping"''', 
               '-a', '''grid_mapping,shelfbmassflux,o,c,"mapping"''', 
               '-a', '''grid_mapping,shelfbtemp,o,c,"mapping"''',
               '{}'.format(final_filename)]
sub.call(ncatted_cmd)


#control run
cmd_ncks = ['ncks', '-O', '-d' , 'time,0',
            final_filename, final_filename_ctrl]
sub.call(cmd_ncks)

