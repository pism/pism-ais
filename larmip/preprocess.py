#!/usr/bin/env python

#Prepare data for LARMIP Antarctica
#Wiki: http://www.climate-cryosphere.org/wiki/index.php?title=InitMIP-Antarctica

#Downloaded LARMIP regions from 
larmip_link = "https://www.pik-potsdam.de/research/earth-system-analysis/models/larmip/larmip-regions"

import os, sys
try:
    import subprocess32 as sub
except:
    import subprocess as sub
from netCDF4 import Dataset as CDF

## this hack is needed to import config.py from the project root
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if project_root not in sys.path: sys.path.append(project_root)
import config as cf; reload(cf)

dataset = "larmip"
data_path = os.path.join(cf.output_data_path, dataset)

data_resolution = 16 # or: 8,1
LM_filename = os.path.join(data_path,'LARMIP_regions_initMIPgrid_'+str(data_resolution)+'.nc')
vername="1p0"
#vername="0p7"

PD_pism_out = cf.initmip_pism_out
try:
  pism_experiment=PD_pism_out.split("/")[-3].split("_")[0] #specific naming
except:
  pism_experiment='climate'
PD_pism_climate = os.path.join(data_path,'pism_'+pism_experiment+'_'+str(data_resolution)+'km.nc')

final_filename_ctrl = os.path.join(data_path,'larmip_'+str(data_resolution)+'km_control'+vername+'.nc')

# if LARMIP region is not yet there in one file, download it
if not os.path.isfile(LM_filename):
  if not os.path.exists(data_path):
    os.makedirs(data_path)

  larmip_file = os.path.join(larmip_link,'larmip-regions-on-initmip-grid-'+str(data_resolution)+'km-resolution/at_download/file')
  cmd_wget = ['wget', larmip_file,'-O', LM_filename]
  sub.call(cmd_wget)

if not os.path.isfile(PD_pism_climate):
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
  sub.call(cmd_ncre)

PD_ocean_file = os.path.join(cf.output_data_path,"schmidtko/schmidtko_"+str(data_resolution)+"km_means_11+12.nc" )

for melt in [1,2,4,8,16,32]:
#for melt in [2]:
  for reg in xrange(6):
    print melt,reg
    final_filename = os.path.join(data_path,'larmip_'+str(data_resolution)+'km_forcing'+vername+'_reg'+str(reg)+'_m'+str(melt)+'.nc')
    if reg==5:
      final_filename = os.path.join(data_path,'larmip_'+str(data_resolution)+'km_forcing'+vername+'_all_m'+str(melt)+'.nc')

    if True:
    #if not os.path.isfile(final_filename):

      #add step forcing
      ca_cmd = ['python','create_anomalies.py','--background_file', PD_pism_climate,
                '--region_file', LM_filename,'--region', str(reg),'--meltrate', str(melt),
                '{}'.format(final_filename)]
      #print ca_cmd
      sub.call(ca_cmd)

      # include schmidtko ocean data
      cmd_ncks = ['ncks', '-A', '-v' , 'salinity_ocean,theta_ocean,basins',
                  PD_ocean_file, final_filename]
      sub.call(cmd_ncks)


      cmd_ncap2 = ['ncap2','-O','-s',"x=float(x);y=float(y)",
                  final_filename, final_filename]
      sub.call(cmd_ncap2)

      #add x and y and mapping
      ncks_cmd = ['ncks', '-A', '-v' , 'x,y',
                  LM_filename, final_filename]
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

print final_filename, "created."
