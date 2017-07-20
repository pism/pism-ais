#!/bin/bash

#Prepare data for initMIP Antarctica
#Wiki: http://www.climate-cryosphere.org/wiki/index.php?title=InitMIP-Antarctica

#Downloaded dBasalMelt and dSMB anomaly fields from 
#ftp searise@cryoftp1.gsfc.nasa.gov initMIP directory /ISMIP6/initMIP/AIS

#Password personal communication with Sophie Nowicki <sophie.nowicki@nasa.gov>

#official email: ismip6 <ismip6@gmail.com>
#cc: helene.seroussi@jpl.nasa.gov


# 1)
# If PISM init file is on different grid than the iniMIP standard grid (1,2,4,8,16,32 km resolution),
# you can rempap the anomaly forcing fields to given native PISM grid and run the experiments...
# Results must be rempapped back to iniMIP grid.

grid=15km
gridfile=../cdo_remapgrids/pism_${grid}.nc
PD_pism_file_orig=/p/tmp/albrecht/pism17/pismOut/forcing/forcing2119_TPSO/results/result_constant_15km_100000yrs.nc
PD_pism_file=pism_smb_f2119_const_15km.nc
IM_outfile=initmip_${grid}

#merge forcing fields to one file: 
mkdir initmip_data
cp dSMB/smb_anomaly_1km.nc initmip_data/initmip_1km_input.nc
ncks -A -v abmb dBasalMelt/basal_melt_anomaly_1km.nc initmip_data/initmip_1km_input.nc

#now run the remap!!!

#postprocessing
'''
#SMB
#fill NaN with 0.0:
#python ../tools/fill_missing.py -v asmb,abmb -i 0.0 -e 10 -f ${IM_outfile}.nc -o ${IM_outfile}_filled.nc
python ../tools/fill_missing.py -v asmb -i 0.0 -e 10 -f ${IM_outfile}.nc -o ${IM_outfile}_filled.nc

#BMB
cdo remapnn,../cdo_remapgrids/pism_${grid}.nc dBasalMelt/basal_melt_anomaly_1km.nc ${IM_outfile}_cdonearest.nc
ncap2 -O -s "abmb=double(abmb)" ${IM_outfile}_cdonearest.nc ${IM_outfile}_cdonearest.nc
ncks -A -v abmb ${IM_outfile}_cdonearest.nc ${IM_outfile}_filled.nc

#init
#get background forcing (result of 100kyr constant forcing):
ncks -A -v effective_climatic_mass_balance,effective_ice_surface_temp,x,y,effective_shelf_base_temperature,effective_shelf_base_mass_flux ${PD_pism_file_orig} ${PD_pism_file}

#ncrename -O -v effective_climatic_mass_balance,climatic_mass_balance pism_smb_f2119_const_15km.nc
#ncrename -O -v effective_ice_surface_temp,ice_surface_temp pism_smb_f2119_const_15km.nc
#ncrename -O -v effective_shelf_base_mass_flux,shelfbmassflux pism_smb_f2119_const_15km.nc
#ncrename -O -v effective_shelf_base_temperature,shelfbtemp pism_smb_f2119_const_15km.nc

#add step forcing
python create_anomalies.py --force_file ${IM_outfile}_filled.nc --background_file ${PD_pism_file} ${IM_outfile}_forcing.nc

#add x and y and mapping
ncks -A -v x,y,mapping ${gridfile} ${IM_outfile}_forcing.nc

#change units
ncatted -a units,climatic_mass_balance,o,c,"kg m-2 year-1" -a standard_name,climatic_mass_balance,o,c,"land_ice_surface_specific_mass_balance" -a grid_mapping,climatic_mass_balance,o,c,"mapping" -a grid_mapping,ice_surface_temp,o,c,"mapping" -a grid_mapping,shelfbmassflux,o,c,"mapping" -a grid_mapping,shelfbtemp,o,c,"mapping" initmip_15km_forcing.nc

#control run
ncks -O -d time,0 initmip_15km_forcing.nc initmip_15km_control.nc
'''