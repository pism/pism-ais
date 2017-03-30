Prepare data for initMIP Antarctica

Wiki: http://www.climate-cryosphere.org/wiki/index.php?title=InitMIP-Antarctica

Downloaded dBasalMelt and dSMB anomaly fields from 
ftp searise@cryoftp1.gsfc.nasa.gov initMIP directory /ISMIP6/initMIP/AIS

Password personal communication with Sophie Nowicki <sophie.nowicki@nasa.gov>

official email: ismip6 <ismip6@gmail.com>
cc: helene.seroussi@jpl.nasa.gov


merge: 
mkdir initmip_data
cp dSMB/smb_anomaly_1km.nc initmip_data/initmip_1km_input.nc
ncks -A -v abmb dBasalMelt/basal_melt_anomaly_1km.nc initmip_data/initmip_1km_input.nc

postprocessing (fill NaN with 0.0):
python ../tools/fill_missing.py -v asmb,abmb -i 0.0 -e 10 -f initmip_15km.nc -o initmip_15km_filled.nc

get background forcing (result of 100kyr constant forcing):
ncks -A -v effective_climatic_mass_balance,effective_ice_surface_temp,x,y,effective_shelf_base_temperature,effective_shelf_base_mass_flux /p/tmp/albrecht/pism17/pismOut/forcing/forcing2119_TPSO/results/result_constant_15km_100000yrs.nc pism_smb_f2119_const_15km.nc

#ncrename -O -v effective_climatic_mass_balance,climatic_mass_balance pism_smb_f2119_const_15km.nc
#ncrename -O -v effective_ice_surface_temp,ice_surface_temp pism_smb_f2119_const_15km.nc
#ncrename -O -v effective_shelf_base_mass_flux,shelfbmassflux pism_smb_f2119_const_15km.nc
#ncrename -O -v effective_shelf_base_temperature,shelfbtemp pism_smb_f2119_const_15km.nc

add step forcing
python create_anomalies.py --force_file initmip_15km_filled.nc --background_file pism_smb_f2119_const_15km.nc initmip_15km_forcing.nc

