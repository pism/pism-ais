#!/bin/bash

#Prepare data for initMIP Antarctica
#Wiki: http://www.climate-cryosphere.org/wiki/index.php?title=InitMIP-Antarctica

#Downloaded dBasalMelt and dSMB anomaly fields from 
#ftp searise@cryoftp1.gsfc.nasa.gov initMIP directory /ISMIP6/initMIP/AIS

#Password personal communication with Sophie Nowicki <sophie.nowicki@nasa.gov>

#official email: ismip6 <ismip6@gmail.com>
#cc: helene.seroussi@jpl.nasa.gov

grid=15km
gridfile=../cdo_remapgrids/pism_${grid}.nc
PD_pism_file_orig=/p/tmp/albrecht/pism17/pismOut/forcing/forcing2119_TPSO/results/result_constant_15km_100000yrs.nc
PD_pism_file=pism_smb_f2119_const_15km.nc
IM_outfile=initmip_${grid}

#merge: 
mkdir initmip_data
cp dSMB/smb_anomaly_1km.nc initmip_data/initmip_1km_input.nc
#ncks -A -v abmb dBasalMelt/basal_melt_anomaly_1km.nc initmip_data/initmip_1km_input.nc

#do the remap
#remapcony for asmb and remapnn for abmb
#python remap.py
#sbatch cdo_remap.sh

#postprocessing 