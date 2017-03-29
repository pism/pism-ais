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

postprocessing:
python ../tools/fill_missing.py -v asmb,abmb -i 0.0 -e 10 -f initmip_15km.nc -o initmip_15km_filled.nc
