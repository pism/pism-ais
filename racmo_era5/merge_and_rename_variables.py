"""
matthias.mengel@pik, torsten.albrecht@pik, ronja.reese@pik, moritz.kreuzer@pik
"""

import os, sys
import subprocess
import netCDF4 as nc
import numpy as np
import shutil

import xarray as xr
import cftime

## this hack is needed to import config.py from the project root
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if project_root not in sys.path: sys.path.append(project_root)
import config as cf;# reload(cf)
import pism_input.pism_input as pi;# reload(pi)

dataset = "racmo_era5"
data_reference = "Van Wessem et al., 2014. Improved representation of East Antarctica surface mass balance in a regional climate model. J. Glac., 60(222), 761-770, doi: 10.3189/2014JoG14J051."
data_reference2 = "Van Wessem, Jan Melchior, Willem Jan Van De Berg, Brice PY Noel, Erik Van Meijgaard, Charles Amory, Gerit Birnbaum, Constantijn L. Jakobs et al. Modelling the climate and surface mass balance of polar ice sheets using racmo2: Part 2: Antarctica (1979-2016). Cryosphere 12, no. 4 (2018): 1479-1498."

data_path = os.path.join(cf.output_data_path, dataset)

if not os.path.exists(data_path): os.makedirs(data_path)

tmp_file = os.path.join(data_path, dataset+"_tmp.nc")
output_file = os.path.join(data_path, dataset+"_input.nc")

source_file = {"t2m": os.path.join(cf.racmo_era5_data_path,"t2m_monthlyA_ANT27_ERA5-3H_RACMO2.3p2_197901_202109.nc"),
                "smb": os.path.join(cf.racmo_era5_data_path,"smb_monthlyS_ANT27_ERA5-3H_RACMO2.3p2_197901_202109.nc"),
                "evap": os.path.join(cf.racmo_era5_data_path,"evap_monthlyS_ANT27_ERA5-3H_RACMO2.3p2_197901_202109.nc"),
                "precip": os.path.join(cf.racmo_era5_data_path,"precip_monthlyS_ANT27_ERA5-3H_RACMO2.3p2_197901_202109.nc")
                }

process_file = {var:os.path.join(data_path, dataset+"_"+var+".nc") for var in ["t2m","smb","evap","precip"]}

for var,fl in process_file.items():
    try:
        os.remove(fl)
    except OSError:
        pass

for var in ["t2m","smb","evap","precip"]:

    subprocess.check_call("ncks -A -v "+var+",lon,lat "+source_file[var]+" "+process_file[var],shell=True)

    subprocess.check_call("ncrename -d rlat,y -d rlon,x -v rlat,y -v rlon,x "+process_file[var],shell=True)
    
    subprocess.check_call("ncatted -O -a reference,"+var+",c,c,'" + data_reference +
                            "' "+process_file[var],shell=True)

    subprocess.check_call("ncatted -O -a grid_mapping,"+var+",d,, "+process_file[var],shell=True)

    subprocess.check_call('ncatted -O -a proj4,global,o,c,"+lon_0=10.0 +ellps=WGS84 +datum=WGS84 +lat_ts=-71.0 +proj=stere +x_0=0.0 +units=m +y_0=0.0 +lat_0=-90.0" '+process_file[var], shell=True)

    subprocess.check_call('ncwa -O -a height '+process_file[var]+" "+process_file[var], shell=True) #delete the height dimension!
    subprocess.check_call('ncks -O -x -v height '+process_file[var]+" "+process_file[var], shell=True) #delete the height dimension!

    subprocess.check_call('ncks -O -C -x -v nblock1 '+process_file[var]+" "+process_file[var], shell=True)
    subprocess.check_call('ncks -O -C -x -v nblock2 '+process_file[var]+" "+process_file[var], shell=True)

# Converting SMB units: Not needed probably.
# subprocess.check_call("ncap2 -O -s 'smb=smb*1000.0/910.0' "+tmp_file+" "+tmp_file,shell=True)

# Set SMB on ocean to missing, so it does not disturb the remapping at the ice sheet boundary.
# subprocess.check_call('cdo -O setrtomiss,-9999,0 '+process_file["smb"]+" "+smb_omask_file, shell=True)

# smb_omask_file = os.path.join(data_path, dataset+"_smb_omask.nc")
# shutil.copyfile(process_file["smb"], smb_omask_file)
# smb_msk_data = nc.Dataset(smb_omask_file,"a")
# # be aware: this is only valid for this dataset and timestep zero.
# mask = smb_msk_data.variables["smb"][0,:,:] < -0.009
# shape = smb_msk_data.variables["smb"].shape
# smb_masked = np.ma.masked_array(smb_msk_data.variables["smb"][:],
#                                 mask=np.tile(mask,(shape[0],1,1)))
# smb_msk_data.variables["smb"][:] = smb_masked
# smb_msk_data.close()

# process_file["smb"] = smb_omask_file

merge_these_files = " ".join([process_file[var] for var in ["t2m","smb","evap","precip"]])

subprocess.check_call('cdo -O merge '+merge_these_files+" "+tmp_file, shell=True)

# make all variables double (some already are).
subprocess.check_call("ncap2 -O -s 't2m=double(t2m);smb=double(smb);evap=double(evap);precip=double(precip)' "+
                      tmp_file+" "+tmp_file,shell=True)

subprocess.check_call("ncrename -v t2m,ice_surface_temp -O "+tmp_file+" "+tmp_file,shell=True)

# Fill the missing SMB field over ocean with the proxy precip - evaporation
ncf = nc.Dataset(tmp_file,"a")
smb = ncf.variables["smb"][:]
# assume: smb[0,0,0] has ocean fill value
mask_ocean = (smb[0,:,:] == smb[0,0,0])
mask_ocean_expanded = np.tile(mask_ocean,(smb.shape[0],1,1))
## be aware: this masking is only valid for this dataset and timestep zero.
#mask_ocean = smb[0,:,:] < -0.009
#mask_ocean_expanded = np.tile(mask_ocean,(smb.shape[0],1,1))
# a proxy for smb
smb_over_ocean = ncf.variables["precip"][:] + ncf.variables["evap"][:]
smb[mask_ocean_expanded] = smb_over_ocean[mask_ocean_expanded]
ncf.variables["smb"][:] = smb
ncf.smb_comment = "SMB is approximated by precip-evap over the ocean."
ncf.description = "RACMO2.3p2 data (ANT27/2) forced by ERA-5 provide yearly mean air temperature (t2m) and surface mass balance (smb) for the years 1979-2021, here averaged over CMIP5 period 1985-2005"
ncf.link= "https://www.dropbox.com/sh/2i21p9p1abbwdb4/AADwdFNztZ4b5aJPWsCENdgEa?dl=0"
ncf.overview = "https://www.projects.science.uu.nl/iceclimate/models/antarctica.php"
ncf.reference = data_reference
ncf.reference2 = data_reference2

ncf.close()

subprocess.check_call('ncatted -a units,smb,o,c,"kg m-2" '+tmp_file,shell=True)
subprocess.check_call("ncrename -v smb,climatic_mass_balance -O "+tmp_file+" "+tmp_file,shell=True)
subprocess.check_call('ncatted -a units,precip,o,c,"kg m-2" '+tmp_file,shell=True)
subprocess.check_call("ncrename -v precip,precipitation -O "+tmp_file+" "+tmp_file,shell=True)
subprocess.check_call("ncap2 -O -s 'air_temp=ice_surface_temp' "+tmp_file+" "+tmp_file,shell=True)


### convert mass units to be time independent  (kg/m2 -> kg/(m2 year))
ds = xr.open_dataset(tmp_file, use_cftime=True)

time_daysinyear = [366 if (t%4 == 0) else 365 for t in ds.time.dt.year]
conversion_per_year = time_daysinyear / ds.time.dt.days_in_month

for v in ['climatic_mass_balance', 'precipitation', 'evap']:
    assert ds[v].attrs['units'] == 'kg m-2', f"units of var {var} is '{ds[v].attrs['units']}' instead of 'kg m-2'"
    attrs = ds[v].attrs
    ds[v] = ds[v] * conversion_per_year
    ds[v].attrs.update(attrs)
    ds[v].attrs['units'] = 'kg m-2 year-1'

ds.to_netcdf(output_file, unlimited_dims = ['time'])



# prepare the input file for cdo remapping
# this step takes a while for high resolution data (i.e. 1km)
# pi.prepare_ncfile_for_cdo(output_file)

subprocess.check_call(f"rm {tmp_file}", shell=True)

print(" RACMO file",output_file,"successfully preprocessed.")
