
"""
matthias.mengel@pik, torsten.albrecht@pik, ronja.reese@pik
"""

import os, sys
import subprocess
import netCDF4 as nc
import numpy as np
import shutil
## this hack is needed to import config.py from the project root
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if project_root not in sys.path: sys.path.append(project_root)
import config as cf; reload(cf)
import pism_input.pism_input as pi; reload(pi)

dataset = "racmo_hadcm3"
data_path = os.path.join(cf.output_data_path, dataset)

if not os.path.exists(data_path): os.makedirs(data_path)

output_file = os.path.join(data_path, dataset+"_input_c20.nc")

source_file = {"tskin": os.path.join(cf.racmo_hadcm3_data_path,"HadCM3_c20_I2S_tskin_Y.nc"),
                "smb": os.path.join(cf.racmo_hadcm3_data_path,"HadCM3_c20_I2S_smb_Y.nc")}
#                "evap": os.path.join(cf.racmo_wessem_data_path,"evap_RACMO2.3p2_yearly_ANT27_1979_2016.nc"),
#                "precip": os.path.join(cf.racmo_wessem_data_path,"precip_RACMO2.3p2_yearly_ANT27_1979_2016.nc")
#                }

#process_file = {var:os.path.join(data_path, dataset+"_"+var+".nc") for var in ["t2m","smb","evap","precip"]}
process_file = {var:os.path.join(data_path, dataset+"_"+var+"_c20.nc") for var in ["tskin","smb"]}

for var,fl in process_file.iteritems():
    try:
        os.remove(fl)
    except OSError:
        pass

#for var in ["t2m","smb","evap","precip"]:
for var in ["tskin","smb"]:

    subprocess.check_call("ncks -A -v "+var+",lon,lat "+source_file[var]+" "+process_file[var],shell=True)

    subprocess.check_call("ncrename -d rlat,y -d rlon,x -v rlat,y -v rlon,x "+process_file[var],shell=True)

    subprocess.check_call("ncatted -O -a grid_mapping,"+var+",d,, "+process_file[var],shell=True)

    subprocess.check_call('ncatted -O -a proj4,global,o,c,"+lon_0=10.0 +ellps=WGS84 +datum=WGS84 +lat_ts=-71.0 +proj=stere +x_0=0.0 +units=m +y_0=0.0 +lat_0=-90.0" '+process_file[var], shell=True)

    subprocess.check_call('ncwa -O -a height '+process_file[var]+" "+process_file[var], shell=True) #delete the height dimension!
    subprocess.check_call('ncks -O -x -v height '+process_file[var]+" "+process_file[var], shell=True) #delete the height dimension!

    subprocess.check_call('ncks -O -C -x -v nblock1 '+process_file[var]+" "+process_file[var], shell=True)
    subprocess.check_call('ncks -O -C -x -v nblock2 '+process_file[var]+" "+process_file[var], shell=True)

# Converting SMB units: Not needed probably.
# subprocess.check_call("ncap2 -O -s 'smb=smb*1000.0/910.0' "+output_file+" "+output_file,shell=True)

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

#merge_these_files = " ".join([process_file[var] for var in ["t2m","smb","evap","precip"]])
merge_these_files = " ".join([process_file[var] for var in ["tskin","smb"]])

subprocess.check_call('cdo -O merge '+merge_these_files+" "+output_file, shell=True)

# make all variables double (some already are).
#subprocess.check_call("ncap2 -O -s 't2m=double(t2m);smb=double(smb);evap=double(evap);precip=double(precip)' "+
#                      output_file+" "+output_file,shell=True)

#subprocess.check_call("ncrename -v t2m,ice_surface_temp -O "+output_file+" "+output_file,shell=True)
subprocess.check_call("ncrename -v tskin,ice_surface_temp -O "+output_file+" "+output_file,shell=True)

# Fill the missing SMB field over ocean with the proxy precip - evaporation
#ncf = nc.Dataset(output_file,"a")
#smb = ncf.variables["smb"][:]
# be aware: this masking is only valid for this dataset and timestep zero.
#mask_ocean = smb[0,:,:] < -0.009
#mask_ocean_expanded = np.tile(mask_ocean,(smb.shape[0],1,1))
# a proxy for smb
#smb_over_ocean = ncf.variables["precip"][:] + ncf.variables["evap"][:]
#smb[mask_ocean_expanded] = smb_over_ocean[mask_ocean_expanded]
#ncf.variables["smb"][:] = smb
#ncf.smb_comment = "SMB is approximated by precip-evap over the ocean."
#ncf.close()

subprocess.check_call('ncatted -a units,smb,o,c,"kg m-2 year-1" '+output_file,shell=True)
subprocess.check_call("ncrename -v smb,climatic_mass_balance -O "+output_file+" "+output_file,shell=True)
#subprocess.check_call('ncatted -a units,precip,o,c,"kg m-2 year-1" '+output_file,shell=True)
#subprocess.check_call("ncrename -v precip,precipitation -O "+output_file+" "+output_file,shell=True)
#subprocess.check_call("ncap2 -O -s 'air_temp=ice_surface_temp' "+output_file+" "+output_file,shell=True)

subprocess.check_call('ncks -O -C -x -v lon_2,lat_2,lon_3,lat_3 '+output_file+" "+output_file, shell=True)
subprocess.check_call("ncap2 -O -s 'x=x*1000.0;y=y*1000.0' "+output_file+" "+output_file,shell=True)

# prepare the input file for cdo remapping
# this step takes a while for high resolution data (i.e. 1km)
pi.prepare_ncfile_for_cdo(output_file)

print " RACMO file",output_file,"successfully preprocessed."
