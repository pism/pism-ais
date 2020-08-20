
"""
matthias.mengel@pik, torsten.albrecht@pik, ronja.reese@pik, julius.garbe@pik
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

dataset = "racmo_cesm2"
data_path = os.path.join(cf.output_data_path, dataset)

if not os.path.exists(data_path): os.makedirs(data_path)

output_file = os.path.join(data_path, dataset+"_input_hist.nc")

source_file = {"tskin": os.path.join(cf.racmo_cesm2_data_path,"hist/tskin_monthlyA_ANT27_CESM2_RACMO2.3p2_195001_201412.nc"),
               "smb": os.path.join(cf.racmo_cesm2_data_path,"hist/smb_monthlyS_ANT27_CESM2_RACMO2.3p2_195001_201412.nc")}
               #"t2m": os.path.join(cf.racmo_cesm2_data_path,"hist/t2m_monthlyA_ANT27_CESM2_RACMO2.3p2_195001_201412.nc"),
               #"precip": os.path.join(cf.racmo_cesm2_data_path,"hist/precip_monthlyS_ANT27_CESM2_RACMO2.3p2_195001_201412.nc"),
               #}

process_file = {var:os.path.join(data_path, dataset+"_"+var+"_hist.nc") for var in ["tskin","smb"]}

for var,fl in process_file.iteritems():
    try:
        os.remove(fl)
    except OSError:
        pass

#for var in ["tskin","smb","t2m","precip"]:
for var in ["tskin","smb"]:

    subprocess.check_call("ncks -A -v "+var+",lon,lat "+source_file[var]+" "+process_file[var],shell=True)

    subprocess.check_call("ncrename -d rlat,y -d rlon,x -v rlat,y -v rlon,x "+process_file[var],shell=True)

    subprocess.check_call("ncatted -O -a grid_mapping,"+var+",d,, "+process_file[var],shell=True)

    # Here any string could be used, if the ob_trans projection is not know by the proj version. For bilinear remapping only the lon lat values count.
    subprocess.check_call('ncatted -O -a proj4,global,o,c,"+lon_0=10.0 +ellps=WGS84 +datum=WGS84 +lat_ts=-71.0 +proj=stere +x_0=0.0 +units=m +y_0=0.0 +lat_0=-90.0" '+process_file[var], shell=True)
    #subprocess.check_call('ncatted -O -a proj4,global,o,c,"-m 57.295779506 +proj=ob_tran +o_proj=latlon +o_lat_p=-180.0 +lon_0=10.0 +x_0=0.0 +units=m +y_0=0.0" '+process_file[var], shell=True)


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

#merge_these_files = " ".join([process_file[var] for var in ["tskin","smb","t2m","precip"]])
merge_these_files = " ".join([process_file[var] for var in ["tskin","smb"]])

subprocess.check_call('cdo -O merge '+merge_these_files+" "+output_file, shell=True)

# make all variables double (some already are).
subprocess.check_call("ncap2 -O -s 'tskin=double(tskin);smb=double(smb);' "+
                      output_file+" "+output_file,shell=True)

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

# attribute fixes
subprocess.check_call('ncatted -a units,x,o,c,"meters" '+output_file,shell=True)
subprocess.check_call('ncatted -a units,y,o,c,"meters" '+output_file,shell=True)
#subprocess.check_call('ncatted -a units,time_bnds,o,c,"days since 1950-01-01 00:00:00.0" '+output_file,shell=True)

subprocess.check_call('ncatted -a units,smb,o,c,"kg m-2 year-1" '+output_file,shell=True)
#subprocess.check_call('ncatted -a units,precip,o,c,"kg m-2 year-1" '+output_file,shell=True)

# rename variables
subprocess.check_call("ncrename -v tskin,ice_surface_temp -O "+output_file+" "+output_file,shell=True)
subprocess.check_call("ncrename -v smb,climatic_mass_balance -O "+output_file+" "+output_file,shell=True)
#subprocess.check_call("ncrename -v t2m,air_temp -O "+output_file+" "+output_file,shell=True)
#subprocess.check_call("ncrename -v precip,precipitation -O "+output_file+" "+output_file,shell=True)


subprocess.check_call('ncks -O -C -x -v lon_2,lat_2,lon_3,lat_3 '+output_file+" "+output_file, shell=True)
#RACMO grid actually comes on a lon lat coordinate, so multipliaction here provides values close to meters, but this is not important here
subprocess.check_call("ncap2 -O -s 'x=x*1.0e5;y=y*1.0e5' "+output_file+" "+output_file,shell=True)
#subprocess.check_call("ncks -C -O -x -v lon,lat "+output_file+" "+output_file,shell=True)


# prepare the input file for cdo remapping
# this step takes a while for high resolution data (i.e. 1km)
pi.prepare_ncfile_for_cdo(output_file)

print " RACMO file",output_file,"successfully preprocessed."
