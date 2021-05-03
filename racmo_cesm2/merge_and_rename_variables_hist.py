
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
data_path = os.path.join(cf.output_data_path, dataset, "input")

if not os.path.exists(data_path): os.makedirs(data_path)

output_file = os.path.join(data_path, dataset+"_hist_input.nc")

source_file = {"alb": os.path.join(cf.racmo_cesm2_data_path,"hist_r567/alb_monthlyA_ANT27_CESM2_RACMO2.3p2_hist_r567_195001_201412.nc"),
               "snowmelt": os.path.join(cf.racmo_cesm2_data_path,"hist_r567/snowmelt_monthlyS_ANT27_CESM2_RACMO2.3p2_hist_r567_195001_201412.nc"),
               #"tskin": os.path.join(cf.racmo_cesm2_data_path,"hist_r567/tskin_monthlyA_ANT27_CESM2_RACMO2.3p2_hist_r567_195001_201412.nc"),
               "t2m": os.path.join(cf.racmo_cesm2_data_path,"hist_r567/t2m_monthlyA_ANT27_CESM2_RACMO2.3p2_hist_r567_195001_201412.nc"),
               "smb": os.path.join(cf.racmo_cesm2_data_path,"hist_r567/smb_monthlyS_ANT27_CESM2_RACMO2.3p2_hist_r567_195001_201412.nc"),
               "refreeze": os.path.join(cf.racmo_cesm2_data_path,"hist_r567/refreeze_monthlyS_ANT27_CESM2_RACMO2.3p2_hist_r567_195001_201412.nc"),
               "runoff": os.path.join(cf.racmo_cesm2_data_path,"hist_r567/runoff_monthlyS_ANT27_CESM2_RACMO2.3p2_hist_r567_195001_201412.nc"),
               "precip": os.path.join(cf.racmo_cesm2_data_path,"hist_r567/precip_monthlyS_ANT27_CESM2_RACMO2.3p2_hist_r567_195001_201412.nc"),
               #"evap": os.path.join(cf.racmo_cesm2_data_path,"hist_r567/evap_monthlyS_ANT27_CESM2_RACMO2.3p2_hist_r567_r567_195001_201412.nc"),
               #"subl": os.path.join(cf.racmo_cesm2_data_path,"hist_r567/subl_monthlyS_ANT27_CESM2_RACMO2.3p2_hist_r567_195001_201412.nc"),
               "sw0d": os.path.join(cf.racmo_cesm2_data_path,"hist_r567/sw0d_monthlyS_ANT27_CESM2_RACMO2.3p2_hist_r567_195001_201412.nc"),
               "swsd": os.path.join(cf.racmo_cesm2_data_path,"hist_r567/swsd_monthlyS_ANT27_CESM2_RACMO2.3p2_hist_r567_195001_201412.nc"),
               #"lwsd": os.path.join(cf.racmo_cesm2_data_path,"hist_r567/lwsd_monthlyS_ANT27_CESM2_RACMO2.3p2_hist_r567_195001_201412.nc"), # needed to compute atmospheric emissivity
               #"clcov": os.path.join(cf.racmo_cesm2_data_path,"hist_r567/clcov_monthlyA_ANT27_CESM2_RACMO2.3p2_hist_r567_195001_201412.nc"),
               }

var_list = ["alb","snowmelt","t2m","smb","refreeze","runoff","precip","sw0d","swsd"]

process_file = {var:os.path.join(data_path, dataset+"_hist_"+var+".nc") for var in var_list}

for var,fl in process_file.iteritems():
    try:
        os.remove(fl)
    except OSError:
        pass

for var in var_list:
    print var

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

    #if var=="smb":
    #    #356day calendar
    #    months=np.array([31.0,28.0,31.0,30.0,31.0,30.0,31.0,31.0,30.0,31.0,30.0,31.0])

    #    dat = nc.Dataset(process_file[var], 'a')
    #    time = dat.variables["time"][:] #days
    #    smbcum = dat.variables["smb"][:] #kg/m2
    #    smbperyear = np.zeros_like(smbcum)
    #    for i,t in enumerate(time):
    #        smbperyear[i]=smbcum[i]*365.0/months[i%12] #kg/m2/yr
    #    dat.variables["smb"][:]=smbperyear[:]  
    #    dat.close()

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

print "Merging files..."

merge_these_files = " ".join([process_file[var] for var in var_list])

subprocess.check_call('cdo -O merge '+merge_these_files+" "+output_file, shell=True)

# make all variables double (some already are).
ncap2str = ''.join([var+"=double("+var+");" for var in var_list]) # join commands for all variables in var_list into common string
subprocess.check_call("ncap2 -O -s '"+ncap2str+"' "+output_file+" "+output_file,shell=True)

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

# convert units of solar flux densities from J m-2 to W m-2 (RACMO sums over a month). Not exactly accurate, as months have different lengths.
subprocess.check_call("ncap2 -O -s 'sw0d=sw0d/(60.*60.*24.*30.4167)' "+output_file+" "+output_file,shell=True)
subprocess.check_call("ncap2 -O -s 'swsd=swsd/(60.*60.*24.*30.4167)' "+output_file+" "+output_file,shell=True)
#subprocess.check_call("ncap2 -O -s 'lwsd=lwsd/(60.*60.*24.*30.4167)' "+output_file+" "+output_file,shell=True)

subprocess.check_call('ncatted -a units,sw0d,o,c,"W m-2" '+output_file,shell=True)
subprocess.check_call('ncatted -a units,swsd,o,c,"W m-2" '+output_file,shell=True)
#subprocess.check_call('ncatted -a units,lwsd,o,c,"W m-2" '+output_file,shell=True)

# convert units of SMB fluxes from kg m-2 to kg m-2 s-1 (RACMO sums over a month). Not exactly accurate, as months have different lengths.
subprocess.check_call("ncap2 -O -s 'snowmelt=snowmelt/(60.*60.*24.*30.4167)' "+output_file+" "+output_file,shell=True)
subprocess.check_call("ncap2 -O -s 'smb=smb/(60.*60.*24.*30.4167)' "+output_file+" "+output_file,shell=True)
subprocess.check_call("ncap2 -O -s 'refreeze=refreeze/(60.*60.*24.*30.4167)' "+output_file+" "+output_file,shell=True)
subprocess.check_call("ncap2 -O -s 'runoff=runoff/(60.*60.*24.*30.4167)' "+output_file+" "+output_file,shell=True)
subprocess.check_call("ncap2 -O -s 'precip=precip/(60.*60.*24.*30.4167)' "+output_file+" "+output_file,shell=True)
#subprocess.check_call("ncap2 -O -s 'evap=evap/(60.*60.*24.*30.4167)' "+output_file+" "+output_file,shell=True)
#subprocess.check_call("ncap2 -O -s 'subl=subl/(60.*60.*24.*30.4167)' "+output_file+" "+output_file,shell=True)

subprocess.check_call('ncatted -a units,snowmelt,o,c,"kg m-2 s-1" '+output_file,shell=True)
subprocess.check_call('ncatted -a units,smb,o,c,"kg m-2 s-1" '+output_file,shell=True)
subprocess.check_call('ncatted -a units,refreeze,o,c,"kg m-2 s-1" '+output_file,shell=True)
subprocess.check_call('ncatted -a units,runoff,o,c,"kg m-2 s-1" '+output_file,shell=True)
subprocess.check_call('ncatted -a units,precip,o,c,"kg m-2 s-1" '+output_file,shell=True)
#subprocess.check_call('ncatted -a units,evap,o,c,"kg m-2 s-1" '+output_file,shell=True)
#subprocess.check_call('ncatted -a units,subl,o,c,"kg m-2 s-1" '+output_file,shell=True)

# further attribute fixes
subprocess.check_call('ncatted -a units,x,o,c,"meters" '+output_file,shell=True)
subprocess.check_call('ncatted -a units,y,o,c,"meters" '+output_file,shell=True)

subprocess.check_call('ncatted -a units,time_bnds,o,c,"days since 1950-01-01 00:00:00.0" '+output_file,shell=True)
subprocess.check_call('ncatted -a coordinates,time_bnds,d,, '+output_file,shell=True)

subprocess.check_call('ncatted -a units,alb,o,c,"-" '+output_file,shell=True)

# rename variables
subprocess.check_call("ncrename -v .alb,albedo -O "+output_file+" "+output_file,shell=True)
subprocess.check_call("ncrename -v .tskin,ice_surface_temp -O "+output_file+" "+output_file,shell=True)
subprocess.check_call("ncrename -v .t2m,air_temp -O "+output_file+" "+output_file,shell=True)
subprocess.check_call("ncrename -v .smb,climatic_mass_balance -O "+output_file+" "+output_file,shell=True)
subprocess.check_call("ncrename -v .precip,precipitation -O "+output_file+" "+output_file,shell=True)
subprocess.check_call("ncrename -v .evap,evaporation -O "+output_file+" "+output_file,shell=True)
subprocess.check_call("ncrename -v .subl,sublimation -O "+output_file+" "+output_file,shell=True)
subprocess.check_call("ncrename -v .sw0d,incoming_shortwave_radiation_TOA -O "+output_file+" "+output_file,shell=True)
subprocess.check_call("ncrename -v .swsd,incoming_shortwave_radiation -O "+output_file+" "+output_file,shell=True)
subprocess.check_call("ncrename -v .lwsd,incoming_longwave_radiation -O "+output_file+" "+output_file,shell=True)
subprocess.check_call("ncrename -v .clcov,cloud_cover -O "+output_file+" "+output_file,shell=True)

subprocess.check_call('ncks -O -C -x -v lon_2,lat_2,lon_3,lat_3 '+output_file+" "+output_file, shell=True)
#RACMO grid actually comes on a lon lat coordinate, so multipliaction here provides values close to meters, but this is not important here
subprocess.check_call("ncap2 -O -s 'x=x*1.0e5;y=y*1.0e5' "+output_file+" "+output_file,shell=True)
#subprocess.check_call("ncks -C -O -x -v lon,lat "+output_file+" "+output_file,shell=True)

## since RACMO doesn't include atmospheric emissivity, we compute it according to the Stefan-Boltzmann law from longwave downward surface radiation and 2m air temperature
#subprocess.check_call("ncap2 -O -s 'emissivity=incoming_longwave_radiation/(5.67e-8 * air_temp^4)' "+output_file+" "+output_file,shell=True)
#subprocess.check_call('ncatted -a long_name,emissivity,o,c,"Atmospheric Emissivity" '+output_file,shell=True)
#subprocess.check_call('ncatted -a standard_name,emissivity,o,c,"atmospheric_emissivity" '+output_file,shell=True)
#subprocess.check_call('ncatted -a units,emissivity,o,c,"-" '+output_file,shell=True)
#subprocess.check_call("ncks -O -x -v incoming_longwave_radiation "+output_file+" "+output_file,shell=True) # delete, since it is not needed anymore

# prepare the input file for cdo remapping
# this step takes a while for high resolution data (i.e. 1km)
pi.prepare_ncfile_for_cdo(output_file)

print " RACMO file",output_file,"successfully preprocessed."
