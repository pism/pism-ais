
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
period = "hist" # "hist" or "ssp585"

data_path = os.path.join(cf.output_data_path, dataset, "input")

if not os.path.exists(data_path): os.makedirs(data_path)

output_file = os.path.join(data_path, dataset+"_"+period+"_input.nc")

if period=="hist":
    source_path = "hist_r567"
    source_file_suffix = "hist_r567_195001_201412"
elif period=="ssp585":
    source_path = "projection_2100"
    source_file_suffix = "SSP585_201501_209912"

### select required variables by (un)commenting respective lines in the dictionary below
source_file = {"alb": os.path.join(cf.racmo_cesm2_data_path,source_path,"alb_monthlyA_ANT27_CESM2_RACMO2.3p2_"+source_file_suffix+".nc"),
               "snowmelt": os.path.join(cf.racmo_cesm2_data_path,source_path,"snowmelt_monthlyS_ANT27_CESM2_RACMO2.3p2_"+source_file_suffix+".nc"),
               "tskin": os.path.join(cf.racmo_cesm2_data_path,source_path,"tskin_monthlyA_ANT27_CESM2_RACMO2.3p2_"+source_file_suffix+".nc"),
               "t2m": os.path.join(cf.racmo_cesm2_data_path,source_path,"t2m_monthlyA_ANT27_CESM2_RACMO2.3p2_"+source_file_suffix+".nc"),
               "smb": os.path.join(cf.racmo_cesm2_data_path,source_path,"smb_monthlyS_ANT27_CESM2_RACMO2.3p2_"+source_file_suffix+".nc"),
               "refreeze": os.path.join(cf.racmo_cesm2_data_path,source_path,"refreeze_monthlyS_ANT27_CESM2_RACMO2.3p2_"+source_file_suffix+".nc"),
               "runoffS": os.path.join(cf.racmo_cesm2_data_path,source_path,"runoffS_monthlyS_ANT27_CESM2_RACMO2.3p2_"+source_file_suffix+".nc"),
               "precip": os.path.join(cf.racmo_cesm2_data_path,source_path,"precip_monthlyS_ANT27_CESM2_RACMO2.3p2_"+source_file_suffix+".nc"),
               #"evap": os.path.join(cf.racmo_cesm2_data_path,source_path,"evap_monthlyS_ANT27_CESM2_RACMO2.3p2_"+source_file_suffix+".nc"),
               #"subl": os.path.join(cf.racmo_cesm2_data_path,source_path,"subl_monthlyS_ANT27_CESM2_RACMO2.3p2_"+source_file_suffix+".nc"),
               "sw0d": os.path.join(cf.racmo_cesm2_data_path,source_path,"sw0d_monthlyS_ANT27_CESM2_RACMO2.3p2_"+source_file_suffix+".nc"),
               "swsd": os.path.join(cf.racmo_cesm2_data_path,source_path,"swsd_monthlyS_ANT27_CESM2_RACMO2.3p2_"+source_file_suffix+".nc"),
               #"lwsd": os.path.join(cf.racmo_cesm2_data_path,source_path,"lwsd_monthlyS_ANT27_CESM2_RACMO2.3p2_"+source_file_suffix+".nc"), # needed to compute atmospheric emissivity
               #"clcov": os.path.join(cf.racmo_cesm2_data_path,source_path,"clcov_monthlyA_ANT27_CESM2_RACMO2.3p2_"+source_file_suffix+".nc"),
               }

var_list = source_file.keys()

process_file = {var:os.path.join(data_path, dataset+"_"+period+"_"+var+".nc") for var in var_list}

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


print "Merging files..."

merge_these_files = " ".join([process_file[var] for var in var_list])

subprocess.check_call('cdo -O merge '+merge_these_files+" "+output_file, shell=True)

# make all variables double (some already are).
ncap2str = ''.join([var+"=double("+var+");" for var in var_list]) # join commands for all variables in var_list into common string
subprocess.check_call("ncap2 -O -s '"+ncap2str+"' "+output_file+" "+output_file,shell=True)

# convert units of solar flux densities from J m-2 to W m-2 (RACMO sums over a month). Not exactly accurate, as months have different lengths.
#for var in ["sw0d","swsd","lwsd"]:
for var in ["sw0d","swsd"]:
    subprocess.check_call("ncap2 -O -s '"+var+"="+var+"/(60.*60.*24.*30.4167)' "+output_file+" "+output_file,shell=True)
    subprocess.check_call('ncatted -a units,'+var+',o,c,"W m-2 month-1" '+output_file,shell=True)

# # Option 1: convert units of SMB fluxes from kg m-2 to kg m-2 s-1 (RACMO sums over a month). Not exactly accurate, as months have different lengths.
# for var in ["snowmelt","smb","refreeze","runoff","precip"]:
#     subprocess.check_call("ncap2 -O -s '"+var+"="+var+"/(60.*60.*24.*30.4167)' "+output_file+" "+output_file,shell=True)
#     subprocess.check_call('ncatted -a units,'+var+',o,c,"kg m-2 s-1" '+output_file,shell=True)

# Option 2: convert units of SMB fluxes from kg m-2 to kg m-2 year-1 (RACMO sums over a month).
dayspermonth = np.array([31.0,28.0,31.0,30.0,31.0,30.0,31.0,31.0,30.0,31.0,30.0,31.0]) # 365_day calendar
for var in ["snowmelt","smb","refreeze","runoffS","precip"]:

    dat = nc.Dataset(output_file, 'a')
    time = dat.variables["time"][:] #days
    varcum = dat.variables[var][:] #kg/m2
    varperyear = np.zeros_like(varcum)
    for i,t in enumerate(time):
        varperyear[i]=varcum[i]*365.0/dayspermonth[i%12] #kg/m2/yr
    dat.variables[var][:]=varperyear[:]
    dat.close()

    subprocess.check_call('ncatted -a units,'+var+',o,c,"kg m-2 year-1" '+output_file,shell=True)

# # Option 3: inherit RACMO's monthly values for the SMB fluxes
# for var in ["snowmelt","smb","refreeze","runoff","precip"]:
#     subprocess.check_call('ncatted -a units,'+var+',o,c,"kg m-2 month-1" '+output_file,shell=True)

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
subprocess.check_call("ncrename -v .runoffS,runoff -O "+output_file+" "+output_file,shell=True)
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
