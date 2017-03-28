"""
matthias.mengel@pik, torsten.albrecht@pik, ronja.reese@pik

INFO to the file:
The file is based on HadCM3_c20_I2S_precip_Y.nc and HadCM3_c20_I2S_t2m_Y.nc which come from data for Frieler et al 2014.
HadCM3_c20_I2S_precip_Y.nc contains precipitation in mm/yr for 1980-1999 over Antarctica from RACMO2 run forced with HadCM3 data. (ATTENTION units are wrong in this file!)
HadCM3_c20_I2S_t2m_Y.nc contains temperature in Kelvin for 1980-1999 over Antarctica from RACMO2 run forced with HadCM3 data.
The variables are extraced from the files, precip is changed in mm/yr ice equivalent and the temporal mean is taken. The variables are renamed and interpolated to a PISM grid.
These actions are documented in the file /home/reese/data/prepareData_forPISM/prepare.for.pism_HadCM3_c20.sh
There are some missing_values which are interpolated using the file fill_missing.py from PISM/util.
This file should be used as PISM input with the atmosphere_given option.
"""

import os, sys
import subprocess
## this hack is needed to import config.py from the project root
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if project_root not in sys.path: sys.path.append(project_root)
import config as cf; reload(cf)

dataset = "racmo_hadcm3_I2S"
data_path = os.path.join(cf.output_data_path, dataset)

if not os.path.exists(data_path): os.makedirs(data_path)

scenarios = ["c20","A1B"]
accumvar = {"c20":"precip","A1B":"smb"}

for scen in scenarios:

    output_file = os.path.join(data_path,"racmo_hadcm3_"+scen+"_input.nc")

    source_t2m_file = os.path.join(cf.racmo_data_path,"HadCM3_"+scen+"_I2S_t2m_Y.nc")
    source_smb_file = os.path.join(cf.racmo_data_path,"HadCM3_"+scen+"_I2S_"+accumvar[scen]+"_Y.nc")
    subprocess.check_call("cp "+source_t2m_file+" "+output_file,shell=True)

    #add temperature
    subprocess.check_call('ncks -A -v t2m '+source_t2m_file+" "+output_file,shell=True)
    subprocess.check_call('ncrename -v t2m,air_temp -O '+output_file+" "+output_file,shell=True)

    #add precip/smb
    subprocess.check_call('ncks -A -v '+accumvar[scen]+",lat,lon,time "+source_smb_file
                          +" "+output_file,shell=True)
    subprocess.check_call('ncrename -v '+accumvar[scen]+",precipitation -O "+output_file+" "+output_file,shell=True)

    subprocess.check_call('ncrename -d rlon,x -d rlat,y -O '+output_file+" "+output_file,shell=True)
    subprocess.check_call('ncrename -v rlon,x -v rlat,y -O '+output_file+" "+output_file,shell=True)
    subprocess.check_call('ncap2 -O -s "x=x*100000.0;y=y*100000.0" '+output_file+" "+output_file,shell=True)
    subprocess.check_call('ncatted -a units,x,o,c,"meters" '+output_file,shell=True)

    # units were water equivalent per year, make it ice equivalent,
    # rho_ice=910 and rho_freshwater=1000 taken from pism config
    subprocess.check_call("ncap2 -O -s 'precipitation=precipitation*1000.0/910.0' "+output_file+" "+output_file,shell=True)
    #catted -a units,climatic_mass_balance,o,c,"mm a-1" +output_file #kg/m2/yr
    subprocess.check_call('ncatted -a units,precipitation,o,c,"kg m-2 year-1" '+output_file,shell=True) #kg/m2/yr
    subprocess.check_call('ncatted -a coordinates,precipitation,o,c,"lon lat" '+output_file,shell=True)
    subprocess.check_call('ncatted -a units,air_temp,o,c,"Kelvin" '+output_file,shell=True)
    subprocess.check_call('ncatted -a coordinates,air_temp,o,c,"lon lat" '+output_file,shell=True)
    subprocess.check_call('ncatted -a units,time,o,c,"years since 01-01-0000" '+output_file,shell=True)
    subprocess.check_call('ncatted -a calendar,time,o,c,"360_day" '+output_file,shell=True)
    subprocess.check_call('ncatted -a dtgstart,time,d,c,1980010100 '+output_file,shell=True)
    subprocess.check_call('ncatted -a axis,time,d,c,T '+output_file,shell=True)
    subprocess.check_call('ncatted -a long_name,time,d,c,"Julian Day" '+output_file,shell=True)

    #time average
    #ncra -F -O -d time,1,20,1 +output_file+" "+output_file

    subprocess.check_call('ncwa -O -a height '+output_file+" "+output_file, shell=True) #delete the height dimension!
    subprocess.check_call('ncks -O -x -v height '+output_file+" "+output_file, shell=True) #delete the height dimension!

    subprocess.check_call('ncatted -O -a projection,global,a,c,"+lon_0=0.0 +ellps=WGS84 +datum=WGS84 +lat_ts=-71.0 +proj=stere +x_0=0.0 +units=m +y_0=0.0 +lat_0=-90.0" '+output_file,shell=True)
    subprocess.check_call('ncatted -a units,time_bnds,o,c,"years since 01-01-0000" '+output_file,shell=True)

    print " RACMO file",output_file,"successfully preprocessed."
