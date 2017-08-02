"""
matthias.mengel@pik, torsten.albrecht@pik, ronja.reese@pik
"""

import os, sys
import subprocess
## this hack is needed to import config.py from the project root
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if project_root not in sys.path: sys.path.append(project_root)
import config as cf; reload(cf)
import pism_input.pism_input as pi; reload(pi)

dataset = "racmo_wessem"
data_path = os.path.join(cf.output_data_path, dataset)

if not os.path.exists(data_path): os.makedirs(data_path)

output_file = os.path.join(data_path, dataset+"_input.nc")

source_file = {"t2m": os.path.join(cf.racmo_wessem_data_path,"t2m_RACMO2.4_yearly_ANT27_1979_2016.nc"),
                "smb": os.path.join(cf.racmo_wessem_data_path,"SMB_RACMO2.4_yearly_ANT27_1979_2016.nc")}
process_file = {var:os.path.join(data_path, dataset+"_"+var+".nc") for var in ["t2m","smb"]}

for var,fl in process_file.iteritems():
    if os.path.exists(fl): subprocess.check_call("rm -v "+fl,shell=True)

for var in ["t2m","smb"]:

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

subprocess.check_call('cdo -O merge '+process_file["smb"]+" "+process_file["t2m"]+" "+ output_file, shell=True)

# make all variables double (some already are).
subprocess.check_call("ncap2 -O -s 't2m=double(t2m);smb=double(smb)' "+
                      output_file+" "+output_file,shell=True)

# prepare the input file for cdo remapping
# this step takes a while for high resolution data (i.e. 1km)
pi.prepare_ncfile_for_cdo(process_file[var])

print " RACMO file",output_file,"successfully preprocessed."
