
"""
matthias.mengel@pik, torsten.albrecht@pik
Regrid Schmidtko data to various grid resolution using cdo remapcony.
"""

import os
import numpy as np
import sys
import netCDF4
import datetime
import subprocess


## this hack is needed to import config.py from the project root
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if project_root not in sys.path: sys.path.append(project_root)
import config as cf; reload(cf)

### Bedmap2   ##########################################################
dataset="schmidtko"
remap_file=cf.output_data_path+"tools/remap.sh"
#input_res=1 #km
input_res='ocean' #km

submit_to_cluster=False
number_of_cpu=4

datapath=os.path.join(cf.output_data_path, "ocean/"+dataset+"_data")
#datapath=os.path.join(cf.output_data_path, dataset+"/"+dataset+"_data")
submit_file="submit_remap_"+dataset+".cmd"
username = subprocess.check_output("whoami", shell=True)
username= username.rstrip("\n")



#for res in [50, 30, 20, 15, 12, 10, 7, 5, 3, 2, 1]:
for res in [15]:
  if submit_to_cluster:

    print "Submit job: Remap "+dataset+" to "+str(res)+"km grid in "+datapath

    os.system("mkdir -p log ")
    os.system("cp "+cf.output_data_path+"tools/remap_submit.cmd ./"+submit_file)

    #slurm settings
    os.system("sed -i 's/job-name=cdo_remap_/job-name=cdo_remap_"+dataset+"/' "+submit_file)
    os.system("sed -i 's/cpus-per-task=8/cpus-per-task="+str(number_of_cpu)+"/' "+submit_file)
    os.system("sed -i 's/mail-user=albrecht/mail-user="+username+"/' "+submit_file)

    cmd="sbatch "+submit_file+" "+dataset+" "+str(res)+" "+str(input_res) #+" "+str(number_of_cpu)
    print cmd
    os.system(cmd)

  else: #interactive

    cmd="bash "+remap_file+" "+dataset+" "+str(res)+" "+str(input_res)
    print cmd
    os.system(cmd)
