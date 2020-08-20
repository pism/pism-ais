"""
matthias.mengel@pik, julius.garbe@pik

A template for time_averaging data.
Applied after remapping. May also be done before it.

"""

import os, sys
import subprocess
import netCDF4 as nc
## this hack is needed to import config.py from the project root
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if project_root not in sys.path: sys.path.append(project_root)
import config as cf; reload(cf)


def run(commands):
    """Run a list of commands (need to be given as strings)."""
    for cmd in commands:
        print cmd
        subprocess.check_call(cmd,shell=True)


dataset = "racmo_cesm2"
period = "hist"
#period = "ssp585"
period_years = "1950_2014" # hist
#period_years = "2015_2099" # ssp585
timestep = "monthly" # time step of RACMO input data ("monthly" or "yearly")
grid = cf.grid_id
tap = [str(t) for t in cf.time_averaging_period]

data_path = os.path.join(cf.output_data_path, dataset)
#inputfile = os.path.join(data_path, dataset+"_"+period+"_"+grid+".nc")
inputfile = os.path.join(data_path, grid, dataset+"_"+period+"_"+period_years+"_"+timestep+"_"+grid+".nc")
#timemean_file = os.path.join(data_path, dataset+"_"+period+"_"+grid+"_mean"+tap[0]+"_"+tap[1]+".nc")
timemean_file = os.path.join(data_path, grid, dataset+"_"+period+"_"+tap[0]+"_"+tap[1]+"_mean_"+grid+".nc")

if timestep=="yearly":
    # If data is already yearly, we just apply time mean over the selected period
    commands = [
                "cdo -O -f nc4c timmean -selyear,{}/{} {} {}".format(tap[0], tap[1], inputfile, timemean_file),
                ]
elif timestep=="monthly":
    # If data is monthly, we have to be cautious:
    # monthly temperatures need to be AVERAGED over the year,
    # while monthly accumulation rates need to be SUMMED over the year.
    # This is done here by creating separate files for both variables and merging them afterwards.
    temp_var = "air_temp"
    accu_var = "precipitation"
    commands = [
                "cdo -O -f nc4c yearmean -select,name={} {} {}".format(temp_var, inputfile, inputfile.replace('monthly','yearly_'+temp_var)),
                "cdo -O -f nc4c yearsum -select,name={} {} {}".format(accu_var, inputfile, inputfile.replace('monthly','yearly_'+accu_var)),
                "cdo -O merge {} {} {}".format(inputfile.replace('monthly','yearly_'+temp_var), inputfile.replace('monthly','yearly_'+accu_var), inputfile.replace('monthly','yearly')),
                "cdo -O -f nc4c timmean -selyear,{}/{} {} {}".format(tap[0], tap[1], inputfile.replace('monthly','yearly'), timemean_file),
                ]

# if you see "cdo: error while loading shared libraries:",
# 'module load cdo'
run(commands)
