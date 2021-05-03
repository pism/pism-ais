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
period = "hist" # historical: "hist" / SSP5-8.5: "ssp585"
period_years = "1950_2014" # hist: "1950_2014" / ssp585: "2015_2099"
timestep = "monthly" # time step of input data ("monthly" or "yearly")
climatology = True # only if timestep is monthly: set to True if you want to have a yearly climatology (e.g. when using PISM with periodic forcing)


grid = cf.grid_id
tap = [str(t) for t in cf.time_averaging_period]

data_path = os.path.join(cf.output_data_path, dataset)
inputfile = os.path.join(data_path, grid, dataset+"_"+period+"_"+period_years+"_"+timestep+"_"+grid+".nc")
timemean_file = os.path.join(data_path, grid, dataset+"_"+period+"_"+tap[0]+"_"+tap[1]+"_mean_"+grid+".nc")


if timestep=="yearly":
    # If data is already yearly, we just apply time mean over the selected period
    commands = [
                "cdo -O -f nc4c timmean -selyear,{}/{} {} {}".format(tap[0], tap[1], inputfile, timemean_file),
                ]
elif timestep=="monthly":
    if not climatology:
        # If data is monthly, we have to be cautious when calculating yearly values:
        # Some variables (albedo, temperatures, ...) need to be AVERAGED over the year, while other variables (SMB rates, insolation rates, ...) need to be SUMMED over the year.
        # This is done here by creating separate files for both variable types and merging them afterwards.
        mean_vars = "albedo,air_temp,ice_surface_temp"
        sum_vars = "climatic_mass_balance,evaporation,incoming_shortwave_flux_surf,incoming_shortwave_flux_toa,precipitation,refreeze,runoff,snowmelt,sublimation"
        commands = [
                    "cdo -O -f nc4c yearmean -selname,{} {} {}".format(mean_vars, inputfile, inputfile.replace('monthly','yearly_MEAN_')),
                    "cdo -O -f nc4c yearsum -selname,{} {} {}".format(sum_vars, inputfile, inputfile.replace('monthly','yearly_SUM_')),
                    "cdo -O merge {} {} {}".format(inputfile.replace('monthly','yearly_MEAN_'), inputfile.replace('monthly','yearly_SUM_'), inputfile.replace('monthly','yearly')),
                    "cdo -O -f nc4c timmean -selyear,{}/{} {} {}".format(tap[0], tap[1], inputfile.replace('monthly','yearly'), timemean_file),
                    ]
    elif climatology:
        commands = [
                    "cdo -O -f nc4c ymonmean -selyear,{}/{} {} {}".format(tap[0], tap[1], inputfile, timemean_file.replace('mean','ymonmean')),
                   ]

# if you see "cdo: error while loading shared libraries:",
# 'module load cdo'
run(commands)
