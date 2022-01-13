#!/usr/bin/env python3
# coding: utf-8

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

import cftime
import os,sys
import datetime

## this hack is needed to import config.py from the project root
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if project_root not in sys.path: sys.path.append(project_root)
import config as cf#; reload(cf)


# settings
year_start = cf.time_averaging_period[0]
year_end   = cf.time_averaging_period[1]
#year_end   = 1985
year_clim  = 1
grid = cf.grid_id
dataset = "racmo_era5"

path = f'/p/projects/pism/kreuzer/pism_input_2020/{dataset}/'
file_input = os.path.join(path, f'{dataset}_{grid}.nc')
file_output = os.path.join(path,f'{dataset}_{grid}_{year_start}_{year_end}_clim_month.nc')

filename=__file__



ds = xr.open_dataset(file_input, use_cftime=True)
ds_tsel = ds.sel(time=slice(f"{year_start}-01-01",f"{year_end}-12-31"))

ds_clim = ds_tsel.groupby("time.month").mean('time', keep_attrs=True)
ds_clim2 = ds_clim.copy(deep=True)

# set new time axis
time_noleap_1yr_month = [cftime.DatetimeNoLeap(year_clim, month, 15) for month in range(1, 13)]
ds_clim2['month'] = time_noleap_1yr_month
ds_clim2 = ds_clim2.rename_dims(dict(month='time'))
ds_clim2 = ds_clim2.rename_vars(dict(month='time'))

ds_clim2['time'].attrs['standard_name'] = 'time'
ds_clim2['time'].attrs['long_name'] = 'time'
ds_clim2['time'].attrs['bounds'] = 'time_bnds'

# set time_bnds
time_bnds = np.zeros((12,2), cftime._cftime.DatetimeNoLeap)
for month in range(1,13):
    time_bnds[month-1][0] = cftime.DatetimeNoLeap(year_clim, month, 1)
    if month // 12 == 0:
        time_bnds[month-1][1] = cftime.DatetimeNoLeap(year_clim, month+1, 1)
    else:
        time_bnds[month-1][1] = cftime.DatetimeNoLeap(year_clim+1, (month+1)%12, 1)
ds_clim2['time_bnds'] = (['time','bnds'], time_bnds)

# add to history attribute
timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
hist_prepend = f'{timestamp}: monthly climatology (year {year_start}-{year_end}) computed by {filename}\n'  
ds_clim2.attrs['history'] = hist_prepend + ds_clim2.attrs['history']

# write to file
ds_clim2.to_netcdf(file_output, 
                   unlimited_dims = ['time'], 
                   encoding={"time":      {"dtype": "float", "units":f"days since {year_clim:04}-01-01 00:00:00"},
                             "time_bnds": {"dtype": "float", "units":f"days since {year_clim:04}-01-01 00:00:00"}})


