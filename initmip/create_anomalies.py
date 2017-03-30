#!/usr/bin/env python
# Copyright (C) 2015 Andy Aschwanden, 2017 Torsten Albrecht

# see https://github.com/pism/initMIP/blob/master/prepare.sh

import os
import numpy as np
from argparse import ArgumentParser
from netCDF4 import Dataset as NC


# Set up the option parser
parser = ArgumentParser()
parser.description = "Create initMIP SMB and BMB anomalies for Antarctica."
parser.add_argument("--force_file", dest="force_file",
                    help='''Anomaly mass balance file''')
parser.add_argument("--background_file", dest="background_file",
                    help='''Background SMB file''')
parser.add_argument("OUTFILE", nargs=1)
options = parser.parse_args()
force_file = options.force_file
background_file = options.background_file

outfile = options.OUTFILE[0]


rho_ice=910.0 #kg/m^3
secperyear=31556926.0 #s/yr.


def get_dims(nc): #from https://github.com/pism/pypismtools/blob/master/pypismtools.py
    '''
    Gets dimensions from netcdf instance
    Parameters:
    -----------
    nc: netCDF instance
    Returns:
    --------
    xdim, ydim, zdim, tdim: dimensions
    '''

    # a list of possible x-dimensions names
    xdims = ['x', 'x1']
    # a list of possible y-dimensions names
    ydims = ['y', 'y1']
    # a list of possible z-dimensions names
    zdims = ['z', 'z1']
    # a list of possible time-dimensions names
    tdims = ['t', 'time']

    xdim = None
    ydim = None
    zdim = None
    tdim = None

    # assign x dimension
    for dim in xdims:
        if dim in list(nc.dimensions.keys()):
            xdim = dim
    # assign y dimension
    for dim in ydims:
        if dim in list(nc.dimensions.keys()):
            ydim = dim
    # assign z dimension
    for dim in zdims:
        if dim in list(nc.dimensions.keys()):
            zdim = dim
    # assign time dimension
    for dim in tdims:
        if dim in list(nc.dimensions.keys()):
            tdim = dim
    return xdim, ydim, zdim, tdim

###########################################################


nc_a = NC(force_file, 'r')
nc_b = NC(background_file, 'r')

try:
    os.remove(outfile)
except OSError:
    pass
nc = NC(outfile, 'w')

xdim_a, ydim_a, zdim_a, tdim_a = get_dims(nc_a)
xdim_b, ydim_b, zdim_b, tdim_b = get_dims(nc_b)

assert xdim_a == xdim_b
assert ydim_a == ydim_b

xdim = xdim_a
ydim = ydim_a
tdim = 'time'

nx = len(nc_a.dimensions[xdim_a])
ny = len(nc_a.dimensions[ydim_b])

start = -5
end = 100
nt = end - start

nc.createDimension(xdim, size = (nx))
nc.createDimension(ydim, size = (ny))
nc.createDimension(tdim)

bnds_var_name = "time_bnds"
# create a new dimension for bounds only if it does not yet exist
bnds_dim = "nb2"
if bnds_dim not in nc.dimensions.keys():
    nc.createDimension(bnds_dim, 2)

time_var = nc.createVariable(tdim, 'float64', dimensions=(tdim))
time_var.bounds = bnds_var_name
time_var.units = 'years'
time_var.axis = 'T'
time_var[:] = range(start, end)

# create time bounds variable
time_bnds_var = nc.createVariable(bnds_var_name, 'd', dimensions=(tdim, bnds_dim))
time_bnds_var[:, 0] = range(start, end)
time_bnds_var[:, 1] = range(start+1, end+1)

smb_background = nc_b.variables['effective_climatic_mass_balance']
temp_background = nc_b.variables['effective_ice_surface_temp']
bmb_background = nc_b.variables['effective_shelf_base_mass_flux']
btemp_background = nc_b.variables['effective_shelf_base_temperature']

#smb_background = nc_b.variables['climatic_mass_balance']
#temp_background = nc_b.variables['ice_surface_temp']
#bmb_background = nc_b.variables['shelfbmassflux']
#btemp_background = nc_b.variables['shelfbtemp']

smb_anomaly = nc_a.variables['asmb']
bmb_anomaly = nc_a.variables['abmb']

smb_var = nc.createVariable('climatic_mass_balance', 'float64', dimensions=(tdim, ydim, xdim))
temp_var = nc.createVariable('ice_surface_temp', 'float64', dimensions=(tdim, ydim, xdim))
bmb_var = nc.createVariable('shelfbmassflux', 'float64', dimensions=(tdim, ydim, xdim))
btemp_var = nc.createVariable('shelfbtemp', 'float64', dimensions=(tdim, ydim, xdim))

temp_var.units = "Kelvin"
temp_var.long_name = "temperature of the ice at the ice surface but below firn processes, as seen by the ice dynamics code"
smb_var.units = "kg m-2 year-1"
smb_var.long_name = "surface mass balance (accumulation/ablation) rate, as seen by the ice dynamics code"
btemp_var.units = "Kelvin"
btemp_var.long_name = "shelf base temperature, as seen by the ice dynamics code"
bmb_var.units = "kg m-2 year-1"
bmb_var.long_name = "shelf base mass flux (positive flux is loss from ice shelf), as seen by the ice dynamics code"

for k in range(nt):
    t = k + start
    temp_var[k,::] = np.squeeze(temp_background[:])
    btemp_var[k,::] = np.squeeze(btemp_background[:])
    if t < 0:
        smb_var[k,::] = np.squeeze(smb_background[:])*secperyear
        bmb_var[k,::] = np.squeeze(bmb_background[:])
    elif (t >= 0) and (t < 40):
        smb_var[k,::] = np.squeeze(smb_background[:])*secperyear + np.squeeze(smb_anomaly[:])*rho_ice * np.floor(t) / 40
        bmb_var[k,::] = np.squeeze(bmb_background[:]) + np.squeeze(bmb_anomaly[:])*rho_ice * np.floor(t) / 40
    else:
        smb_var[k,::] = np.squeeze(smb_background[:])*secperyear + np.squeeze(smb_anomaly[:])*rho_ice
        bmb_var[k,::] = np.squeeze(bmb_background[:]) + np.squeeze(bmb_anomaly[:])*rho_ice
        

nc.close()
nc_a.close()
nc_b.close()
