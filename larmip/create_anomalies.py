#!/usr/bin/env python
# Copyright (C) 2015 Andy Aschwanden, 2017-18 Torsten Albrecht

# see https://github.com/pism/initMIP/blob/master/prepare.sh

import os
import numpy as np
from argparse import ArgumentParser
from netCDF4 import Dataset as NC


# Set up the option parser
parser = ArgumentParser()
parser.description = "Create LARMIP BMB anomalies for Antarctica."
parser.add_argument("--region_file", dest="region_file",
                    help='''LARMIP regions file''')
parser.add_argument("--background_file", dest="background_file",
                    help='''Background mass balance file''')
parser.add_argument("--region", dest="region",
                    help='''LARMIP region index''')
parser.add_argument("--meltrate", dest="meltrate",
                    help='''Melt rate anomaly''')
parser.add_argument("OUTFILE", nargs=1)
options = parser.parse_args()
region_file = options.region_file
background_file = options.background_file

larmip_region = np.float(options.region)
mrate_anomaly = np.float(options.meltrate)

outfile = options.OUTFILE[0]


rho_ice=910.0 #kg/m^3
secperyear=31556926.0 #s/yr. FIXME: which calendar?
secperyear=365.0*24.0*3600.0

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


nc_a = NC(region_file, 'r')
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

dt = 5
start = -5
end = 200
nt = (end - start)/dt

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
time_steps = np.arange(start, end, dt)
time_var[:] = time_steps

# create time bounds variable
time_bnds_var = nc.createVariable(bnds_var_name, 'd', dimensions=(tdim, bnds_dim))
time_bnds_var[:, 0] = time_steps
time_bnds_var[:, 1] = time_steps+dt

smb_background = nc_b.variables['climatic_mass_balance']
temp_background = nc_b.variables['ice_surface_temp']
bmb_background = nc_b.variables['shelfbmassflux']
btemp_background = nc_b.variables['shelfbtemp']


larmip_regs = nc_a.variables['regions']
larmip_regid = np.squeeze(larmip_regs[:]) 

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

larmip_regs = np.squeeze(nc_a.variables['regions'][:])

for k,t in enumerate(time_steps):
#for k in range(nt):
    #t = k + start
    temp_var[k,::] = np.squeeze(temp_background[:])
    btemp_var[k,::] = np.squeeze(btemp_background[:])
    smb_var[k,::] = np.squeeze(smb_background[:])
    bmb_var[k,::] = np.squeeze(bmb_background[:])
    if (t >= 0):
        if larmip_region == 5:
          bmelt_anomaly = np.ones_like(bmb_var[k,::])*mrate_anomaly*rho_ice
        else:
          bmelt_anomaly = np.zeros_like(bmb_var[k,::])
          bmelt_anomaly[ (larmip_regs==larmip_region) ] = mrate_anomaly*rho_ice
        bmb_var[k,::] += bmelt_anomaly
        #bmb_var[k,::][ larmip_regs==larmip_region ] += mrate_anomaly*rho_ice
    if smb_background.units=='kg m-2 s-1':
      smb_var[k,::]*=secperyear
        
nc.close()
nc_a.close()
nc_b.close()
