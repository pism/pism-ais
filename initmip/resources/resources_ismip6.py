#!/usr/bin/env python
# Copyright (C) 2016 Andy Aschwanden (https://github.com/pism/pism-gris/blob/master/initMIP)
# Modified by Torsten Albrecht for initMIP Antarctica

from netCDF4 import Dataset as CDF
import csv
import cf_units
import numpy as np
from pyproj import Proj


class ISMIP6Var(object):
    ismip6_name = None
    pism_name = None
    units = None
    standard_name = None
    def __init__(self, ismip6_name, pism_name, units, standard_name, state, do_mask):
        self.ismip6_name = ismip6_name
        self.pism_name = pism_name
        self.units = units
        self.standard_name = standard_name
        self.state = state
        self.do_mask = do_mask

    def __repr__(self):
        return "ISMIP6 Variable"


pism_mass_to_vol_vars = (
    'discharge_flux'
    )

def get_ismip6_vars_dict(file, dim):
    '''
    Returns a dictionary with ISMIP6 Variables
    '''
        
    ismip6_vars = {}
    with open(file, 'r') as f:
        has_header = csv.Sniffer().has_header(f.read())
        f.seek(0)  # rewind
        incsv = csv.reader(f)
        if has_header:
            next(incsv)  # skip header row
        for row in csv.reader(f, skipinitialspace=True):
            m_list = [x.strip() for x in row]
            ismip6_name = m_list[0]
            pism_name = m_list[1]
            units = m_list[2]
            standard_name = m_list[3]
            dimension = int(m_list[4])
            state = int(m_list[5])
            do_mask = int(m_list[6])
            if dimension == dim:
                ismip6_vars[ismip6_name] = ISMIP6Var(ismip6_name, pism_name, units, standard_name, state, do_mask)
        
    return ismip6_vars


def adjust_time_axis(icesheet,filename):
    '''
    Adjusts the time axis
    '''
    
    nc = CDF(filename, 'a')
    time = nc.variables['time']
    time_bnds_var = time.bounds
    time_bnds = nc.variables[time_bnds_var]
    nt = len(time[:])
    new_timeline = np.linspace(0, 100, nt, endpoint=True)
    time[:] = new_timeline
    time_bnds[:,0] = new_timeline - 5
    time_bnds[:,1] = new_timeline
    if icesheet=="GIS":
      time.units = 'years since 2008-1-1'
    elif icesheet=="AIS":
      time.units = 'years since 2000-1-1'

    nc.close()


def make_spatial_vars_ismip6_conforming(filename, ismip6_vars_dict):
    '''
    Make file ISMIP6 conforming
    '''
    
    # Open file
    nc = CDF(filename, 'a')

    cell_area_var = nc.variables['cell_area']
    cell_area_units = cell_area_var.units
    cell_area = cell_area_var[:]

    pism_to_ismip6_dict = dict((v.pism_name, k) for k, v in ismip6_vars_dict.iteritems())
    
    for pism_var in nc.variables:
        nc_var = nc.variables[pism_var]
        if pism_var in pism_to_ismip6_dict.keys():
            ismip6_var = pism_to_ismip6_dict[pism_var]
            print('Processing {} / {}'.format(pism_var, ismip6_var))
            if not pism_var == ismip6_var:
                print('  Renaming {pism_var} to {ismip6_var}'.format(pism_var=pism_var, ismip6_var=ismip6_var))
                nc.renameVariable(pism_var, ismip6_var)
                nc.sync()
            if pism_var in pism_mass_to_vol_vars:
                print('  Converting {pism_var} from mass to volume'.format(pism_var=pism_var))
                nc_var[:] /= cell_area
                i_units = nc_var.units
                o_units = cf_units.Unit(i_units) / cf_units.Unit(cell_area_units)
                nc_var.units = o_units.format()
            if not nc_var.units == ismip6_vars_dict[ismip6_var].units:
                o_units = ismip6_vars_dict[ismip6_var].units            
                i_units = nc_var.units
                print('  Converting {pism_var} from {i_units} to {o_units}'.format(pism_var=pism_var, i_units=i_units, o_units=o_units))    
                i_f = cf_units.Unit(i_units)
                o_f = cf_units.Unit(o_units)
                nc_var[:] = i_f.convert(nc_var[:], o_f)
                nc_var.units = o_units
                nc_var.standard_name = ismip6_vars_dict[ismip6_var].standard_name
    nc.close()


def make_scalar_vars_ismip6_conforming(filename, ismip6_vars_dict):
    '''
    Make file ISMIP6 conforming
    '''
    
    # Open file
    nc = CDF(filename, 'a')

    pism_to_ismip6_dict = dict((v.pism_name, k) for k, v in ismip6_vars_dict.iteritems())
    
    for pism_var in nc.variables:
        nc_var = nc.variables[pism_var]
        if pism_var in pism_to_ismip6_dict.keys():
            ismip6_var = pism_to_ismip6_dict[pism_var]
            print('Processing {} / {}'.format(pism_var, ismip6_var))
            if not pism_var == ismip6_var:
                print('  Renaming {pism_var} to {ismip6_var}'.format(pism_var=pism_var, ismip6_var=ismip6_var))
                nc.renameVariable(pism_var, ismip6_var)
                nc.sync()
            if not nc_var.units == ismip6_vars_dict[ismip6_var].units:
                o_units = ismip6_vars_dict[ismip6_var].units            
                i_units = nc_var.units
                print('  Converting {pism_var} from {i_units} to {o_units}'.format(pism_var=pism_var, i_units=i_units, o_units=o_units))    
                i_f = cf_units.Unit(i_units)
                o_f = cf_units.Unit(o_units)
                nc_var[:] = i_f.convert(nc_var[:], o_f)
                nc_var.units = o_units
                nc_var.standard_name = ismip6_vars_dict[ismip6_var].standard_name
    nc.close()


def create_searise_grid(icesheet,filename, grid_spacing, **kwargs):
    '''
    Create dummy grid description
    '''

    if 'fileformat' not in kwargs.keys():
        fileformat = 'NETCDF4'
    else:
        fileformat = str.upper(kwargs['fileformat'])
        
    
    xdim = 'x'
    ydim = 'y'

    if icesheet=="GIS":

      # define output grid, these are the extents of the Bamber domain
      e0 = -800000.0
      n0 = -3400000.0
      e1 = 700000.0
      n1 = -600000.0

      # # Shift to cell centers
      # e0 += grid_spacing / 2
      # n0 += grid_spacing / 2
      # e1 -= grid_spacing / 2
      # n1 -= grid_spacing / 2

      projection = "+proj=stere +ellps=WGS84 +datum=WGS84 +lon_0=-39 +lat_0=90 +lat_ts=71 +units=m"


    elif icesheet=="AIS":

      e0 = -3040000.0
      n0 = -3040000.0
      e1 = 3040000.0
      n1 = 3040000.0

      projection = "+proj=stere +ellps=WGS84 +datum=WGS84 +lon_0=0 +lat_0=-90 +lat_ts=-71 +units=m" # +x_0=0.0 +y_0=0.0"



    de = dn = grid_spacing  # m
    M = int((e1 - e0) / de) + 1
    N = int((n1 - n0) / dn) + 1

    easting = np.linspace(e0, e1, M)
    northing = np.linspace(n0, n1, N)
    ee, nn = np.meshgrid(easting, northing)

    # Set up SeaRISE Projection
    #projection = "+proj=stere +ellps=WGS84 +datum=WGS84 +lon_0=-39 +lat_0=90 +lat_ts=71 +units=m"
    proj = Proj(projection)

    lon, lat = proj(ee, nn, inverse=True)

    # number of grid corners
    grid_corners = 4
    # grid corner dimension name
    grid_corner_dim_name = "nv4"

    # array holding x-component of grid corners
    gc_easting = np.zeros((M, grid_corners))
    # array holding y-component of grid corners
    gc_northing = np.zeros((N, grid_corners))
    # array holding the offsets from the cell centers
    # in x-direction (counter-clockwise)
    de_vec = np.array([-de / 2, de / 2, de / 2, -de / 2])
    # array holding the offsets from the cell centers
    # in y-direction (counter-clockwise)
    dn_vec = np.array([-dn / 2, -dn / 2, dn / 2, dn / 2])
    # array holding lat-component of grid corners
    gc_lat = np.zeros((N, M, grid_corners))
    # array holding lon-component of grid corners
    gc_lon = np.zeros((N, M, grid_corners))
    
    for corner in range(0, grid_corners):
        ## grid_corners in x-direction
        gc_easting[:, corner] = easting + de_vec[corner]
        # grid corners in y-direction
        gc_northing[:, corner] = northing + dn_vec[corner]
        # meshgrid of grid corners in x-y space
        gc_ee, gc_nn = np.meshgrid(
            gc_easting[:, corner], gc_northing[:, corner])
        # project grid corners from x-y to lat-lon space
        gc_lon[:, :, corner], gc_lat[:, :, corner] = proj(
            gc_ee, gc_nn, inverse=True)


    nc = CDF(filename, 'w', format=fileformat)

    nc.createDimension(xdim, size=easting.shape[0])
    nc.createDimension(ydim, size=northing.shape[0])
    
    var = xdim
    var_out = nc.createVariable(var, 'f', dimensions=(xdim))
    var_out.axis = xdim
    var_out.long_name = "X-coordinate in Cartesian system"
    var_out.standard_name = "projection_x_coordinate"
    var_out.units = "meters"
    var_out[:] = easting

    var = ydim
    var_out = nc.createVariable(var, 'f', dimensions=(ydim))
    var_out.axis = ydim
    var_out.long_name = "Y-coordinate in Cartesian system"
    var_out.standard_name = "projection_y_coordinate"
    var_out.units = "meters"
    var_out[:] = northing

    var = 'lon'
    var_out = nc.createVariable(var, 'f', dimensions=(ydim, xdim))
    var_out.units = "degrees_east"
    var_out.valid_range = -180., 180.
    var_out.standard_name = "longitude"
    var_out.bounds = "lon_bnds"
    var_out[:] = lon

    var = 'lat'
    var_out = nc.createVariable(var, 'f', dimensions=(ydim, xdim))
    var_out.units = "degrees_north"
    var_out.valid_range = -90., 90.
    var_out.standard_name = "latitude"
    var_out.bounds = "lat_bnds"
    var_out[:] = lat

    nc.createDimension(grid_corner_dim_name, size=grid_corners)

    var = 'lon_bnds'
    # Create variable 'lon_bnds'
    var_out = nc.createVariable(
        var, 'f', dimensions=(ydim, xdim, grid_corner_dim_name))
    # Assign units to variable 'lon_bnds'
    var_out.units = "degreesE"
    # Assign values to variable 'lon_nds'
    var_out[:] = gc_lon
        
    var = 'lat_bnds'
    # Create variable 'lat_bnds'
    var_out = nc.createVariable(
        var, 'f', dimensions=(ydim, xdim, grid_corner_dim_name))
    # Assign units to variable 'lat_bnds'
    var_out.units = "degreesN"
    # Assign values to variable 'lat_bnds'
    var_out[:] = gc_lat

    var = 'dummy'
    var_out = nc.createVariable(
        var,
        'f',
        dimensions=(
            "y",
            "x"),
        fill_value=-2e9)
    var_out.units = "meters"
    var_out.long_name = "Just A Dummy"
    var_out.comment = "This is just a dummy variable for CDO."
    var_out.grid_mapping = "mapping"
    var_out.coordinates = "lon lat"
    var_out[:] = 0.

    mapping = nc.createVariable("mapping", 'c')
    mapping.ellipsoid = "WGS84"
    mapping.false_easting = 0.
    mapping.false_northing = 0.
    mapping.grid_mapping_name = "polar_stereographic"

    if icesheet=="GIS":
          mapping.latitude_of_projection_origin = 90.
          mapping.standard_parallel = 71.
          mapping.straight_vertical_longitude_from_pole = -39.
    elif icesheet=="AIS":
          mapping.latitude_of_projection_origin = -90.
          mapping.standard_parallel = -71.
          mapping.straight_vertical_longitude_from_pole = 0.

    from time import asctime
    historystr = 'Created ' + asctime() + '\n'
    nc.history = historystr
    nc.proj4 = projection
    nc.Conventions = 'CF-1.5'
    nc.close()
