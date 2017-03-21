#!/bin/bash

#This loop calls the create_grid.py script for different resolutions and adds lon/lat information
###################################################################


for res in 50 30 20 15 12 10 7 5 3 2 1;
do

  python create_grid.py -g ${res} pism_${res}km.nc

  python ../tools/nc2cdo.py pism_${res}km.nc

done
