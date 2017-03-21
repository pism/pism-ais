#!/bin/bash

#res=2km #3000x3000
#res=10km #600X600
#res=12km #500X500
#res=15km #400x400
#res=20km #300x300
#res=30km #200x200
res=50km #120x120


# combine bedmap2 and albmap atmosphere date
pisminpath=/p/projects/tumble/pism_input/GitLab
albmapfile=${pisminpath}/albmap/albmap_data/albmap_${res}.nc
bedmapfile=${pisminpath}/bedmap2/bedmap2_data/bedmap2_${res}.nc
pismgridfile=${pisminpath}/GitLab/grids/pism_${res}.nc

pismoutfile=pism_bedmap2_albmap_${res}.nc
ncks -A -v topg,thk,usurf,bedunc,x,y $bedmapfile -o $pismoutfile
ncks -A -v air_temp,bheatflx,precipitation $albmapfile $pismoutfile

ncatted -O -a proj4,global,o,c,"+lon_0=0.0 +ellps=WGS84 +datum=WGS84 +lat_ts=-71.0 +proj=stere +x_0=0.0 +units=m +y_0=0.0 +lat_0=-90.0" $pismoutfile


