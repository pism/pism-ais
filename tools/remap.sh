#!/bin/bash



###################################################################
scriptname="remap.sh"

if [ $# -gt 0 ] ; then  # dataset
  dataset="$2"
  resolution="$1"
  NN="$3"
  #echo "arguments $dataset $resolution $NN"
fi


###### use only MPI if job is submitted
if [ -n "${PISM_ON_CLUSTER:+1}" ]; then  # check if env var is set
  echo "$scriptname: this run was submitted, use MPI"
  echo
  export I_MPI_PMI_LIBRARY=/p/system/slurm/lib/libpmi.so
  #MPIDO="mpiexec.hydra -bootstrap slurm -n"
  MPIDO="srun -n $NN"
  CDOP=" -P $NN"
  echo "$scriptname: MPIDO = $MPIDO"
else
  echo "$scriptname: this is interactive, skip use of MPI"
  echo
  MPIDO=""
  CDOP=""
fi

###################################################################
###################################################################
res=$resolution
data=$dataset

export workpath=/p/projects/tumble/pism_input/GitLab
#export workpath=/home/albrecht/Documents/pism/python/pism_input
#export workpath=../

#tools on PIK cluster
export cdopath=/p/system/packages/cdo/1.7.1/bin
export ncopath=/p/system/packages/nco/4.5.0/bin
#export cdopath=/usr/bin
#export ncopath=/usr/bin

export inputfile=${workpath}/${data}/${data}_data/${data}_1km_input.nc
python ${workpath}/tools/nc2cdo.py $inputfile

mkdir -p ${workpath}/${data}/${data}_weights
echo $inputfile ${workpath}/${data}/${data}_weights

#for res in 50 30 20 15 12 10 7 5 3 2 1;
#do

  #export targetgrid=${workpath}/grids/pism_${res}.0km.nc #albmap grid
  export targetgrid=${workpath}/grids/pism_${res}km.nc #albmap grid
  export mapres=${workpath}/${data}/${data}_data/${data}_${res}km.nc
  export mapweights=${workpath}/${data}/${data}_weights/${data}_${res}km_weights.nc

  echo "Run remapycon and create $mapres"
  echo

  if [ "${data}" == bedmap2 -o albmap ]; then
    $MPIDO ${cdopath}/cdo $CDOP genycon,${targetgrid} ${inputfile} ${mapweights}
    $MPIDO ${cdopath}/cdo $CDOP -b F64 remap,${targetgrid},${mapweights} ${inputfile} ${mapres}
 
  else #velocities etc...
    export REMAP_EXTRAPOLATE=on
    $MPIDO ${cdopath}/cdo $CDOP -b F64 remapbil,${targetgrid} ${inputfile} ${mapres}
 
  fi

  $MPIDO ${ncopath}/ncks -A -v x,y ${targetgrid} ${mapres}

#done


