#!/bin/bash

#SBATCH --qos=medium
#SBATCH --job-name=cdo_remap_bedmap2
#SBATCH --account=ice
#SBATCH --output=log/bedmap2-%j.out
#SBATCH --error=log/bedmap2-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --profile=energy
#SBATCH --acctg-freq=energy=5
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=albrecht@pik-potsdam.de


###################################################################
NN=2  # default number of processors
if [ $# -gt 0 ] ; then  # if user says "remap.sh 8" then NN = 8
  NN="$1"
fi

scriptname='remap_bedmap.sh'

###### use only MPI if job is submitted
if [ -n "${PISM_ON_CLUSTER:+1}" ]; then  # check if env var is set
  echo "$scriptname: this run was submitted, use MPI"
  export I_MPI_PMI_LIBRARY=/p/system/slurm/lib/libpmi.so
  #MPIDO="mpiexec.hydra -bootstrap slurm -n"
  MPIDO="srun -n $NN"
  echo "$scriptname: MPIDO = $MPIDO"
else
  echo "$scriptname: this is interactive, skip use of MPI"
  MPIDO=""
  NN=1
fi

###################################################################
###################################################################

data=bedmap2
#export workpath=/p/projects/tumble/pism_input
export workpath=/home/albrecht/Documents/pism/python/pism_input
#export workpath=../

#tools
#export cdopath=/p/system/packages/cdo/1.7.1/bin
#export ncopath=/p/system/packages/nco/4.5.0/bin
export cdopath=/usr/bin
export ncopath=/usr/bin

export inputfile=${workpath}/${data}/${data}_data/${data}_1km.nc
python ${workpath}/tools/nc2cdo.py $inputfile

mkdir -p ${workpath}/${data}/${data}_weights

#for res in 50 30 20 15 12 10 7 5 3 2;
for res in 15 12;
do
  export targetgrid=${workpath}/grids/pism_${res}.0km.nc #albmap grid
  export mapres=${workpath}/${data}/${data}_data/${data}_${res}km.nc
  export mapweights=${workpath}/${data}/${data}_weights/${data}_${res}km_weights.nc

  echo "remapycon $mapres"
  echo "$MPIDO ${cdopath}/cdo -P $NN genycon,${targetgrid} ${inputfile} ${mapweights}"
  $MPIDO ${cdopath}/cdo -P $NN genycon,${targetgrid} ${inputfile} ${mapweights}
  $MPIDO ${cdopath}/cdo -P $NN remap,${targetgrid},${mapweights} ${inputfile} ${mapres}

  #$MPIDO ${cdopath}/cdo -P 8 remapbil,${targetgrid} ${inputfile} ${mapres}

  $MPIDO ${ncopath}/ncpdq -O --fl_fmt=64bit ${mapres} ${mapres}
  $MPIDO ${ncopath}/ncks -A -v x,y ${targetgrid} ${mapres}

done

    

