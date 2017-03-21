#!/bin/bash

#SBATCH --qos=medium
#SBATCH --job-name=cdo_remap_
#SBATCH --account=ice
#SBATCH --output=log/remapcony-%j.out
#SBATCH --error=log/remapcony-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --profile=energy
#SBATCH --acctg-freq=energy=5
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=albrecht@pik-potsdam.de


#scriptname=remap.sh
scriptname=/p/projects/tumble/pism_input/GitLab/tools/remap.sh

#calls for example "${scriptname} bedmap2 15 1"

if [ $# -gt 0 ] ; then  # dataset
  dataset="$1"
  resolution="$2"
  inputres="$3"
fi


export PISM_ON_CLUSTER=1
echo "${scriptname} ${dataset} ${resolution} ${inputres} $SLURM_NTASKS"
${scriptname} ${dataset} ${resolution} ${inputres} $SLURM_NTASKS >> ./log/remap.out



