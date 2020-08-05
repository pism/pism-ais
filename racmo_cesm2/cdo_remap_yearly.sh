#!/bin/bash


# This script is specific for the PIK cluster, which uses SLURM for job scheduling.
# and the module system.
#SBATCH --qos=priority
#SBATCH --partition=priority
#SBATCH --job-name=cdo_remap_initmip1km
#SBATCH --account=ice
#SBATCH --output=cdo_remap.out
#SBATCH --error=cdo_remap.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=0-10:30:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=albrecht@pik-potsdam.de



cdocmd="srun -n $SLURM_NTASKS cdo -P $SLURM_NTASKS"
OMP_NUM_THREADS=$SLURM_NTASKS


module load cdo
module load nco


export REMAP_EXTRAPOLATE=on

inputdatadir=/p/projects/tumble/albrecht/pism_input/data
filrondir=/p/projects/tumble/albrecht/python/create_filron_setup
folder=racmo_cesm2_yearly

#$cdocmd -b F64 -f nc4c yearmean $inputdatadir/racmo_cesm2/racmo_cesm2_input_hist.nc $filrondir/racmo_cesm2_split/racmo_cesm2_hist_1950-2014_mean.nc
#$cdocmd -b F64 -f nc4c timmean $inputdatadir/racmo_cesm2/racmo_cesm2_input_hist.nc $filrondir/$folder/racmo_cesm2_hist_1950-2014_mean.nc

$cdocmd -b F64 -f nc4c splityear $inputdatadir/racmo_cesm2/racmo_cesm2_input_hist.nc $filrondir/$folder/racmo_cesm2_hist_

$cdocmd -b F64 -f nc4c genbil,$filrondir/merged/bedmap2_orig_mouginotvel.nc $inputdatadir/racmo_cesm2/racmo_cesm2_input_hist.nc $filrondir/$folder/racmo_cesm2_bil_weights.nc

for i in $(seq 1950 2014); 

  do echo $i; 

  $cdocmd -b F64 -f nc4c remap,$filrondir/merged/bedmap2_orig_mouginotvel.nc,$filrondir/$folder/racmo_cesm2_bil_weights.nc $filrondir/$folder/racmo_cesm2_hist_${i}.nc $filrondir/$folder/racmo_cesm2_hist_${i}_bilx.nc

  $cdocmd -b F64 -f nc4c invertlat $filrondir/$folder/racmo_cesm2_hist_${i}_bilx.nc $filrondir/$folder/racmo_cesm2_hist_${i}_bilx_flipped.nc

  ncks -A -v usurf $filrondir/racmo/racmo_hadcm3_a1b_bil_flipped_filled_usurf.nc $filrondir/$folder/racmo_cesm2_hist_${i}_bilx_flipped.nc

  ncks -A -v x,y $filrondir/merged/bedmap2_orig_mouginotvel.nc $filrondir/$folder/racmo_cesm2_hist_${i}_bilx_flipped.nc

  rm $filrondir/$folder/racmo_cesm2_hist_${i}_bilx.nc

done



