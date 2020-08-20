#!/bin/bash


# This script is specific for the PIK cluster, which uses SLURM for job scheduling.
# and the module system.
#SBATCH --qos=priority
#SBATCH --partition=priority
#SBATCH --time=0-00:30:00
#SBATCH --job-name=cdo_remap_initmip8km
#SBATCH --account=ice
#SBATCH --output=cdo_remap.out
#SBATCH --error=cdo_remap.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=garbe@pik-potsdam.de

cdocmd="srun -n $SLURM_NTASKS cdo -P $SLURM_NTASKS"
OMP_NUM_THREADS=$SLURM_NTASKS


module load cdo
module load nco


export REMAP_EXTRAPOLATE=on

datadir=/p/projects/pism/garbe/2018_PISM_Input_Data
griddatadir=$datadir/cdo_remapgrids
inputdatadir=$datadir/racmo_cesm2/input
outputdatadir=$datadir/racmo_cesm2/initmip8km

# bilinear interpolation, an option if conservative methods fail.
$cdocmd -b F64 -f nc4c remapbil,$griddatadir/initmip8km.nc $inputdatadir/racmo_cesm2_spp585.nc $outputdatadir/racmo_cesm2_spp585_initmip8km.nc


# add x and y variables to output file.
ncks -A -v x,y $griddatadir/initmip8km.nc $outputdatadir/racmo_cesm2_spp585_initmip8km.nc

echo "regridded file is $outputdatadir/racmo_cesm2_spp585_initmip8km.nc"
