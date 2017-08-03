#!/bin/bash

#15 km to 16km
#pism_init_file=/p/tmp/albrecht/pism17/pismOut/initmip/initmip2200/initdata/pism_bedmap2_racmo_uplift_velrignot_lgmokill_15km.nc
#pism_exp_outfile=/p/tmp/albrecht/pism17/pismOut/initmip/initmip2200/results/extra_mbcotrl_15km_105yrs.nc
#pism_exp_file=/p/tmp/albrecht/pism17/pismOut/initmip/initmip2200/results/extra_mbcotrl_15km_105yrs_ismip6.nc
#pism_ts_outfile=/p/tmp/albrecht/pism17/pismOut/initmip/initmip2200/results/ts_mbcotrl_15km_105yrs.nc

#python postprocessing_ex.py ${pism_init_file} ${pism_exp_outfile}  -n 1 -e ctrl -r bil -t 16000 --id 5
#-w searise_grid_16000m_bil_weights.nc

# all 16km grid
pism_run=initmip2207
pism_init_file=/p/tmp/albrecht/pism17/pismOut/initmip/${pism_run}/initdata/bedmap2_albmap_racmo_hadcm3_I2S_schmidtko_uplift_lgmokill_fttmask_16km.nc

for res in 16000 # 32000 8000
do
  for exp in asmb abmb ctrl
  do
  
    export cdo_weight=""
    if [ "$exp" == ctrl ]
    then
      echo "Create cdo weight file"
      cdo_weight="-w"
    fi

    pism_exp_outfile=/p/tmp/albrecht/pism17/pismOut/initmip/${pism_run}/results/extra_${exp}_16km_105yrs.nc
    pism_ts_outfile=/p/tmp/albrecht/pism17/pismOut/initmip/${pism_run}/results/ts_${exp}_16km_105yrs.nc
    echo $exp $pism_exp_outfile

    if [ $res == 16000 ]
    then
      echo "Postprocess timeseries"
      python postprocess_ts.py ${pism_ts_outfile} -e ${exp} -t $res --id 3
    fi

    python postprocess_ex.py ${pism_init_file} ${pism_exp_outfile} ${cdo_weight} -n 4 -e ${exp} -r ycon -t $res --id 3

  done
done

