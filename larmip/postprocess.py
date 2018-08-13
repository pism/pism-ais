# Created by julius.garbe@pik-potsdam.de for LARMIP Antarctica 2018

import os
import numpy as np
import netCDF4 as cdf


def get_time():
    ncf = cdf.Dataset(ctrl_file,"r")
    #time = ncf.variables['time'][4:] # crop first 4 years (start at year 2000)
    time = ncf.variables['time'][:] #(start at year 2000)
    time = [y*1./(60*60*24*365) for y in time] # convert time: seconds -> years
    # time = [y+2000.0 for y in time]
    return time

def get_slvol(filename):
    try:
        ncf = cdf.Dataset(filename,"r")
    except IOError as error:
        print filename, "not found."
        raise error
    return ncf.variables['slvol'][:] # (start at year 2000)
    #return ncf.variables['slvol'][4:] # crop first 4 years (start at year 2000)


def compute_anomaly(data):
    return -(data - data[0])

####

def process_slvol(infile,region,drift_corr=True):
    
  print infile
    
  if os.path.exists(infile):    

    ## get slvol variable
    slvol = get_slvol(infile)
    slvol_diff = compute_anomaly(slvol)
    
    ## drift correction
    if drift_corr:
        slvol_ctrl_diff = compute_anomaly(get_slvol(ctrl_file))
        slvol_diff = slvol_diff - slvol_ctrl_diff
    
    ## write variable to file
    varname = region #'sl change from '+region
    var_out = nc.createVariable(varname, 'float64', dimensions=("time"))
    var_out.units = "meters";
    if region=="CTRL":
        var_out.long_name = "sea-level contribution from "+region;
    elif region=="all":
        var_out.long_name = "sea-level contribution from "+region+" regions relative to CTRL";
    else:
        var_out.long_name = "sea-level contribution from "+region+" relative to CTRL";
    #var_out.long_name = "total sea-level relevant ice volume change IN SEA-LEVEL EQUIVALENT";
    var_out[:] = slvol_diff
    
    # print str(np.around(slvol_diff[-1], decimals=3) * 100),"cm"
    return slvol_diff

##################################################################

case='highres'

if case=='equi':

  project_dir = '/home/garbe/PIK_Tmp/LARMIP/'
  experiment_dir = 'pismpik_ant_forcing_larmip8km_58kyr/'
  timeseries_dir = os.path.join(project_dir, experiment_dir, "timeseries")
  res='8'
  dur='200'
  # output_dir = '/p/tmp/albrecht/pism18/pismOut/larmip/larmip2301b/PISMEQUI'+res
  output_dir = "/home/garbe/PIK_Pism/garbe/larmip/postprocess/data/PISMEQUI"+res

elif case=='paleo':

  project_dir = '/p/tmp/albrecht/pism18/pismOut/larmip/'
  experiment_dir = 'larmip2301b/'
  timeseries_dir = os.path.join(project_dir, experiment_dir, "results")
  res='16'
  dur='205'
  # output_dir = "/home/garbe/PIK_Pism/garbe/larmip/postprocess/data/PISMPAL"+res
  output_dir = os.path.join(project_dir, experiment_dir, "PISMPAL"+res)

elif case=='highres':

  project_dir = '/p/projects/pism/albrecht/larmip/'
  experiment_dir = 'f2543k/'
  timeseries_dir = os.path.join(project_dir, experiment_dir, "results")
  res='4'
  dur='200'
  output_dir = os.path.join(project_dir, experiment_dir, "PISMHRES"+res)


os.system("mkdir -p "+output_dir)


#### Changes should not be necessary in this section

#experiments = ["ADD0","ADD1","ADD2","ADD3","ADD4","ADD5","ADD6"]
#experiments = ["ADD4"]
experiments = ["ADD2","ADD3","ADD4","ADD5","ADD6"]


#forcings = ['2','1','4','2','8','32','16'] # do not change!
#forcings = ['8']
forcings = ['4','2','8','32','16']

regions = ["reg1","reg2","reg3","reg4","reg5"]

ctrl_file = os.path.join(timeseries_dir, "ts_control_"+res+"km_"+dur+"yrs.nc")

time = get_time()
ltime = len(time)
####


for i,experiment in enumerate(experiments):
        
    outfile = os.path.join(output_dir, experiment+".nc")
    
    ### FIX-ME: delete existing file first to avoid errors
    os.system("rm -f "+outfile)
    
    nc = cdf.Dataset(outfile,'w',format='NETCDF4_CLASSIC')
    
    print "EXPERIMENT:",experiment
    
    ## create time dimension
    nc.createDimension("time", size=ltime)

    ## add time variable
    var_out_time = nc.createVariable('time', 'float64', dimensions=("time"))
    var_out_time[:] = time
    var_out_time.units = 'years since 2000-01-01 00:00:00'
    
    print "Reading slvol time series from ..."
    
    ## add CTRL variable
    process_slvol(ctrl_file,"CTRL",drift_corr=False)
    

    if experiment=="ADD3":
        for f in forcings:
          infile = os.path.join(timeseries_dir, "ts_abmb_m"+f+"_all_"+res+"km_"+dur+"yrs.nc")
          process_slvol(infile,f,drift_corr=True)
    
    
    else:
        for region in regions:
            infile = os.path.join(timeseries_dir, "ts_abmb_m"+forcings[i]+"_"+region+"_"+res+"km_"+dur+"yrs.nc")
            process_slvol(infile,region,drift_corr=True)
    
            
    
    ## adjust global attributes
    from time import asctime
    nc.history = 'Created on ' + asctime() + ' by Torsten Albrecht, PIK'
    nc.institution = 'Potsdam Institute for Climate Impact Research (PIK), Germany'
    nc.contact = 'torsten.albrecht@pik-potsdam.de and julius.garbe@pik-potsdam.de'
    #nc.source = 'PISM stable 1.0 (https://github.com/talbrecht/pism_pik; branch: pism_pik_1.0_larmip; commit: 3172d9b)'
    nc.source = 'PISM stable 1.0 (https://github.com/talbrecht/pism_pik; branch: dev_larmip; commit: 4402963f)'    
    nc.Conventions = 'CF-1.6'
    
    nc.close()
    print "NetCDF4 file written to",outfile,"\n"
