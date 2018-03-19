import numpy as np
import netCDF4 as nc
import os, sys, glob, csv
import subprocess

## this hack is needed to import config.py from the project root
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if project_root not in sys.path: sys.path.append(project_root)
import config as cf; reload(cf)
#import pism_input.pism_input as pi; reload(pi)

dataset = "raised"
data_path = os.path.join(cf.output_data_path, dataset)

if not os.path.exists(data_path): os.makedirs(data_path)


sourcefile_pattern = os.path.join(cf.raised_data_path,'raised_*1kmmask.txt')

for infile in glob.glob(sourcefile_pattern):

  ncfile = infile.replace('.txt','.nc')
  ncfile = os.path.join(data_path,os.path.basename(ncfile))

  ## get data
  reader = csv.reader(open(infile, "rb"))
  datainfo={}

  for row, record in enumerate(reader):
    if row<6: #header
      dd = np.array(record[0].split())
      datainfo[dd[0]]=np.float(dd[1])

  x = np.arange(datainfo['xllcorner'],datainfo['xllcorner']+datainfo['ncols']*datainfo['cellsize'],datainfo['cellsize']) # 1km grid
  y = np.arange(datainfo['yllcorner'],datainfo['yllcorner']+datainfo['nrows']*datainfo['cellsize'],datainfo['cellsize'])
  print infile, datainfo

  #print x,len(x)
  reader = csv.reader(open(infile, "rb"))
  data = np.empty([len(y),len(x)])
  k=0
  stackdd = np.array([])
  for row, record in enumerate(reader):
    if row>=6:
      #dd = np.array(record[0].split()).astype(np.float)
      dd=np.array(list(record[0])).astype(np.float)
      stackdd = np.hstack([stackdd,dd])
      data[k,:] = stackdd
      k += 1
      stackdd = np.array([])

  data=np.flipud(data)

  print "create netcdf file " + ncfile
  ncout = nc.Dataset(ncfile, 'w', format='NETCDF3_CLASSIC')
  ncout.createDimension('time',size=None)
  ncout.createDimension('x',size=len(x))
  ncout.createDimension('y',size=len(y))
  ncx = ncout.createVariable(varname="x",datatype='float32',dimensions=('x'))
  ncy = ncout.createVariable(varname="y",datatype='float32',dimensions=('y'))
  ncx[:] = x
  ncy[:] = y
  ncx.units = "m"
  ncy.units = "m"

  ncv = ncout.createVariable( varname="mask",datatype='byte',dimensions=('y','x') )
  ncv[:] = data

  ncout.projection = "+proj=stere +lon_0=0 +lat_0=-90 +lat_ts=-71 +ellps=WGS84 +datum=WGS84"
  ncout.data_source = "data based on M.J. Bentley et al., (2014). A community-based geological reconstruction of Antarctic Ice Sheet deglaciation since the Last Glacial Maximum. Quaternary Science Reviews, 100, 1-9."

  print "close"
  ncout.close()
