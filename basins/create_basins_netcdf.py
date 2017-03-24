import numpy as np
import csv
import datetime
import netCDF4 as nc

infile = "./orig_data/ais_basins_imbie.z12.txt"
namefile = "./orig_data/ais_basins_imbie.z12.names.txt"
ncfile = "basiins_data/basins_zwally_5km_input.nc"

## get data
reader = csv.reader(open(infile, "rb"))

x = np.arange(-2800e3,2805e3,5000) # Bamber DEM grid
y = np.arange(-2800e3,2805e3,5000)

data = np.empty([len(x),len(y)])
k=0
stackdd = np.array([])
for row, record in enumerate(reader):
  dd = np.array(record[0].split()).astype(np.float)
  stackdd = np.hstack([stackdd,dd])
  if len(dd) != 9:
    data[k,:] = stackdd
    k += 1
    stackdd = np.array([])

## get names
namereader = csv.reader(open(namefile, "rb"))
basinnames, basinnums = [], []

for row in namereader:
  basinnum, basinname = row[0].split()[0:2]
  print basinnum, basinname
  basinnames.append(basinname)
  basinnums.append(np.int(basinnum))

# add ocean cells to fix resolution
data2 = np.zeros([1200,1200])
for xi in range(len(x)):
	for yi in range(len(y)):
		data2[xi,yi] = data[xi,yi]


print "create netcdf file\n" + ncfile
ncout = nc.Dataset(ncfile, 'w', format='NETCDF4_CLASSIC')
ncout.createDimension('time',size=None)
ncout.createDimension('x',size=1200)#len(x))
ncout.createDimension('y',size=1200)#len(y))
ncx = ncout.createVariable(varname="x",datatype='float32',dimensions=('x'))
ncy = ncout.createVariable(varname="y",datatype='float32',dimensions=('y'))
ncx[:] = np.arange(-2800e3,3200e3,5000)#3195e3,5000)
ncy[:] = np.arange(-2800e3,3200e3,5000)#3195e3,5000)
ncx.units = "m"
ncy.units = "m"

ncv = ncout.createVariable( varname="basins",datatype='float32',dimensions=('y','x') )
ncv[:] = data2

#ncv = {}
for i in range(len(basinnames)):
  varname = basinnames[i] + "_" + str(basinnums[i]).zfill(2)
  ncv = ncout.createVariable( varname=varname,datatype='float32',dimensions=('y','x') )
  ncv[:] = data2 == basinnums[i]
  print ncv
  ncout.sync()

#for row in namereader:
  #ncv[basinname][:] = data == basinnum

now = datetime.datetime.now().strftime("%B %d, %Y")
ncout.created  = "created based on matthias.mengel@pik by reese@pik-potsdam.de at " + now
ncout.data_origin = "downloaded from http://homepages.see.leeds.ac.uk/~earkhb/Basins_page.html"
ncout.regions = "regions follow the definition of 'Ice Flow of the Antarctic Ice Sheet', Rignot et al. 2011, Science"


print "close"
ncout.close()
