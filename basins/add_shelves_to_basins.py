import os, glob
import datetime
import csv
import math
import numpy as np
import numpy.ma as ma
import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pylab as p
from matplotlib.patches import Rectangle
plt.rcParams['font.size']=30
plt.rcParams['legend.fontsize']=30
plt.rcParams['figure.figsize'] = 15, 15
plt.close('all')

"""
Parameters
"""
getData=1

ncfile = "basins_data/basins_zwally_5km_shelves.nc"
basinfile = "basins_data/basins_zwally_5km_input.nc"
maskfile = "../albmap/albmap_data/Antarctica_5km_dev1.0.nc"

print '\nReading mask... '
infile = nc.Dataset(maskfile, 'r')
mask = ma.masked_array(infile.variables['mask'][0,:,:]) 
infile.close()
Mx=np.shape(mask)[0]
My=np.shape(mask)[1]
mask_floating = 2 
mask_grounded = 1  

print 'Mx = ' + str(Mx)
print 'My = ' + str(My)

print '\nReading basins... '
infile = nc.Dataset(basinfile, 'r')
basins_orig = ma.masked_array(infile.variables['basins'][:,:]) 
infile.close()

basins=basins_orig
numberOfBasins = int(basins.max())

"""
Basin id for grounded cells is already set, find the basin id for shelves recursively:
"""
if (getData==1):
    done = 0
    loopcount = 0

    while(done == 0):  
	done = 1
	for i in range(1,Mx-1):
	    for j in range(1,My-1):
		if (mask[i,j] > 0.0 and basins[i,j] == 0.0):
		    done = 0
		    #print i,j
		    for basinid in range(1,numberOfBasins+1):
			#print '\nBasin: ' + str(basinid)
			if (basins[i+1,j] == basinid or basins[i-1,j] == basinid or basins[i,j+1] == basinid or basins[i,j-1] == basinid or basins[i+1,j-1] == basinid or basins[i+1,j+1] == basinid or basins[i-1,j-1] == basinid or basins[i-1,j+1] == basinid):
			    basins[i,j] = basinid
		#if (i==350 and j== 550):
		    #print 'mask = ' + str(mask[i,j]) + ', basinid = ' + str(basins[i,j])
	loopcount = loopcount+1
	print '...loop ' + str(loopcount)
	if (loopcount == 2):
	    done = 1

    
print '\nPlotting... '
fig=plt.figure()
ax = plt.subplot(1,1,1) 
plt.contourf(basins) 

plt.contour(mask,mask_grounded,colors='gray',linewidths=2) # GROUNDING LINE 
plt.contour(mask,mask_floating,colors='gray',linewidths=2) # CALVING FRONT 
#cbar = plt.colorbar(ticks=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19])

ax.set_xticks([]) 
ax.set_yticks([]) 
p.axis('equal')

#plt.show()
os.system("mkdir plots")
plt.savefig('plots/ZwallyBasins.png')
plt.clf()


"""
Create nc-file
"""
print "create netcdf file\n" + ncfile
ncout = nc.Dataset(ncfile, 'w', format='NETCDF4_CLASSIC')
#ncout.createDimension('time',size=None)
ncout.createDimension('x',size=1200)#len(x))
ncout.createDimension('y',size=1200)#len(y))
ncx = ncout.createVariable(varname="x",datatype='float32',dimensions=('x'))
ncy = ncout.createVariable(varname="y",datatype='float32',dimensions=('y'))
ncx[:] = np.arange(-2800e3,3200e3,5000)#3195e3,5000)
ncy[:] = np.arange(-2800e3,3200e3,5000)#3195e3,5000)
ncx.units = "m"
ncy.units = "m"

ncv = ncout.createVariable( varname="basins",datatype='float32',dimensions=('y','x') )
ncv[:] = basins

##ncv = {}
#for i in range(len(basinnames)):
  #varname = basinnames[i] + "_" + str(basinnums[i]).zfill(2)
  #ncv = ncout.createVariable( varname=varname,datatype='float32',dimensions=('y','x') )
  #ncv[:] = basins
  #print ncv
  #ncout.sync()

now = datetime.datetime.now().strftime("%B %d, %Y")
ncout.created  = "created based on reese@pik at " + now
ncout.data_origin = "downloaded from http://homepages.see.leeds.ac.uk/~earkhb/Basins_page.html"
ncout.regions = "regions follow the definition of Zwally et al. 2012"

print "close"
ncout.close()
