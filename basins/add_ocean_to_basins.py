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
mergeBasins=1
findLongitudes=1
adjustAP=1
fillInBasinMask=1

ncfile = "basins_data/basins_zwally_5km_shelves_merged.nc"
ncfile2 = "basins_data/basins_zwally_5km_shelves_and_ocean.nc"
basinfile = "basins_data/basins_zwally_5km_shelves.nc"
maskfile = "../albmap/albmap_data/Antarctica_5km_dev1.0.nc"

print '\nReading mask... '
infile = nc.Dataset(maskfile, 'r')
mask = ma.masked_array(infile.variables['mask'][0,:,:]) 
lon = ma.masked_array(infile.variables['lon'][0,:,:]) 
infile.close()
Mx=np.shape(mask)[0]
My=np.shape(mask)[1]
mask_floating = 2 
mask_grounded = 1  

print 'Mx = ' + str(Mx)
print 'My = ' + str(My)


"""
Merge basins: 
"""
if (mergeBasins==1):
    infile = nc.Dataset(basinfile, 'r')
    basins_orig = ma.masked_array(infile.variables['basins'][:,:]) 
    infile.close()

    basins=basins_orig
    numberOfBasins = int(basins.max())
    for i in range(Mx):
	for j in range(My):
	    #basin 1
	    if (basins[i,j] == 2.0):
		basins[i,j] = 1.0
	    if (basins[i,j] == 3.0):
		basins[i,j] = 1.0
	    #basin 2
	    if (basins[i,j] == 4.0):
		basins[i,j] = 2.0
	    #basin 3
	    if (basins[i,j] == 5.0):
		basins[i,j] = 3.0
	    #basin 4
	    if (basins[i,j] == 6.0):
		basins[i,j] = 4.0
	    #basin 5
	    if (basins[i,j] == 7.0):
		basins[i,j] = 5.0
	    if (basins[i,j] == 8.0):
		basins[i,j] = 5.0
	    #basin 6
	    if (basins[i,j] == 9.0):
		basins[i,j] = 6.0
	    if (basins[i,j] == 10.0):
		basins[i,j] = 6.0
	    if (basins[i,j] == 11.0):
		basins[i,j] = 6.0
	    #basin 7
	    if (basins[i,j] == 12.0):
		basins[i,j] = 7.0
	    #basin 8
	    if (basins[i,j] == 13.0):
		basins[i,j] = 8.0
	    #basin 9
	    if (basins[i,j] == 14.0):
		basins[i,j] = 9.0
	    #basin 10
	    if (basins[i,j] == 15.0):
		basins[i,j] = 10.0
	    #basin 11
	    if (basins[i,j] == 16.0):
		basins[i,j] = 11.0
	    #basin 12
	    if (basins[i,j] == 17.0):
		basins[i,j] = 12.0
	    if (basins[i,j] == 18.0):
		basins[i,j] = 12.0
	    if (basins[i,j] == 19.0):
		basins[i,j] = 12.0
	    #basin 13
	    if (basins[i,j] == 20.0):
		basins[i,j] = 13.0
	    #basin 14
	    if (basins[i,j] == 21.0):
		basins[i,j] = 14.0
	    if (basins[i,j] == 22.0):
		basins[i,j] = 14.0
	    #basin 15
	    if (basins[i,j] == 23.0):
		basins[i,j] = 15.0
	    #basin 16
	    if (basins[i,j] == 24.0):
		basins[i,j] = 16.0
	    #basin 17
	    if (basins[i,j] == 25.0):
		basins[i,j] = 17.0
	    #basin 18
	    if (basins[i,j] == 26.0):
		basins[i,j] = 18.0
	    #basin 19
	    if (basins[i,j] == 27.0):
		basins[i,j] = 19.0

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

    now = datetime.datetime.now().strftime("%B %d, %Y")
    ncout.created  = "created based on reese@pik at " + now
    ncout.data_origin = "downloaded from http://homepages.see.leeds.ac.uk/~earkhb/Basins_page.html"
    ncout.regions = "regions follow the definition of Zwally et al. 2012"

    print "close"
    ncout.close()


'''
Adjust Basin 18/19 (originally 26/27)
'''
print '\nReading basins... '
infile = nc.Dataset(ncfile, 'r')
basins_orig = ma.masked_array(infile.variables['basins'][:,:]) 
infile.close()

basins=basins_orig

if (adjustAP):
	print '\nAdjusting Antarctic Peninsula Basin which was wrong set in add_shelves_to_basins.py'

	for i in range(1, Mx-1):
		for j in range(1,My-1):
			if(basins[i,j]==19.0 and j>=130 and j<=166 and i>=763 and i<=800):
				#if ((i-763.0)/(796-763)> (j-134.0)/(166-134)):
				basins[i,j]=18.0


"""
Find longitude for ice divides
"""
#print '\nReading basins... '
#infile = nc.Dataset(ncfile, 'r')
#basins_orig = ma.masked_array(infile.variables['basins'][:,:]) 
#infile.close()

#basins=basins_orig
numberOfBasins = int(basins.max())

if (findLongitudes == 1):
    print '\nFinding longitude - boundaries... '
    longitudeVec=np.zeros((20))
    for i in range(1,Mx-1):
	for j in range(1,My-1):
	    for basinid in range(2,numberOfBasins+1):
		if (basins[i,j] == basinid and (basins[i+1,j] == basinid-1 or basins[i-1,j] == basinid-1 or basins[i,j+1] == basinid-1 or basins[i,j-1] == basinid-1 or basins[i+1,j-1] == basinid-1 or basins[i+1,j+1] == basinid-1 or basins[i-1,j-1] == basinid-1 or basins[i-1,j+1] == basinid-1) and (mask[i+1,j] == 0.0 or mask[i-1,j] == 0.0 or mask[i,j+1] == 0.0 or mask[i,j-1] == 0.0 or mask[i+1,j-1] == 0.0 or mask[i+1,j+1] == 0.0 or mask[i-1,j-1] == 0.0 or mask[i-1,j+1] == 0.0)):
		  longitudeVec[basinid-2] = lon[i,j]
	    if (basins[i,j] == 1.0 and (basins[i+1,j] == 19 or basins[i-1,j] == 19 or basins[i,j+1] == 19 or basins[i,j-1] == 19 or basins[i+1,j-1] == 19 or basins[i+1,j+1] == 19 or basins[i-1,j-1] == 19 or basins[i-1,j+1] == 19) and (mask[i+1,j] == 0.0 or mask[i-1,j] == 0.0 or mask[i,j+1] == 0.0 or mask[i,j-1] == 0.0 or mask[i+1,j-1] == 0.0 or mask[i+1,j+1] == 0.0 or mask[i-1,j-1] == 0.0 or mask[i-1,j+1] == 0.0)):
	      longitudeVec[-1] = lon[i,j]

    print longitudeVec
    np.save('basins_data/longitudeVec.npy', longitudeVec)



if (fillInBasinMask) :
	print '\nFilling basin mask... '
	longitudeVec = np.load('basins_data/longitudeVec.npy')
	longitudeVec = longitudeVec
	for i in range(1,Mx-1):
    	    for j in range(1,My-1):
		if (basins[i,j] == 0.0):
		    for longitude_i in range(len(longitudeVec)-3): # new! vector too long
			if (longitudeVec[longitude_i]==0.0):
			    longitudeVec[longitude_i] = -180.0
			if (longitudeVec[longitude_i+1]==0.0):
			    longitudeVec[longitude_i+1] = 180.0
			if (lon[i,j] >= longitudeVec[longitude_i] and lon[i,j] < longitudeVec[longitude_i+1]):
			    basins[i,j] = longitude_i+2
		    # fill in basin number 1 which is not possible with the meachnism above (..+2)
		    if (lon[i,j] >= longitudeVec[-3] and lon[i,j] < longitudeVec[0]):
			    basins[i,j] = 1.0
		    # fill in at the border between -180/+180 degrees, here between basin 11 and 12
		    if (lon[i,j] >= -180.0 and lon[i,j] < longitudeVec[11]): 
			    basins[i,j] = 12.0 
		    if (j>=560 and j< 684 and i > 221 and i < 322): # for basin 11
			    basins[i,j] = 11.0
		    if (j>=560 and j< 627 and i < 322): # for basin 12
			if ((j-560)<i/5.0):
			    basins[i,j] = 12.0 
		    if (j>=560 and j< 720 and i <= 221): # rest in this region should be basin 10
			if ((j-560)>=i/5.0):
			    basins[i,j] = 10.0 
		    #Amery wider at bottom
		    if (i>=683 and i< 750 and j > 987): 
			if ( (i-683.0)/(750-683) < (j-987.0)/(1200-987)):
				basins[i,j] = 6.0 # replaced
		    # region of Pine Island Glacier
		    if (i>=482 and i <= 494 and j>= 217 and j<=244):#filling 'hole' in PIG region
			basins[i,j]=14.0 #replaced
		    if (i>=430 and i <= 480 and j>= 170 and j<=211): #extending basin 11 downwards
			if((i-430.0)/(480-430)> (j-170.0)/(211-170)):
				basins[i,j] = 15.0 
		    if (i>=430 and i <= 480 and j<= 170 ):
			basins[i,j] = 15.0 # replaced
		    # Antarctic Peninsula between 17 and 18 		    
		    if (j< 66 and i > 764):
			basins[i,j] = 17.0 	
		    if (j> 65 and j< 167 and i > 764 and basins[i,j]!=17.0):
			basins[i,j] = 18.0 
		    #small insulas which are basin 17 and should be 18...
		    if (j> 156 and j< 166 and  i > 770 and i< 780 and basins[i,j]==17.0):
			basins[i,j] = 18.0 
		    if (j> 75 and j< 190 and  i > 810 and i< 825 and basins[i,j]==17.0):
			basins[i,j] = 18.0 	
		    if (j == 65 and i==875):
			basins[i,j] = 18.0
		    #Adjust Antarctic Peninsula-FRIS-border
		    if (j <= 263 and j >= 167 and i >= 720 ): 
			if ((i-720.0)/(795-720) + (j-167.0)/(263-167)>= 0.5):
				basins[i,j]= 19.0 
		    #Adjust FRIS, small edge which is AP
		    if (j >= 263 and j <= 271 and i >= 715 and i<= 721 ): 
			basins[i,j]= 1.0 



"""
Create nc-file
"""
print "create netcdf file\n" + ncfile2
ncout = nc.Dataset(ncfile2, 'w', format='NETCDF4_CLASSIC')
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

now = datetime.datetime.now().strftime("%B %d, %Y")
ncout.created  = "created based on reese@pik at " + now
ncout.data_origin = "downloaded from http://homepages.see.leeds.ac.uk/~earkhb/Basins_page.html"
ncout.regions = "regions follow the definition of Zwally et al. 2012"

print "close"
ncout.close()


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
plt.savefig('plots/ZwallyBasinsWithOcean.png')
plt.clf()
