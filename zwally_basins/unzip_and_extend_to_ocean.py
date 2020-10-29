"""
matthias.mengel@pik, torsten.albrecht@pik, ronja.reese@pik, ricarda.winkelmann@pik
Download Zwally basin data and save to (5km) netcdf file.
"""

import os, glob
import numpy as np
import numpy.ma as ma
import sys, csv, datetime
import netCDF4 as nc
import datetime, math

## this hack is needed to import config.py from the project root
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if project_root not in sys.path: sys.path.append(project_root)
import config as cf; reload(cf)
import pism_input.pism_input as pi; reload(pi)

dataset="zwally_basins"

"""
Parameters
"""
getData=1
mergeBasins=1
findLongitudes=1
adjustAP=1
fillInBasinMask=1


### Zwally basins   ##########################################################
# Documentation of the data: http://homepages.see.leeds.ac.uk/~earkhb/Basins_page.html
basins_link="http://homepages.see.leeds.ac.uk/~earkhb/ais_basins_imbie_ascii.zip"

## such file definitions should go to config.py, so that other functions can access them.
basins_data_path = os.path.join(cf.output_data_path, dataset)
basins_ascii_path = os.path.join(basins_data_path, 'basins_ascii')
infile = os.path.join(basins_ascii_path, 'ais_basins_imbie.z12.txt')
namefile = os.path.join(basins_ascii_path, 'ais_basins_imbie.z12.names.txt')
maskfile = os.path.join(cf.output_data_path, "albmap/Antarctica_5km_dev1.0.nc")

#ncout_name = os.path.join(basins_data_path, 'basins_data/basins_zwally_5km_input.nc')
ncout_raw = os.path.join(basins_data_path, 'basins_zwally_5km_drainage.nc')
ncout_shelves = os.path.join(basins_data_path, 'basins_zwally_5km_shelves.nc')
ncout_merged = os.path.join(basins_data_path, 'basins_zwally_5km_merged.nc')
ncout_ocean = os.path.join(basins_data_path, 'basins_zwally_5km_input.nc')


# if data is not yet extracted in basins_ascii
if not os.path.exists(basins_ascii_path):
  #print "Downloading Zwally basin ascii data."
  print "Unzipping Zwally basin ascii data."
  os.system("mkdir -p " + basins_ascii_path)
  #os.system("wget -N " + basins_link + " -P " + basins_ascii_path)
  os.system("unzip ais_basins_imbie_ascii.zip -d "+basins_ascii_path)
  #os.system("cd "+basins_ascii_path+" && unzip ais_basins_imbie_ascii.zip")

# if os.path.isfile(ncout_name):
#   print ncout_name, "already written, do nothing."


#########################################################################################

# if folder basins is not yet created
if not os.path.exists(os.path.join(basins_data_path,"basins_data")):
  os.system("mkdir -v " + basins_data_path)

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
  #print basinnum, basinname
  basinnames.append(basinname)
  basinnums.append(np.int(basinnum))

# add ocean cells to fix resolution
data2 = np.zeros([1200,1200])
for xi in range(len(x)):
  for yi in range(len(y)):
    data2[xi,yi] = data[xi,yi]


print "\nCreate netcdf file " + ncout_raw

ncout = nc.Dataset(ncout_raw, 'w', format='NETCDF4_CLASSIC')
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
  #print ncv
  ncout.sync()

#for row in namereader:
  #ncv[basinname][:] = data == basinnum

now = datetime.datetime.now().strftime("%B %d, %Y")
ncout.created  = "created based on matthias.mengel@pik by reese@pik-potsdam.de at " + now
ncout.data_origin = "downloaded from http://homepages.see.leeds.ac.uk/~earkhb/Basins_page.html"
ncout.regions = "regions follow the definition of 'Ice Flow of the Antarctic Ice Sheet', Rignot et al. 2011, Science"

ncout.close()

##############################################################################################################

print '\nExtend basins to ice shelves... '


infile = nc.Dataset(maskfile, 'r')
mask = ma.masked_array(infile.variables['mask'][0,:,:])
lon = ma.masked_array(infile.variables['lon'][0,:,:])
infile.close()
Mx=np.shape(mask)[0]
My=np.shape(mask)[1]
mask_floating = 2
mask_grounded = 1

#print 'Mx = ' + str(Mx)
#print 'My = ' + str(My)


infile = nc.Dataset(ncout_raw, 'r')
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

      loopcount = loopcount+1
      print '...loop ' + str(loopcount)
      if (loopcount == 2):
        done = 1


# # TODO: remove plotting from data processing scripts.
# print '\nPlotting... '

# import matplotlib.pyplot as plt
# import matplotlib.cm as cm
# import pylab as p

# from matplotlib.patches import Rectangle
# plt.rcParams['font.size']=30
# plt.rcParams['legend.fontsize']=30
# plt.rcParams['figure.figsize'] = 15, 15
# plt.close('all')


# fig=plt.figure()
# ax = plt.subplot(1,1,1)
# plt.contourf(basins)

# plt.contour(mask,mask_grounded,colors='gray',linewidths=2) # GROUNDING LINE
# plt.contour(mask,mask_floating,colors='gray',linewidths=2) # CALVING FRONT
# #cbar = plt.colorbar(ticks=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19])

# ax.set_xticks([])
# ax.set_yticks([])
# p.axis('equal')

# #plt.show()
# if not os.path.exists(os.path.join(basins_data_path,"plots")):
#   os.system("mkdir "+os.path.join(basins_data_path,"plots"))
# plt.savefig(os.path.join(basins_data_path,"plots/ZwallyBasins.png"))
# plt.clf()


"""
Create nc-file
"""
#print "create netcdf file\n" + ncout_shelves
ncout = nc.Dataset(ncout_shelves, 'w', format='NETCDF4_CLASSIC')
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

#print "close"
ncout.close()



##############################################################################################################

print '\nMerge basins... '


#infile = nc.Dataset(maskfile, 'r')
#mask = ma.masked_array(infile.variables['mask'][0,:,:])
#lon = ma.masked_array(infile.variables['lon'][0,:,:])
#infile.close()
#Mx=np.shape(mask)[0]
#My=np.shape(mask)[1]
#mask_floating = 2
#mask_grounded = 1

#print 'Mx = ' + str(Mx)
#print 'My = ' + str(My)


"""
Merge basins:
"""
if (mergeBasins==1):
    print '\nMerge basins... '
    infile = nc.Dataset(ncout_shelves, 'r')
    basins = ma.masked_array(infile.variables['basins'][:,:])
    infile.close()

    basins_dict = {0.0:0.0, 1.0:1.0, 2.0:1.0 ,3.0:1.0, 4.0:2.0, 5.0:3.0, 6.0:4.0, 7.0:5.0,
                  8.0:5.0, 9.0:6.0, 10.0:6.0, 11.0:6.0, 12.0:7.0, 13.0:8.0, 14.0:9.0,
                  15.0:10.0, 16.0:11.0, 17.0:12.0, 18.0:12.0, 19.0:12.0, 20.0:13.0,
                  21.0:14.0, 22.0:14.0, 23.0:15.0, 24.0:16.0, 25.0:17.0, 26.0:18.0, 27.0:19.0}


    numberOfBasins = int(basins.max())

    for i in range(Mx):
      for j in range(My):
        #print i,j,basins[i,j],basins_dict[basins[i,j]]
        basins[i,j] = basins_dict[basins[i,j]]

    """
    Create nc-file
    """
    print "\nCreate netcdf file\n" + ncout_merged
    ncout = nc.Dataset(ncout_merged, 'w', format='NETCDF4_CLASSIC')
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



##############################################################################################################


'''
Adjust Basin 18/19 (originally 26/27)
'''
#print '\nReading basins... '
infile = nc.Dataset(ncout_merged, 'r')
basins = ma.masked_array(infile.variables['basins'][:,:])
infile.close()

if (adjustAP):
  print '\nAdjusting Antarctic Peninsula Basin which was wrong set in add_shelves_to_basins.py'

  for i in range(1, Mx-1):
    for j in range(1,My-1):
      if(basins[i,j]==19.0 and j>=130 and j<=166 and i>=763 and i<=800):
        #if ((i-763.0)/(796-763)> (j-134.0)/(166-134)):
        basins[i,j]=18.0


##############################################################################################################

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
    np.save(os.path.join(basins_data_path,"longitudeVec.npy"), longitudeVec)


##############################################################################################################


if (fillInBasinMask) :
  print '\nFilling basin mask... '
  longitudeVec = np.load(os.path.join(basins_data_path,"longitudeVec.npy"))
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


##############################################################################################################

"""
Create nc-file
"""
print "\nCreate netcdf input file for pism \n" + ncout_ocean
ncout = nc.Dataset(ncout_ocean, 'w', format='NETCDF4_CLASSIC')
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
#ncout.created  = "created by reese@pik at " + now
#ncout.data_origin = "downloaded from http://homepages.see.leeds.ac.uk/~earkhb/Basins_page.html"
ncout.Descricption = "Antarctic drainage basins mapped by NASA and modified." ;
ncout.Reference = "Zwally, H. Jay, Mario B. Giovinetto, Matthew A. Beckley, and Jack L. Saba, 2012, Antarctic and Greenland Drainage  Systems, GSFC Cryospheric Sciences Laboratory, at http://icesat4.gsfc.nasa.gov/cryo_data/ant_grn_drainage_systems.php." ;
ncout.regions = "regions follow the definition of Zwally et al. 2012"
ncout.proj4 = "+lon_0=0.0 +ellps=WGS84 +datum=WGS84 +lat_ts=-71.0 +proj=stere +x_0=0.0 +units=m +y_0=0.0 +lat_0=-90.0"
ncout.comment  = cf.authors+" created netcdf basins file at " + now

ncout.close()

# prepare the input file for cdo remapping
# this step takes a while for high resolution data (i.e. 1km)
pi.prepare_ncfile_for_cdo(ncout_ocean)


# print '\nPlotting... '
# fig=plt.figure()
# ax = plt.subplot(1,1,1)
# plt.contourf(basins)

# plt.contour(mask,mask_grounded,colors='gray',linewidths=2) # GROUNDING LINE
# plt.contour(mask,mask_floating,colors='gray',linewidths=2) # CALVING FRONT
# #cbar = plt.colorbar(ticks=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19])

# ax.set_xticks([])
# ax.set_yticks([])
# p.axis('equal')

# #plt.show()
# plt.savefig(os.path.join(basins_data_path,"plots/ZwallyBasinsWithOcean.png"))
# plt.clf()

# print "\nDone"

