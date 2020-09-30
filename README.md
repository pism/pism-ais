Preprocessing for PISM
=============


This is a collection of preprocessing routines for input data to be used in PISM
for Antarctic ice sheet simulations.
The preprocessing provides the steps from downloading the data to
making the data PISM-ready, i.e. no other steps in preprocessing should be needed
to use the data with PISM.

Initial conditions for PISM include bed topography, ice thickness and
ice temperature distribution.

Boundary conditions include ice surface temperature and ice surface mass balance or
atmospheric temperature at the snow surface and mass balance at the snow surface.
PISM needs ocean melt rates and ice shelf basal temperatures as boundary conditions
at the ice shelf base.
PISM includes simple ocean models that can provide such data from open
ocean properties (i.e. the PICO ocean model).


## Structure

Scripts for dataset-specific preprocessing are in the folder with the name
of the dataset. They should be run in the order listed below.
The `grids` folder holds scripts for the creation of target grids to which the
input data is regridded. The `pism_input` folder holds general program code.
The `merging` folder holds code to merge single datasets into
files that hold all the necessary initial and boundary conditions to
drive PISM.
Here is the graph:


```
pism_input
+-- grids
|   +-- create_cdo_targetgrids.py
+-- bedmap2
|   +-- download_and_extract_to_nc.py
|   +-- remap.py
+-- bedmachine
|   +-- remap.py
+-- albmap
|   +-- download_and_rename_variables.py
|   +-- remap.py
+-- racmo_ice2sea
|   +-- merge_and_rename_variables.py
|   +-- remap.py
|   +-- time_mean.py
+-- racmo_wessem
|   +-- merge_and_rename_variables.py
|   +-- remap.py
|   +-- time_mean.py
+-- racmo_hadcm3
|   +-- merge_and_rename_variables_c20/a1b.py
+-- racmo_cesm2
|   +-- merge_and_rename_variables_hist/ssp585.py
|   +-- remap.py
|   +-- time_mean.py
|   +-- cdo_remap_hist/ssp858.sh
|   +-- cdo_remap_yearly.sh (example for year-wise processing)
+-- zwally_basins
|   +-- unzip_and_extend_to_ocean.py
|   +-- remap.py (with local cdo rempnn)
+-- icesat4_basins
|   +-- download_and_write_to_netcdf.py
|   +-- remap.py (FIXME: cdo rempnn)
+-- schmidtko
|   +-- download_and_write_to_netcdf.py
|   +-- calculate_potential_temps.py
|   +-- remap.py
|   +-- compute_basin_means.py
|   +-- fill_schmidtko_initmip8km_means_to_basins_on_different_resolution.py 
+-- accum
|   +-- download_and_extract_to_nc.py
|   +-- remap.py
+-- bheatflx_martos
|   +-- download_and_write_to_netcdc.py
|   +-- remap.py
+-- bheatflx_an
|   +-- download_and_extract.py
|   +-- remap.py
+-- bheatflx_purucker
|   +-- download_and_interpolate_to_netcdc.py
+-- bheatflx_shen
|   +-- interpolate_to_netcdc.py
+-- vel_rignot
|   +-- preprocess_netcdf.py
|   +-- remap.py
+-- vel_mouginot17_annual
|   +-- preprocess_netcdf.py
+-- vel_mouginot19
|   +-- preprocess_netcdf.py
|   +-- remap.py
+-- raised
|   +-- convert_to_netcdf.py
|   +-- remap.py
+-- sealevel
|   +--download_specmap_and_rename.py
|   +--convert_ice6g_to_netcdf.py
|   +--download_bintanja_and_convert_to_nc.py
|   +--convert_lambeck_to_netcdf.py
+-- wdc
|   +--convert_tempdata_to_netcdf.py
|   +--convert_accumdatato_netcdf.py
+-- edc
|   +--download_tempdata_and_convert_to_netcdf.py
+-- initmip
|   +-- preprocess.py (incl. create_anomalies.py)
|   +-- remap.py (if not on 32,16,8,1km grid)
|   +-- postprocessing (make data ISMIP6 compatible)
+-- larmip
|   +-- preprocess.py (incl. create_anomalies.py)
|   +-- postprocessing (extract slvol anomalies)
+-- merge_dataset.py
+-- pism_input
    +-- pism_input.py (tools)
```


## Data sources



### Surface elevation

#### Albmap data

5km resolution, NetCDF files of Antarctica gathered from various data sources and interpolated, when necessary, onto the same grid using polar stereographic projection

Documentation: <http://websrv.cs.umt.edu/isis/index.php/Present_Day_Antarctica>

Citation:  
- Bamber, J. L., Gomez-Dans, J. L., & Griggs, J. A. (2009). A new 1 km digital elevation model of the Antarctic derived from combined satellite radar and laser data-Part 1: Data and methods. The Cryosphere, 3(1), 101.

- Griggs, J. A., & Bamber, J. L. (2009). A new 1 km digital elevation model of Antarctica derived from combined radar and laser data-Part 2: Validation and error estimates. The Cryosphere, 3(1), 113.

#### Bedmap2 data

Dataset describing surface elevation, ice-thickness and the sea ﬂoor and subglacial bed elevation of the Antarctic south of 60◦ S. Dataset was derived using data from a variety of sources, including many substantial surveys completed since the original Bedmap compilation (Bedmap1) in 2001. The Bedmap2 ice thickness grid is made from 25 million measurements, over two orders of magnitude more than were used in Bedmap1.

Documentation: <https://www.bas.ac.uk/project/bedmap-2/>

Citation: 
- Fretwell, P., Pritchard, H. D., Vaughan, D. G., Bamber, J. L., Barrand, N. E., Bell, R., Bianchi, C., Bingham, R. G., Blankenship, D. D., Casassa, G., Catania, G., Callens, D., Conway, H., Cook, A. J., Corr, H. F. J., Damaske, D., Damm, V., Ferraccioli, F., Forsberg, R., Fujita, S., Gim, Y., Gogineni, P., Griggs, J. A., Hindmarsh, R. C. A., Holmlund, P., Holt, J. W., Jacobel, R. W., Jenkins, A., Jokat, W., Jordan, T., King, E. C., Kohler, J., Krabill, W., Riger-Kusk, M., Langley, K. A., Leitchenkov, G., Leuschen, C., Luyendyk, B. P., Matsuoka, K., Mouginot, J., Nitsche, F. O., Nogi, Y., Nost, O. A., Popov, S. V., Rignot, E., Rippin, D. M., Rivera, A., Roberts, J., Ross, N., Siegert, M. J., Smith, A. M., Steinhage, D., Studinger, M., Sun, B., Tinto, B. K., Welch, B. C., Wilson, D., Young, D. A., Xiangbin, C., and Zirizzotti, A.: Bedmap2: improved ice bed, surface and thickness datasets for Antarctica, The Cryosphere, 7, 375-393, doi:10.5194/tc-7-375-2013, 2013.

### Bedrock topography

#### Bedmap2 data

Dataset describing surface elevation, ice-thickness and the sea ﬂoor and subglacial bed elevation of the Antarctic south of 60◦ S. Dataset was derived using data from a variety of sources, including many substantial surveys completed since the original Bedmap compilation (Bedmap1) in 2001. The Bedmap2 ice thickness grid is made from 25 million measurements, over two orders of magnitude more than were used in Bedmap1.

Documentation: <https://www.bas.ac.uk/project/bedmap-2/>

Citation: 
- Fretwell, P., Pritchard, H. D., Vaughan, D. G., Bamber, J. L., Barrand, N. E., Bell, R., Bianchi, C., Bingham, R. G., Blankenship, D. D., Casassa, G., Catania, G., Callens, D., Conway, H., Cook, A. J., Corr, H. F. J., Damaske, D., Damm, V., Ferraccioli, F., Forsberg, R., Fujita, S., Gim, Y., Gogineni, P., Griggs, J. A., Hindmarsh, R. C. A., Holmlund, P., Holt, J. W., Jacobel, R. W., Jenkins, A., Jokat, W., Jordan, T., King, E. C., Kohler, J., Krabill, W., Riger-Kusk, M., Langley, K. A., Leitchenkov, G., Leuschen, C., Luyendyk, B. P., Matsuoka, K., Mouginot, J., Nitsche, F. O., Nogi, Y., Nost, O. A., Popov, S. V., Rignot, E., Rippin, D. M., Rivera, A., Roberts, J., Ross, N., Siegert, M. J., Smith, A. M., Steinhage, D., Studinger, M., Sun, B., Tinto, B. K., Welch, B. C., Wilson, D., Young, D. A., Xiangbin, C., and Zirizzotti, A.: Bedmap2: improved ice bed, surface and thickness datasets for Antarctica, The Cryosphere, 7, 375-393, doi:10.5194/tc-7-375-2013, 2013.


### BedMachine Antarctica

Dataset based on mass conservation, version = "24-Jan-2019 (v1.27), retrieved from Mathieu Morlighem, personal communication (not published yet, probably mid of 2019)

Website: <https://sites.uci.edu/morlighem/bedmachine-antarctica/>

Citation:
- M. Morlighem, E. Rignot, J. Mouginot, H. Seroussi and E. Larour, Deeply incised submarine glacial valleys beneath the Greenland Ice Sheet, Nat. Geosci., 7, 418-422, 2014, doi:10.1038/ngeo2167



### Accumulation

#### Arthern accumulation

Antarctic surface accumulation map on the same grid, with the same projection and in the same file formats as bedmap2.

Documentation: <https://secure.antarctica.ac.uk/data/bedmap2/resources/Arthern_accumulation/Arthern_accumulation_bedmap2_grid_readme.rtf>

Citation: 
- Arthern, R. J., D. P. Winebrenner, and D. G. Vaughan (2006), Antarctic snow accumulation mapped using polarization of 4.3-cm wavelength microwave emission, J. Geophys. Res., 111, D06107, doi:10.1029/2004JD005667.


#### Albmap data

5km resolution, NetCDF files of Antarctica gathered from various data sources and interpolated, when necessary, onto the same grid using polar stereographic projection

Documentation: <http://websrv.cs.umt.edu/isis/index.php/Present_Day_Antarctica>

Citation: 
- Vaughan, D. G., Bamber, J. L., Giovinetto, M., Russell, J., & Cooper, A. P. R. (1999). Reassessment of net surface mass balance in Antarctica. Journal of climate, 12(4), 933-946.



### Surface mass balance

#### Racmo data (v.2.1)

HadCM3_c20_I2S_precip_Y.nc contains precipitation in mm/yr for 1980-1999 over Antarctica from RACMO2 run forced with HadCM3 data. HadCM3_c20_I2S_t2m_Y.nc contains temperature in Kelvin for 1980-1999 over Antarctica from RACMO2 run forced with HadCM3 data. (ANT55/HAD)

Citation: 
- Ligtenberg, S. R. M., van de Berg, W. J., van den Broeke, M. R., Rae, J. G. L. & van Meijgaard, E. (2013). Future surface mass balance of the Antarctic Ice Sheet and its influence on sea level change, simulated by a regional atmospheric climate model. Clim. Dynam. 41, 867–884

###### Note:
Data was given for use in this paper:
- Frieler Katja; Clark Peter U.; He Feng; Buizert Christo; Reese Ronja; Ligtenberg Stefan R. M.; van den Broeke Michiel R.; Winkelmann Ricarda & Levermann Anders Consistent evidence of increasing Antarctic accumulation with warming. Nature Clim. Change, 2015, 5, 348-352



#### Racmo data (v.2.3p2)

The latest RACMO2.3p2 data (ANT27/2) forced by *ERA-Interim* provide yearly mean air temperature (t2m) and surface mass balance (smb) for the years 1979-2016. 


Citation:
- Van Wessem, Jan Melchior, Willem Jan Van De Berg, Brice PY Noël, Erik Van Meijgaard, Charles Amory, Gerit Birnbaum, Constantijn L. Jakobs et al. "Modelling the climate and surface mass balance of polar ice sheets using racmo2: Part 2: Antarctica (1979-2016)." Cryosphere 12, no. 4 (2018): 1479-1498.

download link: https://www.projects.science.uu.nl/iceclimate/publications/data/2018/vwessem2018_tc/RACMO_Yearly/

see also overview: https://www.projects.science.uu.nl/iceclimate/models/antarctica.php


RACMO2.3p2 data (ANT27/2) forced by *HadCM3* climate model provide yearly mean air temperature (tskin) and surface mass balance (smb) for the historical period 1980-1999 and the A1B projection period 2000-2200.

RACMO2.3p2 data (ANT27/2) forced by *CESM2* climate model provide monthly mean air temperature (tskin) and surface mass balance (smb) for the historical period 1950-2014 and the SSP5-8.5 projection period 2015-2100.

not yet published, contact J.M.vanWessem[at]uu.nl




### Mean annual temperature

#### Albmap data

5km resolution, NetCDF files of Antarctica gathered from various data sources and interpolated, when necessary, onto the same grid using polar stereographic projection

Documentation: <http://websrv.cs.umt.edu/isis/index.php/Present_Day_Antarctica>

Citation: 
- Comiso, J. C. (2000). Variability and trends in Antarctic surface temperatures from in situ and satellite infrared measurements. Journal of Climate, 13(10), 1674-1696.


### Geothermal flux

overview white paper:
- https://www.scar.org/scar-news/serce-news/scar-serce-ghf

#### Albmap data

5km resolution, NetCDF files of Antarctica gathered from various data sources and interpolated, when necessary, onto the same grid using polar stereographic projection

Documentation: <http://websrv.cs.umt.edu/isis/index.php/Present_Day_Antarctica>

Citation: 
- Shapiro, N. M., & Ritzwoller, M. H. (2004). Inferring surface heat flux distributions guided by a global seismic model: particular application to Antarctica. Earth and Planetary Science Letters, 223(1), 213-224.

- Maule, C. F., Purucker, M. E., Olsen, N., & Mosegaard, K. (2005). Heat flux anomalies in Antarctica revealed by satellite magnetic data. Science, 309(5733), 464-467.

#### Martos data

Documentation: https://doi.pangaea.de/10.1594/PANGAEA.882503

Citation: 
- Martos, Y. M., Catalán, M., Jordan, T. A., Golynsky, A., Golynsky, D., Eagles, G., & Vaughan, D. G. (2017). Heat flux distribution of Antarctica unveiled. Geophysical Research Letters, 44(22).

#### Purucker data

Documentation: https://core2.gsfc.nasa.gov/research/purucker/heatflux_updates.html

Citation: 
- Purucker, M. E.: Geothermal heat flux data set based on low res- olution observations collected by the CHAMP satellite between 2000 and 2010, and produced from the MF-6 model following the technique described in Fox Maule et al. (2005), 2013


#### An data

Documentation: http://www.seismolab.org/model/antarctica/lithosphere/#an1-hf

Citation:
- An, M., Wiens, D.A., Zhao, Y., Feng, M., Nyblade, A., Kanao, M., Li, Y., Maggi, A. and Lévêque, J.J., 2015. Temperature, lithosphere‐asthenosphere boundary, and heat flux beneath the Antarctic Plate inferred from seismic velocities. Journal of Geophysical Research: Solid Earth, 120(12), pp.8720-8742.


#### Shen data

Documentation: https://sites.google.com/view/weisen/research-products?authuser=0

Citation:
- Shen, W., Wiens, D., Lloyd, A. and Nyblade, A., 2020. A Geothermal heat flux map of Antarctica empirically constrained by seismic structure. Geophysical Research Letters, https://doi.org/10.1029/2020GL086955



### Southern Ocean

#### Schmidtko dataset.

Citation: 
- Schmidtko, Sunke, Heywood, Karen J., Thompson, Andrew F. and Aoki, Shigeru (2014) Multidecadal warming of Antarctic waters Science, 346 (6214). pp. 1227-1231.



### Drainage regions

#### Zwally dataset

Documentation: <http://homepages.see.leeds.ac.uk/~earkhb/Basins_page.html>

Citation: 
- Zwally, H. J., Giovinetto, M. B., Beckley, M. A., & Saba, J. L. (2012). Antarctic and Greenland drainage systems. GSFC Cryospheric Sciences Laboratory.


#### Icesat4 dataset

Documentation: https://icesat4.gsfc.nasa.gov/cryo_data/ant_grn_drainage_systems.php

Citation:
- Zwally, H. J., Giovinetto, M. B., Beckley, M. A., & Saba, J. L. (2012). Antarctic and Greenland drainage systems. GSFC Cryospheric Sciences Laboratory.



### Surface Velocity

#### Rignot data (2017)

Citation: 
- Rignot, E., J. Mouginot, and B. Scheuchl. 2017. MEaSUREs InSAR-Based Antarctica Ice Velocity Map, Version 2. [Indicate subset used]. Boulder, Colorado USA. NASA National Snow and Ice Data Center Distributed Active Archive Center. doi: http://dx.doi.org/10.5067/D7GK8F5J8M8R. [Date Accessed].

More infos: http://nsidc.org/data/nsidc-0484



#### Mouginot annual data (2017)

Citation: 
- Mouginot, J., E. Rignot, B. Scheuchl, and R. Millan. 2017. Comprehensive Annual Ice Sheet Velocity Mapping Using Landsat-8, Sentinel-1, and RADARSAT-2 Data, Remote Sensing. 9. Art. #364. https://doi.org/10.3390/rs9040364

More infos: https://nsidc.org/data/NSIDC-0720/



#### Mouginot data (2019)

Citation: 
- Mouginot, J., E. Rignot, and B. Scheuchl. 2019. MEaSUREs Phase-Based Antarctica Ice Velocity Map, Version 1. [Indicate subset used]. doi: https://doi.org/10.5067/PZ3NJ5RXRH10. [Date Accessed].

More infos: https://nsidc.org/data/NSIDC-0754/


### Grounding line

#### RAISED

Citation: 
- Bentley, Michael J., et al. "A community-based geological reconstruction of Antarctic Ice Sheet deglaciation since the Last Glacial Maximum." Quaternary Science Reviews 100 (2014): 1-9.


Vertices downloaded and converted to csv as publication [supplement](https://www.sciencedirect.com/science/article/pii/S0277379114002546#appsec1) from:

https://ars.els-cdn.com/content/image/1-s2.0-S0277379114002546-mmc1.xlsx



### Temperature forcing

#### EDC

Citation: Jouzel, J., Masson-Delmotte, V., Cattani, O., Dreyfus, G., Falourd, S., Hoffmann, G., ... & Fischer, H. (2007). Orbital and millennial Antarctic climate variability over the past 800,000 years. science, 317(5839), 793-796.

Anomaly over the last 800kyr provided as difference from the average of the last 1000 years, downloaded from 

ftp://ftp.ncdc.noaa.gov/pub/data/paleo/icecore/antarctica/epica_domec/edc3deuttemp2007.txt


#### WDC

Citation: Cuffey, K.M., G.D. Clow, E.J. Steig, C. Buizert, T.J. Fudge, M. Koutnik, E.D. Waddington, R.B. Alley, and J.P. Severinghaus (2016). Deglacial temperature history of West Antarctica. Proc. Natl. Acad. Sci. 113(50): 14249-14254. doi:10.1073/pnas.1609132113.

Download txt.data from http://www.usap-dc.org/view/dataset/600377 and calculate anomaly to year 0 BP over the last 67kyr.



### Accumulation forcing

#### WDC

Citation: Fudge, T. J., et al. "Variable relationship between accumulation and temperature in West Antarctica for the past 31,000 years." Geophysical Research Letters 43.8 (2016): 3795-3803.

Download xls file from http://www.usap-dc.org/view/dataset/601004 and convert sheet "50year" to csv.



### Sea level forcing

#### ice6g

Citation: Stuhne, G. R., and W. R. Peltier. "Reconciling the ICE‐6G_C reconstruction of glacial chronology with ice sheet dynamics: The cases of Greenland and Antarctica." Journal of Geophysical Research: Earth Surface 120.9 (2015): 1841-1865.

Anomaly over last 122kyr. Data not public (contact Dick Peltier).


#### specmap

Citation: Imbrie, J. D., and A. McIntyre. "SPECMAP time scale developed by Imbrie et al., 1984 based on normalized planktonic records (normalized O-18 vs. time, specmap. 017)." Earth System Science Data (2006).

Anomaly over last 405kyr, data extracted from [Albmap dataset](http://websrv.cs.umt.edu/isis/images/4/4d/Antarctica_5km_dev1.0.nc).


#### bintanja08

Citation: Bintanja, R., & Van de Wal, R. S. W. (2008). North American ice-sheet dynamics and the onset of 100,000-year glacial cycles. Nature, 454(7206), 869.

Anomaly data over last 3,000kyr downloaded from: ftp://ftp.ncdc.noaa.gov/pub/data/paleo/contributions_by_author/bintanja2008/bintanja2008.txt


#### lambeck14

Citation: Lambeck, K., Rouby, H., Purcell, A., Sun, Y., & Sambridge, M. (2014). Sea level and global ice volumes from the Last Glacial Maximum to the Holocene. Proceedings of the National Academy of Sciences, 111(43), 15296-15303.

Anomaly data over 35kyr extracted from Table 3 in pdf supplement http://www.pnas.org/content/suppl/2014/10/08/1411762111.DCSupplemental




### Model intercomparisons

#### InitMIP

Wiki: http://www.climate-cryosphere.org/wiki/index.php?title=InitMIP-Antarctica


#### Larmip

Wiki: https://www.pik-potsdam.de/research/earth-system-analysis/models/larmip


#### Abumip

Wiki: http://www.climate-cryosphere.org/wiki/index.php?title=ABUMIP-Antarctica


## License

pism-ais is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; version 2 of the License, or (at your option) any later version.

psim-ais is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You find a copy of the GNU General Public License along with pism-ais in the file LICENSE.txt;
