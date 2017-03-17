## Preprocessing for PISM


This repository contains scripts to preprocess input data for PISM.
Ideally, the scripts provide all steps from downloading the data to
making the PISM-ready, i.e. no other steps in preprocessing should be needed
to use the data with PISM.

Initial conditions for PISM include bed topography, ice thickness and ice
ice temperature distribution.

Boundary conditions include ice surface temperature and ice surface mass balance or
atmospheric temperature at the snow surface and mass balance at the snow surface.
Ocean melt rates and ice shelf basal temperatures need to be provided from the
ocean. PISM includes simple ocean models that can provide such data from open
ocean properties (i.e. the SIMPEL ocean model).


### Structure

Each dataset should get a separate folder. Such folder should contain scripts
that are named by what they do. The `merging` folder contains scripts that
merge several datasets to one. The `regridding` folder provides
more general regridding options that are not specific to a certain dataset.
Here is the graph:

```
pism_input
+-- bedmap2
|   +-- download_and_extract_to_nc.py
|   +-- regrid.py
|   +-- rename_to_pism_variables.py
+-- albmap
|   +-- download.py
|   +-- regrid.py
|   +-- rename_to_pism_variables.py
+-- merging
+-- regridding
```
