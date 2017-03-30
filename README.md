## Preprocessing for PISM


This is a collection of preprocessing routines for input data to be used in PISM
for Antarctic ice sheet simulations.
The preprocessing provides the steps from downloading the data to
making the PISM-ready, i.e. no other steps in preprocessing should be needed
to use the data with PISM.

Initial conditions for PISM include bed topography, ice thickness and
ice temperature distribution.

Boundary conditions include ice surface temperature and ice surface mass balance or
atmospheric temperature at the snow surface and mass balance at the snow surface.
Ocean melt rates and ice shelf basal temperatures need to be provided from the
ocean. PISM includes simple ocean models that can provide such data from open
ocean properties (i.e. the SIMPEL ocean model).


### Structure

Scripts for dataset-specific preprocessing are in the folder with the name
of the dataset.
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
+-- albmap
|   +-- download_and_rename_variables.py
|   +-- remap.py
+-- racmo
|   +-- preprocess.sh (TODO: flipud)
|   +-- remap_racmo.py
+-- basins
|   +-- create_basins_netcdf.py
|   +-- add_shelves_to_basins.py
|   +-- add_ocean_to_basins.py
|   +-- remap_basins.py (with local cdo rempnn)
+-- schmidtko
|   +-- create_NetCDF.py
|   +-- calculate_potential_temps.py
|   +-- remap_schmidtko.py
|   +-- compute_basin_means.py
+-- accum
|   +-- download_and_extract_to_nc.py
|   +-- remap.py
+-- initmip
|   +-- preprocess.sh
|   +-- remap.py (asmb)
|   +-- postprocess.sh (incl. cdo remapnn abmb)
+-- merging
+-- pism_input
    +-- pism_input.py
```
