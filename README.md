## Collection of pre- and postprocessing scripts for PISM input.

### Input

Initial conditions for PISM include bed topography, ice thickness,
ice temperature distribution.

Boundary conditions include ice surface temperature and ice surface mass balance or
atmospheric temperature at the snow surface and mass balance at the snow surface.
Ocean melt rates and ice shelf basal temperatures need to be provided from the
ocean. PISM includes simple ocean models that can provide such data from open
ocean properties (i.e. the SIMPEL ocean model).

##
Structure

```
pism_input
+-- download
+-- preprocess
+-- regridding
+-- merging
+-- postprocessing
|   +-- example
|   +--
+-- example
```
