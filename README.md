# GEDI Waveform Structural Complexity Index (WSCI)

This repository contains the models and code used to develop the GEDI WSCI data product (de Conto et al. 2024).

## Description 

The scripts provided here enable performing 2 tasks:
- (R) Calculate the 3D Canopy Entropy (CExyz) (Liu et al. 2022) from LiDAR point cloud samples.
- (python) Apply the WSCI models to a set of RH metrics.  

All software versions listed in this repository were used for developing and testing the scripts.

#### R code

The R packages necessary to run the `als_ce_xyz.R` tool are listed in the `r_session.txt` file. 

This code was originally used to extract CExyz measures from Airborle laser Scanning (ALS) data matched to GEDI footprints - same location and size (25m diameter plots).

##### Usage

```terminal 

Rscript als_ce_xyz.R
    [-[-help|h]]                        print help
    [-[-input|i] <character>]           path with point cloud LiDAR files (las/laz)
    [-[-output|o] <character>]          path to write structural complexity metrics - a directory when processing wall-to-wall <w2w> or a geospatial vector format otherwie (gpkg, shp, geojson etc.)
    [-[-plots_path|f] <character>]      (optional) input file with plot locations (e.g. GEDI footprints) in geospatial vector format (gpkg, shp, geojson etc.)
    [-[-n_plots|n] <integer>]           (optional) if no input plot locations are provided, how many random samples to draw? [default = 100]
    [-[-plot_size|s] <integer>]         (optional) plot diameter, in point cloud units [default = 25]
    [-[-gridded|g]]                     (optional) sample from a regular grid instead of randomly
    [-[-w2w|w]]                         (optional) wall-to-wall processing - calculate metrics for all pixels at <plot_size> resolution
    [-[-las_epsg|e] <integer>]          (optional) EPSG code of input LiDAR point clouds, if not encoded in the las/laz file
    [-[-voxel|v] <double>]              (optional) voxel size flter to apply to LiDAR files before processing 
    [-[-cores|c] <integer>]             (optional) number of CPU cores to use
```

##### Output
|    | column   | dtype   | description |
|---:|:---------|:--------|:------------|
|  0 | XY       | float   | horizontal entropy component |
|  1 | XZ       | float   | vertical entropy component |
|  2 | YZ       | float   | vertical entropy component |
|  3 | XYZ      | float   | 3D canopy entropy|
|  4 | p        | float   | p value from Mann-Kendall trend test|
|  5 | vox      | float   | voxel size |

When processing wall-to-wall (`w2w`) generates raster files with bands corresponding to the metrics in `column`.

#### Python code

The python libraries necessary to run the WSCI models are listed in the `python_requirements.txt` file.

This code was originally used to generate WSCI estimates from GEDI RH metrics extracted from the L2A product (Dubayah et al. 2020) at different Plant Functional Types (PFTs).

##### Usage

```terminal 
gedi_wsci.py [-h] -i INPUT -o OUTPUT [-p PFT] [-r RH_COLS] [-x INDEX [INDEX ...]]

Apply Waveform Structural Complexity (WSCI) model to input RH metrics

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        input RH metrics [cm] as columns (from 0 to 100)
  -o OUTPUT, --output OUTPUT
                        output file
  -p PFT, --pft PFT     PFT code or column name
  -r RH_COLS, --rh-cols RH_COLS
                        string pattern to select RH columns
  -x INDEX [INDEX ...], --index INDEX [INDEX ...]
                        columns to copy from input file into the output

```

The `PFT` codes used to call different WSCI models are listed below:

| PFT | short name  | long name                 |
|---: |:--------|:--------------------------|
|  1  | ENT     | Evergreen Needleaf Trees  |
|  2  | EBT     | Evergreen Broadleaf Trees |
|  4  | DBT     | Deciduous Broadleaf Trees |
|  5  | G       | Grassland                 |
|  6  | S       | Shrubland                 |
| 11  | W       | Woodland                  |
| -1  | -       | Other                     |

##### Output

|    | column   | dtype   | description |
|---:|:---------|:--------|:------------|
|  0 | wsci       | float   | Waveform Structural Complexity Index estimate |
|  1 | wsci_pi       | float   | WSCI prediction interval at 95% probability |
|  2 | wsci_xy       | float   | Estimated horizontal complexity component | 
|  3 | wsci_xy_pi      | float   | Horizontal component prediction interval at 95% probability |
|  4 | wsci_z        | float   | Estimated vertical complexity |
|  5 | wsci_z_pi      | float   | Vertical component prediction interval at 95% probability |


### References

De Conto, T., Armston, J. & Dubayah, R. O. Global Ecosystem Dynamics Investigation (GEDI)GEDI L4C Footprint Level Waveform Structural Complexity Index, Version 2. 0 MB Preprint at https://doi.org/10.3334/ORNLDAAC/2338 (2024).

Dubayah, R. et al. GEDI L2A Elevation and Height Metrics Data Global Footprint Level V002. NASA EOSDIS Land Processes Distributed Active Archive Center https://doi.org/10.5067/GEDI/GEDI02_A.002 (2021).

Liu, X. et al. A novel entropy-based method to quantify forest canopy structural complexity from multiplatform lidar point clouds. Remote Sensing of Environment 282, 113280 (2022).
