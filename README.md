# GEDI Waveform Structural Complexity Index (WSCI)

This repository contains the models and code used to develop the GEDI WSCI data product (de Conto et al. 2024).

## Description 

The scripts provided here enable performing 2 tasks:
- (R) Calculate the 3D Canopy Entropy (CExyz) (Liu et al. 2023) from LiDAR point cloud samples.
- (python) Apply the WSCI models to a set of RH metrics.  

All software versions listed in this repository were used for developing and testing the scripts.

#### R code

The R packages necessary to run the `als_ce_xyz.R` tool are listed in the `r_session.txt` file. 

This code was originally used to extract CExyz measures from Airborle laser Scanning (ALS) data matched to GEDI footprints - same location and size (25m diameter plots).

##### Usage

```terminal 
Rscript als_ce_xyz.R -h

Usage: 
    als_ce_xyz.R 
    [-[-help|h]]                        print help
    [-[-input|i] <character>]           path with point cloud LiDAR files (las/laz)
    [-[-output|o] <character>]          path to write structural complexity metrics in geospatial vector format (gpkg, shp, geojson etc.)
    [-[-plots_path|f] <character>]      (optional) input file with plot locations (e.g. GEDI footprints) in geospatial vector format (gpkg, shp, geojson etc.)
    [-[-n_plots|n] <integer>]           (optional) if no input plot locations are provided, how many random samples to draw? [default = 100]
    [-[-plot_size|s] <integer>]         (optional) plot diameter, in point cloud units [default = 25]
    [-[-gridded|g]]                     (optional) sample from a regular grid instead of randomly
    [-[-las_epsg|e] <integer>]          (optional) EPSG code of input LiDAR point clouds, if not encoded in the las/laz file
```


#### Python code

The python libraries necessary to run the WSCI models are listed in the `python_requirements.txt` file.

This code was originally used to generate WSCI estiamtes from GEDI RH metrics extracted from teh L2A product (Dubayah et al. 2020) at different Plant Functional Types (PFTs).

##### Usage

```terminal 
python gedi_wsci.py -h

usage: gedi_wsci.py [-h] -i INPUT -o OUTPUT [-p PFT] [-r RH_COLS] [-x INDEX [INDEX ...]]

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

### References

