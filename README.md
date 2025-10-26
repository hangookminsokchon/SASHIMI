# SASHIMI
SASHIMI (Spatial Analysis for Segmented Histopathology Images using Machine Intelligence) is a tool for capturing spatial cell-cell interactions in histopathology image using spatial summary statistics.

## Description
**SASHIMI** is a computational framework designed to extract quantitative, statistically robust summaries of cell-cell interactions within the tumor tissue architecture. Leveraging spatial statistics and point pattern analysis, SASHIMI translates raw cell coordinate data from histopathological images into meaningful scalar and functional descriptors that characterize tissue organization and cellular architecture.

### Key Features

- **Spatial Summary Statistics**: Spatial summary statistics are functions that analyze the spatial characteristics of points based on their locations and relationships. These functions characterize properties such as clustering, regularity, and inter-type spatial relationships.

- **Areal Data Indices**: This representation enables quantification of spatial structure through statistics that assess cell-type aggregation, dispersion, and compositional similarity across tissue regions.

- **Topological Features**: Topology, a branch of mathematics concerned with the qualitative properties of geometric structures, forms the basis of topological data analysis (TDA). It characterizes spatial patterns through geometric features such as the number of connected components and loops. 


## System Requirements
This framework requires the use of following software dependencies. 

R Environment (v4.5.1+)
```{r}
library(spatstat)  #(version >= 3.3)
library(dispRity)  #(version >= 1.9)
library(spdep)     #(version >= 1.3)
library(dplyr)     #(version >= 1.1.4)
```

Python Environment (v3.11+)
```{python}
import pandas as pd                              #(version >= 2.3)
import numpy as np                               #(version >= 2.2.0)

import matplotlib.pyplot as plt                  #(version >= (3.10)
 from matplotlib.colors import ListedColormap

import scipy import ndimage                      #(version >= (1.16)
 from scipy.ndimage import distance_transform_bf
 from scipy.stats import gaussian_kde

import gudhi                                     #(version >= (3.11)

```
## Web page
### ** To be added ** 


## Input/Output Specifications
### Standardized Input Format

All feature modules accept a unified CSV format with the following structure:

| Column | Type | Description |
|--------|------|-------------|
| `x` | float | X-coordinate of cell centroid |
| `y` | float | Y-coordinate of cell centroid |
| `type` | string | Cell classification (`immune`, `stromal`, `tumor`, `other`) |

**File Requirements:**
- Format: CSV (comma-separated values)
- Maximum size: 4MB

### Output Specifications
All types of computed features are downloadable in .csv format, directly from the web.

#### Areal Data  
- **Format**: 1 × m DataFrame
- **Content**: Autocorrelation/Similarity indicies in scalar values.
- **Example metrics**: Moran's I, Geary's C, Cosine similarity, etc...
  
#### Functional Data  
- **Format**: 500 × 3 DataFrame
- **Content**: Functional data representing spatial functions, computed at 500 bin width points.
- **Example metrics**: K-function, Pair correlation function, G-fucntion, etc...

#### Topological Data  
- **Format**: 1 × m DataFrame
- **Content**: Summary statistics of persistence diagram
- **Example metrics**: min/max/std of Betti 0, 1 numbers

## Repository Structure

```
SASHIMI/
├── src/                        # Source code and computation modules
│   ├── areal_features.R           # Computes a suite of Areal Data Indices
│   ├── functional_features.R      # Computes a suite of Spatial Summary Statistics
│   ├── helperfunctions.R          # Helper functions for image batch normalization & cell_type regularization
│   ├── topological_features.py    # Computes Topological Features
│   ├── workflow_example.R                 # Data pipeline example of Area Data Indices/Spatial Summary Statistics
│   └── workflow_example_topological.ipynb # Data pipeline example of Topological Features
|
├── data/                       # Example input datasets
|   ├── example_imageA.png        # .png image of example A point pattern data
|   ├── example_imageB.png        # .png image of example B point pattern data
│   ├── example_point_patternA.csv  #  .csv input data of example A     
│   └── example_point_patternB.csv  #  .csv input data of example B      
├── example/                    # Example outputs
|
|   ├── output_example_areal.csv     # Example output of Areal Data Indices, computed from example A
|   ├── output_exampleA_functional_K # Example output of Spatial Summary Statistic(K-function), computed from example A
|   ├── output_exampleB_topological1.csv # Example output of Topological Feature, computed from example B
│   └── output_exampleB_topological2.csv # Example output of Topological Feature, computed from example B, by dimension
|
├── README.MD
├── LICENSE
└── DESCRIPTION.txt
```


## DEVLOGS

## Comments for Developers

The SASHIMI framework consists of three feature modules:

- **Spatial Summary Statistics** *(currently available)*
- **Areal/Autocorrelation/Similarity Indices** *(currently available)*
- **Topological** *(currently available)*

### Updates
June 3, 2025: 

Initial release with Areal feature module implementation. Functional and Topological feature modules scheduled for mid-to-late July release.

⸻

July 15, 2025:

- Areal features module completed (additional features may be added in future releases)
- Functional features module implemented
- File naming updates for improved consistency:

 *areal_data_features → areal_feature*
 
 *functional_data_features → functional_features*
 
 *compute_features → feature_computation* 
 
⸻

Aug 29, 2025
New Feature: Added topological_features.py.
- Dependencies: Updated README.md and DESCRIPTION.txt to include required Python packages and their versions.
- Computation: For efficiency, the exampleB image is now used for computing topological features (different from exampleA, which is still used for areal and functional features).
- Current Scope: Only scalar topological features are supported at present. A functional version will be added in a future release.

Repository Structure
The repository now consists of three main folders:

/src – source code and example workflow scripts

/data – example data, including point pattern images and raw .csv files

/example – output data generated from /data inputs

Reminders:
Column name specifications are standardized across all features.
- Example A: x, y, Z_cell
- Example B: x, y, class
- Cell type names are also fixed.
  
A more flexible implementation for column and cell-type names is planned, including additional helper functions.

⸻

Sep 14, 2025
Bug Fixes:
- Fixed KeyError in read_img (line 59, topological_features.py).
- Fixed FutureError in compute_cubical_complex_pair (lines 148–149, topological_features.py).

Input Standardization:
- All feature types (functional, areal, topological) now use a unified input format: three columns (x, y, type) where type = {'immune', 'stromal', 'tumor', 'other'}.
- The column name does not strictly have to be "type"; custom naming is supported.
