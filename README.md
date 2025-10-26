# SASHIMI: Spatial Analysis for Segmented Histopathology Images using Machine Intelligence
SASHIMI is a tool for capturing spatial cell-cell interactions in histopathology image using spatial summary statistics.

## Description
**SASHIMI** is a web-based framework developed on R and Python for the extraction, vi-
sualization, and computation of spatial features in AI-segmented histopathology images,
enabling real-time analysis. The framework serves as an exploratory and feature-extraction
tool, producing two types of outputs from marked point pattern data: (i) graphical visual-
izations of functional statistics and (ii) scalar-valued indices. Functional outputs capture
distance-based spatial dynamics (e.g., Ripley’s K-function, pair correlation function), while
scalar outputs summarize spatial autocorrelation and similarity across tissue slides (e.g.,
Moran’s I, Jaccard index, Summary statistics of persistence diagram).

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
├── src/                            # Core computational modules
│   ├── features/
│   │   ├── areal.R                # Areal Data Indices computation
│   │   ├── functional.R           # Spatial Summary Statistics
│   │   └── topological.py         # Topological Features
│   └── utils/
│       └── helpers.R              # Helper functions (normalization, regularization)
│
├── data/                          # Input datasets
│   ├── example_point_patternA.csv # Raw example data
│   ├── example_point_patternB.csv
│   ├── example_imageA.png
│   └── example_imageB.svg         # Visualization of example data
│    
├── examples/                      # Usage examples and workflows
│   ├── workflows/
│   │   ├── 01_areal_functional_workflow.R
│   │   └── 02_topological_workflow.ipynb
│   └── outputs/                   # Example outputs
|       ├── example_A_areal.csv
|       ├── example_A_K_function.csv
|       ├── example_B_topological.csv
|       └── example_B_topological_by_dimension.csv       
│
├── README.md
├── LICENSE
└── DESCRIPTION                  
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
