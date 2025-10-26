# SASHIMI
SASHIMI (Spatial Analysis for Segmented Histopathology Images using Machine Intelligence) is a tool for capturing spatial cell-cell interactions in histopathology image using spatial summary statistics.

## Description
**SASHIMI** is a computational framework designed to extract quantitative, statistically robust summaries of cell-cell interactions within the tumor tissue architecture. Leveraging spatial statistics and point pattern analysis, SASHIMI translates raw cell coordinate data from histopathological images into meaningful scalar and functional descriptors that characterize tissue organization and cellular architecture.

### Key Features

- **Spatial Summary Statistics**: Capture local and global spatial patterns of cellular distributions
- **Areal/Autocorrelation/Similarity Indices**: Quantify spatial dependencies and clustering behaviors
- **Topological Data Analysis**: Extract persistent homological features representing tissue architecture


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

#### Areal Data  
- **Format**: 1 × m DataFrame
- **Content**: Areal/Autocorrelation/Similarity indicies in scalar values.
- **Example metrics**: Moran's I, Geary's C, Cosine similarity
  
#### Functional Data  
- **Format**: 500 × 3 DataFrame
- **Content**: Functional curves representing spatial relationships
- **Example metrics**: functional data of K-function, Pair correlation function

#### Topological Data  
- **Format**: 1 × m DataFrame
- **Content**: Summary statistics of persistence diagram
- **Example metrics**: min/max/std of Betti 0, 1 numbers

## Repository Structure

```
SASHIMI/
├── src/                      # Source code and computation modules
│   ├── areal_features.R
│   ├── functional_features.R
│   ├── helperfunctions.R
│   ├── topological_features.py
│   ├── workflow_example.R
│   └── workflow_example_topological.ipynb
├── data/                     # Example input datasets
|   ├── example_imageA.png
|   ├── example_imageB.png
│   ├── example_point_patternA.csv         
│   └── example_point_patternB.csv         
├── example/                  # Example outputs
|   ├── output_example_areal.csv
|   ├── output_exampleA_functional_K
|   ├── output_exampleB_topological1.csv
│   └── output_exampleB_topological2.csv
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
