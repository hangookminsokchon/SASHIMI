# SASHIMI
SASHIMI (Spatial Analysis for Segmented Histopathology Images using Machine Intelligence) is a tool for capturing spatial cell-cell interactions in histopathology image using spatial summary statistics.

## Description
**SASHIMI** is a computational framework designed to extract quantitative, statistically robust summaries of cell-cell interactions within the Tumor Microenvironment (TME). Leveraging spatial statistics and point pattern analysis, SASHIMI translates raw cell coordinate data from histopathological images into meaningful scalar and functional descriptors that characterize tissue organization and cellular architecture.

## Dependencies
This framework requires the use of following software dependencies. 

R packages
```{r}
# R version 4.5.1
library(spatstat)  #(version >= 3.3)
library(dispRity)  #(version >= 1.9)
library(spdep)     #(version >= 1.3)
library(dplyr)     #(version >= 1.1.4)
```

Ptyhon packages
```{python}
# Python version 3.11
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

## Comments for Developers

The SASHIMI framework consists of three feature modules:

- **Functional** *(currently available)*
- **Areal** *(currently available)*
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


## Web Interface Specifications

#### Areal Data  
**Input**: `n × 3` CSV file with columns:
- `x`, `y`: coordinates  
- `type`: cell type  
(*Max file size: 4MB. Example available in `/data` folder.*)

**Output**: `1 × m` DataFrame of scalar summary values extracted from the spatial pattern.

#### Functional Data  
**Input**: `n × 3` CSV file with columns:
- `x`, `y`: coordinates  
- `type`: cell type  
(*Max file size: 4MB. Example available in `/data` folder.*)

**Output**: '500 x 3' DataFrame of functional data, which can be ploted using **plot()** function.

#### Topological Data  
**Input**: `n × 3` CSV file with columns:
- `x`, `y`: coordinates  
- `type`: cell type  
(*Max file size: 4MB. Example available in `/data` folder.*)

**Output**: `1 × m` DataFrame of scalar summary values extracted from the spatial pattern.

