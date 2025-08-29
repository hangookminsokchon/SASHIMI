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


July 15, 2025:

- Areal features module completed (additional features may be added in future releases)
- Functional features module implemented
- File naming updates for improved consistency:

 *areal_data_features → areal_feature*
 
 *functional_data_features → functional_features*
 
 *compute_features → feature_computation* 

Aug 29, 2025:

Updated new feature: topolofical_features.py
Since new feature requires Python and the corresponding packages, package dependencies and versions are also updates on both README.md
and DESCRIPTION.txt. For computational efficiency, 'exampleB' image was used to compute the topological features, which is different from the exampleA image
that was used to compute areal, functional data. Current version only supports scalar topological features, however functional version will be added in the future.

Updated file directory:
Now the repo is consist of 3 folders
/src: source codes, example workflow codes
/data: example data, point pattern data images and raw .csv files
/example: output data from those corresponding inputs from /data

Reminder: For all 3 features, column name specification is fixed now
(ex: x, y, Z_cell for exampleA, x, y, class for exampleB), so are the cell type names.
Flexible implementation for col names and cell-type names, I'm planning to add more helper functions.


### Web Interface Specifications

#### Areal Data  
**Input**: `n × 3` CSV file with columns:
- `x`, `y`: coordinates  
- `z`: cell type  
(*Max file size: 4MB. Example available in `/example` folder.*)

**Output**: `1 × m` DataFrame of scalar summary values extracted from the spatial pattern.


#### Functional Data  
**Input**: `n × 3` CSV file with columns:
- `x`, `y`: coordinates  
- `z`: cell type  
(*Max file size: 4MB. Example available in `/example` folder.*)

**Output**: '500 x 3' DataFrame of functional data, which can be ploted using **plot()** function.


#### Topological Data  
**Input**: `n × 3` CSV file with columns:
- `x`, `y`: coordinates  
- `z`: cell type  
(*Max file size: 4MB. Example available in `/example` folder.*)

**Output**: `1 × m` DataFrame of scalar summary values extracted from the spatial pattern.

