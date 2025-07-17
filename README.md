# SASHIMI
SASHIMI (Spatial Analysis for Segmented Histopathology Images using Machine Intelligence) is a tool for capturing spatial cell-cell interactions in histopathology image using spatial summary statistics.

## Description
**SASHIMI** is a computational framework designed to extract quantitative, statistically robust summaries of cell-cell interactions within the Tumor Microenvironment (TME). Leveraging spatial statistics and point pattern analysis, SASHIMI translates raw cell coordinate data from histopathological images into meaningful scalar descriptors that characterize tissue organization and cellular architecture.

## Dependencies
This framework requires the use of following R packages. 

```{r}
# R version 4.5.1
library(spatstat)  #(version >= 3.3)
library(dispRity)  #(version >= 1.9)
library(spdep)     #(version >= 1.3)
library(dplyr)     #(version >= 1.1.4)
```

## Web page
### ** To be added ** 

## Comments for Developers

The SASHIMI framework consists of three feature modules:

- **Functional** *(currently available)*
- **Areal** *(currently available)*
- **Topological**

### Updates
June 3, 2024: Initial release with Areal feature module implementation. Functional and Topological feature modules scheduled for mid-to-late July release.

July 15, 2024:

- Areal features module completed (additional features may be added in future releases)
- Functional features module implemented
- File naming updates for improved consistency:

 *areal_data_features → areal_feature*
 
 *functional_data_features → functional_features*
 
 *compute_features → feature_computation* 



Additional functional features planned for future releases

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

**Output**: '500 x 2' DataFrame of functional data, which can be ploted using **plot()** function.

#### Topological Data  
**To be added.**

