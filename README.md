# SASHIMI
SASHIMI (Spatial Analysis for Segmented Histopathology Images using Machine Intelligence) is a tool for capturing spatial cell-cell interactions in histopathology image using spatial summary statistics.

## Description
**SASHIMI** is a computational framework designed to extract quantitative, statistically robust summaries of cell-cell interactions within the Tumor Microenvironment (TME). Leveraging spatial statistics and point pattern analysis, SASHIMI translates raw cell coordinate data from histopathological images into meaningful scalar descriptors that characterize tissue organization and cellular architecture.

## Dependencies
This program requires the use of following R packages. 

```{r}
library(spatstat)
library(dispRity)
library(spdep)
```

## Web page
### ** To be added ** 

## Comments for Developers

The SASHIMI framework consists of three feature modules:

- **Functional**
- **Areal** *(currently available)*
- **Topological**

As of **June 3**, only the **Areal** feature module is implemented. Functional and Topological features are in development and expected to be released by **mid-to-late July**.

### Web Interface Specifications

#### Areal Data  
**Input**: `n × 3` CSV file with columns:
- `x`, `y`: coordinates  
- `z`: cell type  
(*Max file size: 4MB. Example available in `/example` folder.*)

**Output**: `1 × m` DataFrame of scalar summary values extracted from the spatial pattern.

#### Functional Data  
**To be added.**

#### Topological Data  
**To be added.**

