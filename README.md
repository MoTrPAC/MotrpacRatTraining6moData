# MotrpacRatTraining6moData

<!-- badges: start -->
![R package version](https://img.shields.io/github/r-package/v/MoTrPAC/MotrpacRatTraining6moData?label=R%20package)
[![R-CMD-check](https://github.com/MoTrPAC/MotrpacRatTraining6moData/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/MoTrPAC/MotrpacRatTraining6moData/actions/workflows/R-CMD-check.yaml)
![Last commit](https://img.shields.io/github/last-commit/MoTrPAC/MotrpacRatTraining6moData/main)
<!-- badges: end -->

**IMPORTANT: Watch this repository to get email notifications about new releases, which
can include changes to the data. In the top-right corner of [this page](https://github.com/MoTrPAC/MotrpacRatTraining6moData), 
choose `Watch` > `Custom` > `Releases`.**  

***

## Table of Contents
* [Overview](#overview)
  * [About this package](#about-this-package)
  * [About MoTrPAC](#about-motrpac)
* [Installation](#installation)
  * [Troubleshooting](#troubleshooting)
  * [Last resort](#last-resort)
* [Explore available data objects](#explore-available-data-objects)
* [Access epigenomics data through Google Cloud Storage](#access-epigenomics-data-through-google-cloud-storage)
* [Getting help](#getting-help)
* [Acknowledgements](#acknowledgements)
* [Data use agreement](#data-use-agreement)
* [Citing MoTrPAC data](#citing-motrpac-data)

## Overview

### About this package 
This package provides convenient access to the processed data and downstream
analysis results presented in the main paper for the first 
large-scale multi-omic multi-tissue endurance exercise training study conducted 
in young adult rats by the Molecular Transducers of Physical Activity Consortium 
(MoTrPAC). Find the [preprint on bioRxiv](https://www.biorxiv.org/content/10.1101/2022.09.21.508770v2).
**See the [vignette](https://motrpac.github.io/MotrpacRatTraining6moData/articles/MotrpacRatTraining6moData.html) for examples of how to use this package.** 

While the data in this package can be used by themselves, the 
[MotrpacRatTraining6mo](https://motrpac.github.io/MotrpacRatTraining6mo/)
R package relies heavily on this package and provides many functions to help
retrieve and explore the data. See examples in the 
[MotrpacRatTraining6mo vignette](https://motrpac.github.io/MotrpacRatTraining6mo/articles/MotrpacRatTraining6mo.html). 

### About MoTrPAC
MoTrPAC is a national research consortium designed to discover and perform 
preliminary characterization of the range of molecular transducers (the 
"molecular map") that underlie the effects of physical activity in humans. 
The program's goal is to study the molecular changes that occur during and after 
exercise and ultimately to advance the understanding of how physical activity 
improves and preserves health. The six-year program is the largest targeted NIH 
investment of funds into the mechanisms of how physical activity improves health 
and prevents disease. See [motrpac.org](https://www.motrpac.org/) and 
[motrpac-data.org](https://motrpac-data.org/) for more details. 

## Installation
Install this package with `devtools`:
```r
if (!require("devtools", quietly = TRUE)){
  install.packages("devtools")
}
devtools::install_github("MoTrPAC/MotrpacRatTraining6moData")
```

The output for a successful installation looks something like this. 
Note that the `*** moving datasets to lazyload DB` step takes the longest (~5 minutes):
```
Downloading GitHub repo MoTrPAC/MotrpacRatTraining6moData@HEAD
✓  checking for file ‘.../MoTrPAC-MotrpacRatTraining6moData-1c6478a/DESCRIPTION’ ...
─  preparing ‘MotrpacRatTraining6moData’:
✓  checking DESCRIPTION meta-information
─  checking for LF line-endings in source and make files and shell scripts
─  checking for empty or unneeded directories
─  building ‘MotrpacRatTraining6moData_1.3.2.tar.gz’ (1.3s)
   
* installing *source* package ‘MotrpacRatTraining6moData’ ...
** using staged installation
** R
** data
*** moving datasets to lazyload DB
** inst
** byte-compile and prepare package for lazy loading
** help
*** installing help indices
** building package indices
** testing if installed package can be loaded from temporary location
** testing if installed package can be loaded from final location
** testing if installed package keeps a record of temporary installation path
* DONE (MotrpacRatTraining6moData)
```

### Troubleshooting
If you get this error:  
```
Downloading GitHub repo MoTrPAC/MotrpacRatTraining6moData@HEAD
Error in utils::download.file(url, path, method = method, quiet = quiet,  : 
  download from 'https://api.github.com/repos/MoTrPAC/MotrpacRatTraining6moData/tarball/HEAD' failed
```
Try extending the timeout:  
```r
options(timeout=1e5)
devtools::install_github("MoTrPAC/MotrpacRatTraining6moData")
```

If you get this error after extending the `timeout`:  
```
Downloading GitHub repo MoTrPAC/MotrpacRatTraining6moData@HEAD
Error in utils::download.file(url, path, method = method, quiet = quiet,  : 
  download from 'https://api.github.com/repos/MoTrPAC/MotrpacRatTraining6moData/tarball/HEAD' failed
Error in `action()`:
! `class` is absent but must be supplied.
Run `rlang::last_error()` to see where the error occurred.
```
...this seems to be an intermittent issue seen only on Mac, not Linux or Windows. 
This was resolved by [installing the newest version of R](https://cran.rstudio.com/). 

### Last resort
If you can't get `devtools::install_github("MoTrPAC/MotrpacRatTraining6moData")` to work, try this:  

1. Go to <https://api.github.com/repos/MoTrPAC/MotrpacRatTraining6moData/tarball/HEAD>, which will automatically start downloading this repository in a tarball 
2. Install the package from source: 
   ```r
   install.packages("~/Downloads/MoTrPAC-MotrpacRatTraining6moData-0729e2e.tar.gz", 
     repos = NULL, 
     type = "source")
   library(MotrpacRatTraining6moData)
   ```

## Explore available data objects
Find the index of available objects [here](https://motrpac.github.io/MotrpacRatTraining6moData/reference/index.html).
Click on the name of an object in this index to see its documentation. 

Alternatively, explore data objects in R:
```r
# load the package
library(MotrpacRatTraining6moData)

# list available objects
data(package = "MotrpacRatTraining6moData")

# get the documentation for a given data object
?TRNSCRPT_LIVER_DA

# load a data object into your environment using `data()`
data(TRNSCRPT_LIVER_DA)
# or just use it directly, e.g., 
head(TRNSCRPT_LIVER_DA)
```

## Access epigenomics data through Google Cloud Storage

Due to file size, only normalized sample-level data and differential analysis results 
corresponding to **training-regulated features** (5% IHW FDR) are contained in this package
for chromatin accessibility (ATAC) and DNA methylation (METHYL). The full sets of epigenetic results
may be downloaded either with R functions in the [MotrpacRatTraining6mo package](https://motrpac.github.io/MotrpacRatTraining6mo/) 
or through the following public URLs. 

**Note that clicking on a [Link]() automatically starts a download.**
To instead copy the URL for an object, right-click the [Link]() and select `Copy Link Address`. 

To download and load epigenetic data within R, use one of the following functions in the
[MotrpacRatTraining6mo package](https://motrpac.github.io/MotrpacRatTraining6mo/): 

* [load_sample_data()](https://motrpac.github.io/MotrpacRatTraining6mo/reference/load_sample_data.html): For sample-level data from a single tissue and ome.  
* [combine_normalized_data()](https://motrpac.github.io/MotrpacRatTraining6mo/reference/combine_normalized_data.html): For sample-level data from multiple tissues or omes. Use `include_epigen = TRUE`.    
* [combine_da_results()](https://motrpac.github.io/MotrpacRatTraining6mo/reference/combine_da_results.html): For differential analysis results from multiple tissues or omes. Use `include_epigen = TRUE`.  
* Several other functions specifically for loading epigenetic data are documented [here](https://motrpac.github.io/MotrpacRatTraining6mo/reference/index.html#load-epigenetic-data).  

Note that the size in this table is the compressed size. Each object occupies several times more
memory when loaded into R. The total compressed size for all of these objects is 8.68 GiB (~9.32 GB).  

Type|Assay|Tissue|Object|Size (MiB)|Click to download|
---|---|---|---|---|---|
Differential analysis results|ATAC|BAT|ATAC_BAT_DA|329.71|[Link](https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/ATAC_BAT_DA.rda)
Differential analysis results|ATAC|HEART|ATAC_HEART_DA|238.99|[Link](https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/ATAC_HEART_DA.rda)
Differential analysis results|ATAC|HIPPOC|ATAC_HIPPOC_DA|223.06|[Link](https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/ATAC_HIPPOC_DA.rda)
Differential analysis results|ATAC|KIDNEY|ATAC_KIDNEY_DA|272.7|[Link](https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/ATAC_KIDNEY_DA.rda)
Differential analysis results|ATAC|LIVER|ATAC_LIVER_DA|202.48|[Link](https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/ATAC_LIVER_DA.rda)
Differential analysis results|ATAC|LUNG|ATAC_LUNG_DA|259.19|[Link](https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/ATAC_LUNG_DA.rda)
Differential analysis results|ATAC|SKM-GN|ATAC_SKMGN_DA|320.42|[Link](https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/ATAC_SKMGN_DA.rda)
Differential analysis results|ATAC|WAT-SC|ATAC_WATSC_DA|359.23|[Link](https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/ATAC_WATSC_DA.rda)
Normalized sample-level data|ATAC|BAT|ATAC_BAT_NORM_DATA|43.11|[Link](https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/ATAC_BAT_NORM_DATA.rda)
Normalized sample-level data|ATAC|HEART|ATAC_HEART_NORM_DATA|31.34|[Link](https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/ATAC_HEART_NORM_DATA.rda)
Normalized sample-level data|ATAC|HIPPOC|ATAC_HIPPOC_NORM_DATA|33.52|[Link](https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/ATAC_HIPPOC_NORM_DATA.rda)
Normalized sample-level data|ATAC|KIDNEY|ATAC_KIDNEY_NORM_DATA|40.03|[Link](https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/ATAC_KIDNEY_NORM_DATA.rda)
Normalized sample-level data|ATAC|LIVER|ATAC_LIVER_NORM_DATA|27.53|[Link](https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/ATAC_LIVER_NORM_DATA.rda)
Normalized sample-level data|ATAC|LUNG|ATAC_LUNG_NORM_DATA|34.54|[Link](https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/ATAC_LUNG_NORM_DATA.rda)
Normalized sample-level data|ATAC|SKM-GN|ATAC_SKMGN_NORM_DATA|43.31|[Link](https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/ATAC_SKMGN_NORM_DATA.rda)
Normalized sample-level data|ATAC|WAT-SC|ATAC_WATSC_NORM_DATA|45.98|[Link](https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/ATAC_WATSC_NORM_DATA.rda)
Raw sample-level counts|ATAC|BAT|ATAC_BAT_RAW_COUNTS|49.38|[Link](https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/ATAC_BAT_RAW_COUNTS.rda)
Raw sample-level counts|ATAC|HEART|ATAC_HEART_RAW_COUNTS|47.68|[Link](https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/ATAC_HEART_RAW_COUNTS.rda)
Raw sample-level counts|ATAC|HIPPOC|ATAC_HIPPOC_RAW_COUNTS|50.63|[Link](https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/ATAC_HIPPOC_RAW_COUNTS.rda)
Raw sample-level counts|ATAC|KIDNEY|ATAC_KIDNEY_RAW_COUNTS|52.6|[Link](https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/ATAC_KIDNEY_RAW_COUNTS.rda)
Raw sample-level counts|ATAC|LIVER|ATAC_LIVER_RAW_COUNTS|47.25|[Link](https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/ATAC_LIVER_RAW_COUNTS.rda)
Raw sample-level counts|ATAC|LUNG|ATAC_LUNG_RAW_COUNTS|48.92|[Link](https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/ATAC_LUNG_RAW_COUNTS.rda)
Raw sample-level counts|ATAC|SKM-GN|ATAC_SKMGN_RAW_COUNTS|50.96|[Link](https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/ATAC_SKMGN_RAW_COUNTS.rda)
Raw sample-level counts|ATAC|WAT-SC|ATAC_WATSC_RAW_COUNTS|48.91|[Link](https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/ATAC_WATSC_RAW_COUNTS.rda)
Feature annotation|ATAC|all|ATAC_FEATURE_ANNOT|27.73|[Link](https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/ATAC_FEATURE_ANNOT.rda)
Differential analysis results|METHYL|BAT|METHYL_BAT_DA|281.5|[Link](https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/METHYL_BAT_DA.rda)
Differential analysis results|METHYL|HEART|METHYL_HEART_DA|244.44|[Link](https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/METHYL_HEART_DA.rda)
Differential analysis results|METHYL|HIPPOC|METHYL_HIPPOC_DA|194.39|[Link](https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/METHYL_HIPPOC_DA.rda)
Differential analysis results|METHYL|KIDNEY|METHYL_KIDNEY_DA|138.43|[Link](https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/METHYL_KIDNEY_DA.rda)
Differential analysis results|METHYL|LIVER|METHYL_LIVER_DA|218.21|[Link](https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/METHYL_LIVER_DA.rda)
Differential analysis results|METHYL|LUNG|METHYL_LUNG_DA|340.03|[Link](https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/METHYL_LUNG_DA.rda)
Differential analysis results|METHYL|SKM-GN|METHYL_SKMGN_DA|278.02|[Link](https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/METHYL_SKMGN_DA.rda)
Differential analysis results|METHYL|WAT-SC|METHYL_WATSC_DA|220.35|[Link](https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/METHYL_WATSC_DA.rda)
Normalized sample-level data|METHYL|BAT|METHYL_BAT_NORM_DATA|76.76|[Link](https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/METHYL_BAT_NORM_DATA.rda)
Normalized sample-level data|METHYL|HEART|METHYL_HEART_NORM_DATA|59.43|[Link](https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/METHYL_HEART_NORM_DATA.rda)
Normalized sample-level data|METHYL|HIPPOC|METHYL_HIPPOC_NORM_DATA|56.76|[Link](https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/METHYL_HIPPOC_NORM_DATA.rda)
Normalized sample-level data|METHYL|KIDNEY|METHYL_KIDNEY_NORM_DATA|38.66|[Link](https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/METHYL_KIDNEY_NORM_DATA.rda)
Normalized sample-level data|METHYL|LIVER|METHYL_LIVER_NORM_DATA|55.02|[Link](https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/METHYL_LIVER_NORM_DATA.rda)
Normalized sample-level data|METHYL|LUNG|METHYL_LUNG_NORM_DATA|83.75|[Link](https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/METHYL_LUNG_NORM_DATA.rda)
Normalized sample-level data|METHYL|SKM-GN|METHYL_SKMGN_NORM_DATA|71.99|[Link](https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/METHYL_SKMGN_NORM_DATA.rda)
Normalized sample-level data|METHYL|WAT-SC|METHYL_WATSC_NORM_DATA|57.15|[Link](https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/METHYL_WATSC_NORM_DATA.rda)
Raw data|METHYL|BAT|BAT_RAW_DATA|394.79|[Link](https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/METHYL_BAT_RAW_DATA.rda)
Raw data|METHYL|HEART|HEART_RAW_DATA|359.21|[Link](https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/METHYL_HEART_RAW_DATA.rda)
Raw data|METHYL|HIPPOC|HIPPOC_RAW_DATA|361.81|[Link](https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/METHYL_HIPPOC_RAW_DATA.rda)
Raw data|METHYL|KIDNEY|KIDNEY_RAW_DATA|363.32|[Link](https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/METHYL_KIDNEY_RAW_DATA.rda)
Raw data|METHYL|LIVER|LIVER_RAW_DATA|349.14|[Link](https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/METHYL_LIVER_RAW_DATA.rda)
Raw data|METHYL|LUNG|LUNG_RAW_DATA|349.14|[Link](https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/METHYL_LUNG_RAW_DATA.rda)
Raw data|METHYL|SKM-GN|SKMGN_RAW_DATA|378.84|[Link](https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/METHYL_SKMGN_RAW_DATA.rda)
Raw data|METHYL|WAT-SC|WATSC_RAW_DATA|366.74|[Link](https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/METHYL_WATSC_RAW_DATA.rda)
Raw sample-level counts|METHYL|BAT|METHYL_BAT_RAW_COUNTS|71.6|[Link](https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/METHYL_BAT_RAW_COUNTS.rda)
Raw sample-level counts|METHYL|HEART|METHYL_HEART_RAW_COUNTS|56.99|[Link](https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/METHYL_HEART_RAW_COUNTS.rda)
Raw sample-level counts|METHYL|HIPPOC|METHYL_HIPPOC_RAW_COUNTS|52.56|[Link](https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/METHYL_HIPPOC_RAW_COUNTS.rda)
Raw sample-level counts|METHYL|KIDNEY|METHYL_KIDNEY_RAW_COUNTS|36.7|[Link](https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/METHYL_KIDNEY_RAW_COUNTS.rda)
Raw sample-level counts|METHYL|LIVER|METHYL_LIVER_RAW_COUNTS|52.54|[Link](https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/METHYL_LIVER_RAW_COUNTS.rda)
Raw sample-level counts|METHYL|LUNG|METHYL_LUNG_RAW_COUNTS|77.17|[Link](https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/METHYL_LUNG_RAW_COUNTS.rda)
Raw sample-level counts|METHYL|SKM-GN|METHYL_SKMGN_RAW_COUNTS|67.99|[Link](https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/METHYL_SKMGN_RAW_COUNTS.rda)
Raw sample-level counts|METHYL|WAT-SC|METHYL_WATSC_RAW_COUNTS|53.91|[Link](https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/METHYL_WATSC_RAW_COUNTS.rda)
Feature annotation|METHYL|all|METHYL_FEATURE_ANNOT|152.77|[Link](https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/METHYL_FEATURE_ANNOT.rda)

## Getting help 
**See the [vignette](https://motrpac.github.io/MotrpacRatTraining6moData/articles/MotrpacRatTraining6moData.html) 
for examples of how to use this package.**
Still have questions? For questions, bug reporting, and data requests for this package, please 
[submit a new issue](https://github.com/MoTrPAC/MotrpacRatTraining6moData/issues) 
and include as many details as possible. 

If the concern is related to functions provided in the 
[MotrpacRatTraining6mo](https://github.com/MoTrPAC/MotrpacRatTraining6mo)
package, please submit an issue 
[here](https://github.com/MoTrPAC/MotrpacRatTraining6mo/issues) instead. 

## Acknowledgements 
MoTrPAC is supported by the National Institutes of Health (NIH) Common
Fund through cooperative agreements managed by the National Institute of Diabetes and
Digestive and Kidney Diseases (NIDDK), National Institute of Arthritis and Musculoskeletal
Diseases (NIAMS), and National Institute on Aging (NIA). 
Specifically, the MoTrPAC Study is supported by NIH grants U24OD026629 (Bioinformatics Center), 
U24DK112349, U24DK112342, U24DK112340, U24DK112341, U24DK112326, U24DK112331, U24DK112348 (Chemical Analysis Sites), 
U01AR071133, U01AR071130, U01AR071124, U01AR071128, U01AR071150, U01AR071160, U01AR071158 (Clinical Centers), 
U24AR071113 (Consortium Coordinating Center), U01AG055133, U01AG055137 and U01AG055135 (PASS/Animal Sites).

## Data use agreement 
Recipients and their Agents agree that in publications using **any** data from MoTrPAC public-use data sets 
they will acknowledge MoTrPAC as the source of data, including the version number of the data sets used, e.g.:

* Data used in the preparation of this article were obtained from the Molecular Transducers of Physical Activity 
Consortium (MoTrPAC) database, which is available for public access at motrpac-data.org. 
Specific datasets used are [version numbers].

* Data used in the preparation of this article were obtained from the Molecular Transducers of Physical Activity 
Consortium (MoTrPAC) MotrpacRatTraining6moData R package [version number]. 

## Citing MoTrPAC data 
MoTrPAC Study Group. 2022. Temporal dynamics of the multi-omic response to endurance exercise training across tissues. 
bioRxiv doi: 10.1101/2022.09.21.508770
