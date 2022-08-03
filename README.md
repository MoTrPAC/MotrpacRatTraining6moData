# MotrpacRatTraining6moData
Data for analysis of the MoTrPAC endurance exercise training study in 6-month-old rats

## Installation
```r
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
─  building ‘MotrpacRatTraining6moData_1.0.0.tar.gz’ (1.3s)
   
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
   install.packages("~/Downloads/MoTrPAC-MotrpacRatTraining6moData-0729e2e.tar.gz", repos = NULL, type = "source")
   library(MotrpacRatTraining6moData)
   ```

## Explore available data objects
List available data objects after loading the package:
```r
library(MotrpacRatTraining6moData)
data(package = "MotrpacRatTraining6moData")
```
Get the documentation for a given data object:
```r
?TRNSCRPT_LIVER_DA
```
Load a data object into your environment using `data()`:
```r
data(TRNSCRPT_LIVER_DA)
```

List of current sets:
```
ACETYL_HEART_DA                                Differential analysis of proteomics datasets
ACETYL_HEART_NORM_DATA                         Normalized protein acetylation data
ACETYL_LIVER_DA                                Differential analysis of proteomics datasets
ACETYL_LIVER_NORM_DATA                         Normalized protein acetylation data
ASSAY_ABBREV                                   Assay or "ome" abbreviations
ASSAY_ABBREV_TO_CODE                           Assay abbreviation-to-code mapping
ASSAY_CODE_TO_ABBREV                           Assay code-to-abbreviation mapping
ASSAY_COLORS                                   Assay colors
ASSAY_ORDER                                    Assay order
ATAC_BAT_NORM_DATA_05FDR                       Normalized ATAC-seq data for training-regulated features
ATAC_HEART_NORM_DATA_05FDR                     Normalized ATAC-seq data for training-regulated features
ATAC_HIPPOC_NORM_DATA_05FDR                    Normalized ATAC-seq data for training-regulated features
ATAC_KIDNEY_NORM_DATA_05FDR                    Normalized ATAC-seq data for training-regulated features
ATAC_LIVER_NORM_DATA_05FDR                     Normalized ATAC-seq data for training-regulated features
ATAC_LUNG_NORM_DATA_05FDR                      Normalized ATAC-seq data for training-regulated features
ATAC_META                                      ATAC-seq metadata and QC
ATAC_SKMGN_NORM_DATA_05FDR                     Normalized ATAC-seq data for training-regulated features
ATAC_WATSC_NORM_DATA_05FDR                     Normalized ATAC-seq data for training-regulated features
FEATURE_TO_GENE                                Feature-to-gene map
GENE_UNIVERSES                                 Gene-centric universes
GRAPH_COMPONENTS                               Graph components
GRAPH_PW_ENRICH                                Graph pathway enrichment results
GRAPH_STATES                                   'repfdr' state assignments
GROUP_COLORS                                   Intervention group colors
IMMUNO_ADRNL_DA                                Differential analysis of multiplexed immunoassays
IMMUNO_BAT_DA                                  Differential analysis of multiplexed immunoassays
IMMUNO_COLON_DA                                Differential analysis of multiplexed immunoassays
IMMUNO_CORTEX_DA                               Differential analysis of multiplexed immunoassays
IMMUNO_HEART_DA                                Differential analysis of multiplexed immunoassays
IMMUNO_HIPPOC_DA                               Differential analysis of multiplexed immunoassays
IMMUNO_KIDNEY_DA                               Differential analysis of multiplexed immunoassays
IMMUNO_LIVER_DA                                Differential analysis of multiplexed immunoassays
IMMUNO_LUNG_DA                                 Differential analysis of multiplexed immunoassays
IMMUNO_META                                    Multiplexed immunoassay metadata and QC
IMMUNO_NORM_DATA_FLAT                          Combined immunoassay data used for visualization
IMMUNO_NORM_DATA_NESTED                        Processed immunoassay data used for differential analysis
IMMUNO_OVARY_DA                                Differential analysis of multiplexed immunoassays
IMMUNO_PLASMA_DA                               Differential analysis of multiplexed immunoassays
IMMUNO_SKMGN_DA                                Differential analysis of multiplexed immunoassays
IMMUNO_SKMVL_DA                                Differential analysis of multiplexed immunoassays
IMMUNO_SMLINT_DA                               Differential analysis of multiplexed immunoassays
IMMUNO_SPLEEN_DA                               Differential analysis of multiplexed immunoassays
IMMUNO_TESTES_DA                               Differential analysis of multiplexed immunoassays
IMMUNO_WATSC_DA                                Differential analysis of multiplexed immunoassays
METAB_ADRNL_DA                                 Differential analysis of merged metabolomics datasets
METAB_ADRNL_DA_METAREG                         Meta-regression of metabolomics differential analysis results
METAB_BAT_DA                                   Differential analysis of merged metabolomics datasets
METAB_BAT_DA_METAREG                           Meta-regression of metabolomics differential analysis results
METAB_COLON_DA                                 Differential analysis of merged metabolomics datasets
METAB_COLON_DA_METAREG                         Meta-regression of metabolomics differential analysis results
METAB_CORTEX_DA                                Differential analysis of merged metabolomics datasets
METAB_CORTEX_DA_METAREG                        Meta-regression of metabolomics differential analysis results
METAB_HEART_DA                                 Differential analysis of merged metabolomics datasets
METAB_HEART_DA_METAREG                         Meta-regression of metabolomics differential analysis results
METAB_HIPPOC_DA                                Differential analysis of merged metabolomics datasets
METAB_HIPPOC_DA_METAREG                        Meta-regression of metabolomics differential analysis results
METAB_HYPOTH_DA                                Differential analysis of merged metabolomics datasets
METAB_HYPOTH_DA_METAREG                        Meta-regression of metabolomics differential analysis results
METAB_KIDNEY_DA                                Differential analysis of merged metabolomics datasets
METAB_KIDNEY_DA_METAREG                        Meta-regression of metabolomics differential analysis results
METAB_LIVER_DA                                 Differential analysis of merged metabolomics datasets
METAB_LIVER_DA_METAREG                         Meta-regression of metabolomics differential analysis results
METAB_LUNG_DA                                  Differential analysis of merged metabolomics datasets
METAB_LUNG_DA_METAREG                          Meta-regression of metabolomics differential analysis results
METAB_OVARY_DA                                 Differential analysis of merged metabolomics datasets
METAB_OVARY_DA_METAREG                         Meta-regression of metabolomics differential analysis results
METAB_PLASMA_DA                                Differential analysis of merged metabolomics datasets
METAB_PLASMA_DA_METAREG                        Meta-regression of metabolomics differential analysis results
METAB_SAMPLE_DATA                              Nested metabolomics data used for differential analysis
METAB_SKMGN_DA                                 Differential analysis of merged metabolomics datasets
METAB_SKMGN_DA_METAREG                         Meta-regression of metabolomics differential analysis results
METAB_SKMVL_DA                                 Differential analysis of merged metabolomics datasets
METAB_SKMVL_DA_METAREG                         Meta-regression of metabolomics differential analysis results
METAB_SMLINT_DA                                Differential analysis of merged metabolomics datasets
METAB_SMLINT_DA_METAREG                        Meta-regression of metabolomics differential analysis results
METAB_SPLEEN_DA                                Differential analysis of merged metabolomics datasets
METAB_SPLEEN_DA_METAREG                        Meta-regression of metabolomics differential analysis results
METAB_TESTES_DA                                Differential analysis of merged metabolomics datasets
METAB_TESTES_DA_METAREG                        Meta-regression of metabolomics differential analysis results
METAB_VENACV_DA                                Differential analysis of merged metabolomics datasets
METAB_VENACV_DA_METAREG                        Meta-regression of metabolomics differential analysis results
METAB_VIZ_DATA                                 
METAB_WATSC_DA                                 Differential analysis of merged metabolomics datasets
METAB_WATSC_DA_METAREG                         Meta-regression of metabolomics differential analysis results
METHYL_BAT_NORM_DATA_05FDR                     Normalized DNA methylation data for training-regulated features
METHYL_HEART_NORM_DATA_05FDR                   Normalized DNA methylation data for training-regulated features
METHYL_HIPPOC_NORM_DATA_05FDR                  Normalized DNA methylation data for training-regulated features
METHYL_KIDNEY_NORM_DATA_05FDR                  Normalized DNA methylation data for training-regulated features
METHYL_LIVER_NORM_DATA_05FDR                   Normalized DNA methylation data for training-regulated features
METHYL_LUNG_NORM_DATA_05FDR                    Normalized DNA methylation data for training-regulated features
METHYL_META                                    RRBS metadata and QC
METHYL_SKMGN_NORM_DATA_05FDR                   Normalized DNA methylation data for training-regulated features
METHYL_WATSC_NORM_DATA_05FDR                   Normalized DNA methylation data for training-regulated features
OUTLIERS                                       Sample outliers
PATHWAY_PARENTS                                KEGG and Reactome parent pathways
PHENO                                          Phenotypic data
PHOSPHO_CORTEX_DA                              Differential analysis of proteomics datasets
PHOSPHO_CORTEX_NORM_DATA                       Normalized protein phosphorylation data
PHOSPHO_HEART_DA                               Differential analysis of proteomics datasets
PHOSPHO_HEART_NORM_DATA                        Normalized protein phosphorylation data
PHOSPHO_KIDNEY_DA                              Differential analysis of proteomics datasets
PHOSPHO_KIDNEY_NORM_DATA                       Normalized protein phosphorylation data
PHOSPHO_LIVER_DA                               Differential analysis of proteomics datasets
PHOSPHO_LIVER_NORM_DATA                        Normalized protein phosphorylation data
PHOSPHO_LUNG_DA                                Differential analysis of proteomics datasets
PHOSPHO_LUNG_NORM_DATA                         Normalized protein phosphorylation data
PHOSPHO_SKMGN_DA                               Differential analysis of proteomics datasets
PHOSPHO_SKMGN_NORM_DATA                        Normalized protein phosphorylation data
PHOSPHO_WATSC_DA                               Differential analysis of proteomics datasets
PHOSPHO_WATSC_NORM_DATA                        Normalized protein phosphorylation data
PROT_CORTEX_DA                                 Differential analysis of proteomics datasets
PROT_CORTEX_NORM_DATA                          Normalized protein expression data
PROT_HEART_DA                                  Differential analysis of proteomics datasets
PROT_HEART_NORM_DATA                           Normalized protein expression data
PROT_KIDNEY_DA                                 Differential analysis of proteomics datasets
PROT_KIDNEY_NORM_DATA                          Normalized protein expression data
PROT_LIVER_DA                                  Differential analysis of proteomics datasets
PROT_LIVER_NORM_DATA                           Normalized protein expression data
PROT_LUNG_DA                                   Differential analysis of proteomics datasets
PROT_LUNG_NORM_DATA                            Normalized protein expression data
PROT_SKMGN_DA                                  Differential analysis of proteomics datasets
PROT_SKMGN_NORM_DATA                           Normalized protein expression data
PROT_WATSC_DA                                  Differential analysis of proteomics datasets
PROT_WATSC_NORM_DATA                           Normalized protein expression data
RAT_TO_HUMAN_GENE                              Rat-to-human gene ortholog map
REPEATED_FEATURES                              Repeated feature info
REPFDR_INPUTS                                  'repfdr' inputs
REPFDR_RES                                     'repfdr' results
SEX_COLORS                                     Sex colors
TISSUE_ABBREV                                  Tissue abbreviations
TISSUE_ABBREV_TO_CODE                          Tissue abbreviation-to-code mapping
TISSUE_CODE_TO_ABBREV                          Tissue code-to-abbreviation mapping
TISSUE_COLORS                                  Tissue colors
TISSUE_ORDER                                   Tissue order
TRAINING_REGULATED_FEATURES                    Training-regulated features
TRNSCRPT_ADRNL_DA                              Differential analysis of RNA-seq datasets
TRNSCRPT_ADRNL_NORM_DATA                       Normalized RNA-seq data
TRNSCRPT_ADRNL_RAW_COUNTS                      RNA-seq raw counts
TRNSCRPT_BAT_DA                                Differential analysis of RNA-seq datasets
TRNSCRPT_BAT_NORM_DATA                         Normalized RNA-seq data
TRNSCRPT_BAT_RAW_COUNTS                        RNA-seq raw counts
TRNSCRPT_BLOOD_DA                              Differential analysis of RNA-seq datasets
TRNSCRPT_BLOOD_NORM_DATA                       Normalized RNA-seq data
TRNSCRPT_BLOOD_RAW_COUNTS                      RNA-seq raw counts
TRNSCRPT_COLON_DA                              Differential analysis of RNA-seq datasets
TRNSCRPT_COLON_NORM_DATA                       Normalized RNA-seq data
TRNSCRPT_COLON_RAW_COUNTS                      RNA-seq raw counts
TRNSCRPT_CORTEX_DA                             Differential analysis of RNA-seq datasets
TRNSCRPT_CORTEX_NORM_DATA                      Normalized RNA-seq data
TRNSCRPT_CORTEX_RAW_COUNTS                     RNA-seq raw counts
TRNSCRPT_HEART_DA                              Differential analysis of RNA-seq datasets
TRNSCRPT_HEART_NORM_DATA                       Normalized RNA-seq data
TRNSCRPT_HEART_RAW_COUNTS                      RNA-seq raw counts
TRNSCRPT_HIPPOC_DA                             Differential analysis of RNA-seq datasets
TRNSCRPT_HIPPOC_NORM_DATA                      Normalized RNA-seq data
TRNSCRPT_HIPPOC_RAW_COUNTS                     RNA-seq raw counts
TRNSCRPT_HYPOTH_DA                             Differential analysis of RNA-seq datasets
TRNSCRPT_HYPOTH_NORM_DATA                      Normalized RNA-seq data
TRNSCRPT_HYPOTH_RAW_COUNTS                     RNA-seq raw counts
TRNSCRPT_KIDNEY_DA                             Differential analysis of RNA-seq datasets
TRNSCRPT_KIDNEY_NORM_DATA                      Normalized RNA-seq data
TRNSCRPT_KIDNEY_RAW_COUNTS                     RNA-seq raw counts
TRNSCRPT_LIVER_DA                              Differential analysis of RNA-seq datasets
TRNSCRPT_LIVER_NORM_DATA                       Normalized RNA-seq data
TRNSCRPT_LIVER_RAW_COUNTS                      RNA-seq raw counts
TRNSCRPT_LUNG_DA                               Differential analysis of RNA-seq datasets
TRNSCRPT_LUNG_NORM_DATA                        Normalized RNA-seq data
TRNSCRPT_LUNG_RAW_COUNTS                       RNA-seq raw counts
TRNSCRPT_META                                  RNA-seq metadata and QC
TRNSCRPT_OVARY_DA                              Differential analysis of RNA-seq datasets
TRNSCRPT_OVARY_NORM_DATA                       Normalized RNA-seq data
TRNSCRPT_OVARY_RAW_COUNTS                      RNA-seq raw counts
TRNSCRPT_SKMGN_DA                              Differential analysis of RNA-seq datasets
TRNSCRPT_SKMGN_NORM_DATA                       Normalized RNA-seq data
TRNSCRPT_SKMGN_RAW_COUNTS                      RNA-seq raw counts
TRNSCRPT_SKMVL_DA                              Differential analysis of RNA-seq datasets
TRNSCRPT_SKMVL_NORM_DATA                       Normalized RNA-seq data
TRNSCRPT_SKMVL_RAW_COUNTS                      RNA-seq raw counts
TRNSCRPT_SMLINT_DA                             Differential analysis of RNA-seq datasets
TRNSCRPT_SMLINT_NORM_DATA                      Normalized RNA-seq data
TRNSCRPT_SMLINT_RAW_COUNTS                     RNA-seq raw counts
TRNSCRPT_SPLEEN_DA                             Differential analysis of RNA-seq datasets
TRNSCRPT_SPLEEN_NORM_DATA                      Normalized RNA-seq data
TRNSCRPT_SPLEEN_RAW_COUNTS                     RNA-seq raw counts
TRNSCRPT_TESTES_DA                             Differential analysis of RNA-seq datasets
TRNSCRPT_TESTES_NORM_DATA                      Normalized RNA-seq data
TRNSCRPT_TESTES_RAW_COUNTS                     RNA-seq raw counts
TRNSCRPT_VENACV_DA                             Differential analysis of RNA-seq datasets
TRNSCRPT_VENACV_NORM_DATA                      Normalized RNA-seq data
TRNSCRPT_VENACV_RAW_COUNTS                     RNA-seq raw counts
TRNSCRPT_WATSC_DA                              Differential analysis of RNA-seq datasets
TRNSCRPT_WATSC_NORM_DATA                       Normalized RNA-seq data
TRNSCRPT_WATSC_RAW_COUNTS                      RNA-seq raw counts
UBIQ_HEART_DA                                  Differential analysis of proteomics datasets
UBIQ_HEART_NORM_DATA                           Normalized protein ubiquitynation data
UBIQ_LIVER_DA                                  Differential analysis of proteomics datasets
UBIQ_LIVER_NORM_DATA                           Normalized protein ubiquitynation data
```
