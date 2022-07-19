# MotrpacRatTraining6moData
Data for analysis of the MoTrPAC endurance exercise training study in 6-month-old rats

## Installation
```r
devtools::install_github("MoTrPAC/MotrpacRatTraining6moData")
```
Note that the `*** moving datasets to lazyload DB` step takes the longest (~5 minutes). 

### Troubleshooting
If you get this error:
```
Downloading GitHub repo MoTrPAC/MotrpacRatTraining6moData@HEAD
Error in utils::download.file(url, path, method = method, quiet = quiet,  : 
  download from 'https://api.github.com/repos/MoTrPAC/MotrpacRatTraining6moData/tarball/HEAD' failed
```
Try extending the timeout: 
```r
devtools::install_github("MoTrPAC/MotrpacRatTraining6moData", timeout=1e5)
```

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
ATAC_META                                      ATAC-seq metadata and QC
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
IMMUNO_NORM_DATA                               Processed immunoassay data used for differential analysis
IMMUNO_OVARY_DA                                Differential analysis of multiplexed immunoassays
IMMUNO_PLASMA_DA                               Differential analysis of multiplexed immunoassays
IMMUNO_SKMGN_DA                                Differential analysis of multiplexed immunoassays
IMMUNO_SKMVL_DA                                Differential analysis of multiplexed immunoassays
IMMUNO_SMLINT_DA                               Differential analysis of multiplexed immunoassays
IMMUNO_SPLEEN_DA                               Differential analysis of multiplexed immunoassays
IMMUNO_TESTES_DA                               Differential analysis of multiplexed immunoassays
IMMUNO_VIZ_DATA                                Combined immunoassay data used for visualization
IMMUNO_WATSC_DA                                Differential analysis of multiplexed immunoassays
METAB_ADRNL_DA                                 Differential analysis of merged metabolomics datasets
METAB_BAT_DA                                   Differential analysis of merged metabolomics datasets
METAB_COLON_DA                                 Differential analysis of merged metabolomics datasets
METAB_CORTEX_DA                                Differential analysis of merged metabolomics datasets
METAB_HEART_DA                                 Differential analysis of merged metabolomics datasets
METAB_HIPPOC_DA                                Differential analysis of merged metabolomics datasets
METAB_HYPOTH_DA                                Differential analysis of merged metabolomics datasets
METAB_KIDNEY_DA                                Differential analysis of merged metabolomics datasets
METAB_LIVER_DA                                 Differential analysis of merged metabolomics datasets
METAB_LUNG_DA                                  Differential analysis of merged metabolomics datasets
METAB_OVARY_DA                                 Differential analysis of merged metabolomics datasets
METAB_PLASMA_DA                                Differential analysis of merged metabolomics datasets
METAB_SAMPLE_DATA                              Normalized metabolomics data
METAB_SKMGN_DA                                 Differential analysis of merged metabolomics datasets
METAB_SKMVL_DA                                 Differential analysis of merged metabolomics datasets
METAB_SMLINT_DA                                Differential analysis of merged metabolomics datasets
METAB_SPLEEN_DA                                Differential analysis of merged metabolomics datasets
METAB_TESTES_DA                                Differential analysis of merged metabolomics datasets
METAB_VENACV_DA                                Differential analysis of merged metabolomics datasets
METAB_VIZ_DATA                                 Combined metabolomics data
METAB_WATSC_DA                                 Differential analysis of merged metabolomics datasets
METHYL_META                                    RRBS metadata and QC
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
