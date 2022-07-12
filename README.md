# MotrpacRatTraining6moData
Data for analysis of the MoTrPAC endurance exercise training study in 6-month-old rats

## Installation
This takes a few minutes: 
```r
devtools::install_github("MoTrPAC/MotrpacRatTraining6moData")
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
ACETYL_HEART_DA                                    Differential analysis of proteomics datasets
ACETYL_LIVER_DA                                    Differential analysis of proteomics datasets
ACETYL_SAMPLE_DATA                                 Protein acetylation normalized sample-level data
ASSAY_ABBREV                                       Assay or "ome" abbreviations
ASSAY_ABBREV_TO_CODE                               Assay abbreviation-to-code mapping
ASSAY_CODE_TO_ABBREV                               Assay code-to-abbreviation mapping
ASSAY_COLORS                                       Assay colors
ASSAY_ORDER                                        Assay order
ATAC_META                                          ATAC-seq metadata and QC
ATAC_SAMPLE_DATA                                   ATAC-seq normalized sample-level data for training-regulated features
FEATURE_TO_GENE                                    Feature-to-gene map
GENE_UNIVERSES                                     Gene-centric universes
GRAPH_COMPONENTS                                   Graph components
GRAPH_PW_ENRICH                                    Graph pathway enrichment results
GRAPH_STATES                                       'repfdr' state assignments
GROUP_COLORS                                       Intervention group colors
IMMUNO_ADRNL_DA                                    Differential analysis of multiplexed immunoassays
IMMUNO_BAT_DA                                      Differential analysis of multiplexed immunoassays
IMMUNO_COLON_DA                                    Differential analysis of multiplexed immunoassays
IMMUNO_CORTEX_DA                                   Differential analysis of multiplexed immunoassays
IMMUNO_HEART_DA                                    Differential analysis of multiplexed immunoassays
IMMUNO_HIPPOC_DA                                   Differential analysis of multiplexed immunoassays
IMMUNO_KIDNEY_DA                                   Differential analysis of multiplexed immunoassays
IMMUNO_LIVER_DA                                    Differential analysis of multiplexed immunoassays
IMMUNO_LUNG_DA                                     Differential analysis of multiplexed immunoassays
IMMUNO_OVARY_DA                                    Differential analysis of multiplexed immunoassays
IMMUNO_PLASMA_DA                                   Differential analysis of multiplexed immunoassays
IMMUNO_SAMPLE_DATA                                 Immunoassay normalized sample-level data
IMMUNO_SKMGN_DA                                    Differential analysis of multiplexed immunoassays
IMMUNO_SKMVL_DA                                    Differential analysis of multiplexed immunoassays
IMMUNO_SMLINT_DA                                   Differential analysis of multiplexed immunoassays
IMMUNO_SPLEEN_DA                                   Differential analysis of multiplexed immunoassays
IMMUNO_TESTES_DA                                   Differential analysis of multiplexed immunoassays
IMMUNO_WATSC_DA                                    Differential analysis of multiplexed immunoassays
METAB_ADRNL_DA                                     Differential analysis of merged metabolomics datasets
METAB_BAT_DA                                       Differential analysis of merged metabolomics datasets
METAB_COLON_DA                                     Differential analysis of merged metabolomics datasets
METAB_CORTEX_DA                                    Differential analysis of merged metabolomics datasets
METAB_HEART_DA                                     Differential analysis of merged metabolomics datasets
METAB_HIPPOC_DA                                    Differential analysis of merged metabolomics datasets
METAB_HYPOTH_DA                                    Differential analysis of merged metabolomics datasets
METAB_KIDNEY_DA                                    Differential analysis of merged metabolomics datasets
METAB_LIVER_DA                                     Differential analysis of merged metabolomics datasets
METAB_LUNG_DA                                      Differential analysis of merged metabolomics datasets
METAB_OVARY_DA                                     Differential analysis of merged metabolomics datasets
METAB_PLASMA_DA                                    Differential analysis of merged metabolomics datasets
METAB_SAMPLE_DATA                                  Metabolomics normalized sample-level data
METAB_SKMGN_DA                                     Differential analysis of merged metabolomics datasets
METAB_SKMVL_DA                                     Differential analysis of merged metabolomics datasets
METAB_SMLINT_DA                                    Differential analysis of merged metabolomics datasets
METAB_SPLEEN_DA                                    Differential analysis of merged metabolomics datasets
METAB_TESTES_DA                                    Differential analysis of merged metabolomics datasets
METAB_VENACV_DA                                    Differential analysis of merged metabolomics datasets
METAB_WATSC_DA                                     Differential analysis of merged metabolomics datasets
METHYL_META                                        RRBS metadata and QC
METHYL_SAMPLE_DATA                                 RRBS normalized sample-level data for training-regulated features
PATHWAY_PARENTS                                    KEGG and Reactome parent pathways
PHENO                                              Phenotypic data
PHOSPHO_CORTEX_DA                                  Differential analysis of proteomics datasets
PHOSPHO_HEART_DA                                   Differential analysis of proteomics datasets
PHOSPHO_KIDNEY_DA                                  Differential analysis of proteomics datasets
PHOSPHO_LIVER_DA                                   Differential analysis of proteomics datasets
PHOSPHO_LUNG_DA                                    Differential analysis of proteomics datasets
PHOSPHO_SAMPLE_DATA                                Protein phosphorylation normalized sample-level data
PHOSPHO_SKMGN_DA                                   Differential analysis of proteomics datasets
PHOSPHO_WATSC_DA                                   Differential analysis of proteomics datasets
PROT_CORTEX_DA                                     Differential analysis of proteomics datasets
PROT_HEART_DA                                      Differential analysis of proteomics datasets
PROT_KIDNEY_DA                                     Differential analysis of proteomics datasets
PROT_LIVER_DA                                      Differential analysis of proteomics datasets
PROT_LUNG_DA                                       Differential analysis of proteomics datasets
PROT_SAMPLE_DATA                                   Global proteomics normalized sample-level data
PROT_SKMGN_DA                                      Differential analysis of proteomics datasets
PROT_WATSC_DA                                      Differential analysis of proteomics datasets
RAT_TO_HUMAN_GENE                                  Rat-to-human gene ortholog map
REPEATED_FEATURES                                  Repeated feature info
REPFDR_INPUTS                                      'repfdr' inputs
REPFDR_RES                                         'repfdr' results
SEX_COLORS                                         Sex colors
TISSUE_ABBREV                                      Tissue abbreviations
TISSUE_ABBREV_TO_CODE                              Tissue abbreviation-to-code mapping
TISSUE_CODE_TO_ABBREV                              Tissue code-to-abbreviation mapping
TISSUE_COLORS                                      Tissue colors
TISSUE_ORDER                                       Tissue order
TRAINING_REGULATED_FEATURES                        Training-regulated features
TRNSCRPT_ADRNL_DA                                  Differential analysis of RNA-seq datasets
TRNSCRPT_ADRNL_RAW_COUNTS                          RNA-seq raw counts
TRNSCRPT_BAT_DA                                    Differential analysis of RNA-seq datasets
TRNSCRPT_BAT_RAW_COUNTS                            RNA-seq raw counts
TRNSCRPT_BLOOD_DA                                  Differential analysis of RNA-seq datasets
TRNSCRPT_BLOOD_RAW_COUNTS                          RNA-seq raw counts
TRNSCRPT_COLON_DA                                  Differential analysis of RNA-seq datasets
TRNSCRPT_COLON_RAW_COUNTS                          RNA-seq raw counts
TRNSCRPT_CORTEX_DA                                 Differential analysis of RNA-seq datasets
TRNSCRPT_CORTEX_RAW_COUNTS                         RNA-seq raw counts
TRNSCRPT_HEART_DA                                  Differential analysis of RNA-seq datasets
TRNSCRPT_HEART_RAW_COUNTS                          RNA-seq raw counts
TRNSCRPT_HIPPOC_DA                                 Differential analysis of RNA-seq datasets
TRNSCRPT_HIPPOC_RAW_COUNTS                         RNA-seq raw counts
TRNSCRPT_HYPOTH_DA                                 Differential analysis of RNA-seq datasets
TRNSCRPT_HYPOTH_RAW_COUNTS                         RNA-seq raw counts
TRNSCRPT_KIDNEY_DA                                 Differential analysis of RNA-seq datasets
TRNSCRPT_KIDNEY_RAW_COUNTS                         RNA-seq raw counts
TRNSCRPT_LIVER_DA                                  Differential analysis of RNA-seq datasets
TRNSCRPT_LIVER_RAW_COUNTS                          RNA-seq raw counts
TRNSCRPT_LUNG_DA                                   Differential analysis of RNA-seq datasets
TRNSCRPT_LUNG_RAW_COUNTS                           RNA-seq raw counts
TRNSCRPT_META                                      RNA-seq metadata and QC
TRNSCRPT_OVARY_DA                                  Differential analysis of RNA-seq datasets
TRNSCRPT_OVARY_RAW_COUNTS                          RNA-seq raw counts
TRNSCRPT_SAMPLE_DATA                               RNA-seq normalized sample-level data
TRNSCRPT_SKMGN_DA                                  Differential analysis of RNA-seq datasets
TRNSCRPT_SKMGN_RAW_COUNTS                          RNA-seq raw counts
TRNSCRPT_SKMVL_DA                                  Differential analysis of RNA-seq datasets
TRNSCRPT_SKMVL_RAW_COUNTS                          RNA-seq raw counts
TRNSCRPT_SMLINT_DA                                 Differential analysis of RNA-seq datasets
TRNSCRPT_SMLINT_RAW_COUNTS                         RNA-seq raw counts
TRNSCRPT_SPLEEN_DA                                 Differential analysis of RNA-seq datasets
TRNSCRPT_SPLEEN_RAW_COUNTS                         RNA-seq raw counts
TRNSCRPT_TESTES_DA                                 Differential analysis of RNA-seq datasets
TRNSCRPT_TESTES_RAW_COUNTS                         RNA-seq raw counts
TRNSCRPT_VENACV_DA                                 Differential analysis of RNA-seq datasets
TRNSCRPT_VENACV_RAW_COUNTS                         RNA-seq raw counts
TRNSCRPT_WATSC_DA                                  Differential analysis of RNA-seq datasets
TRNSCRPT_WATSC_RAW_COUNTS                          RNA-seq raw counts
UBIQ_HEART_DA                                      Differential analysis of proteomics datasets
UBIQ_LIVER_DA                                      Differential analysis of proteomics datasets
UBIQ_SAMPLE_DATA                                   Protein ubiquitynation normalized sample-level data
```
