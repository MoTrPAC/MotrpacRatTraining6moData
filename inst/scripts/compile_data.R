#!/bin/R
# Author: Nicole Gay
# 6 July 2022

# Set things up -----------------------------------------------------------------

library(MotrpacBicQC)
library(devtools)
library(usethis)
library(data.table)

setwd("/oak/stanford/groups/smontgom/nicolerg/src/MOTRPAC/MotrpacRatTraining6moData")
data_dir = "/oak/stanford/groups/smontgom/shared/motrpac/internal_releases/motrpac-data-freeze-pass/v1.1/analysis"

# Define tissue codes 
tissue_codes = unique(bic_animal_tissue_code$tissue_name_release)
tissue_codes = tissue_codes[tissue_codes != ""]

# Differential analysis results -------------------------------------------------

# Define directories
dea_directories = c('TRNSCRPT' = sprintf("%s/transcriptomics/transcript-rna-seq/dea", data_dir),
                    'PROT' = sprintf("%s/proteomics-untargeted/prot-pr/dea", data_dir),
                    'PHOSPHO' = sprintf("%s/proteomics-untargeted/prot-ph/dea", data_dir),
                    'ACETYL' = sprintf("%s/proteomics-untargeted/prot-ac/dea", data_dir),
                    'UBIQ' = sprintf("%s/proteomics-untargeted/prot-ub/dea", data_dir),
                    #'METHYL' = sprintf("%s/epigenomics/epigen-rrbs/dea", data_dir),
                    #'ATAC' = sprintf("%s/epigenomics/epigen-atac-seq/dea", data_dir),
                    'IMMUNO' = sprintf("%s/proteomics-targeted/immunoassay-luminex/dea", data_dir),
                    'METAB' = sprintf("%s/metabolomics-named-merged/dea", data_dir))

dea_patterns = c('TRNSCRPT' = 'timewise-dea-fdr',
                 'PROT' = 'prot-pr_timewise-dea-fdr',
                 'PHOSPHO' = 'prot-ph_timewise-dea-fdr',
                 'ACETYL' = 'prot-ac_timewise-dea-fdr',
                 'UBIQ' = 'prot-ub-protein-corrected_timewise-dea-fdr',
                 #'METHYL' = '',
                 #'ATAC' = '',
                 'IMMUNO' = 'pass1b-06_immunoassay_timewise-dea-fdr_20211005.txt', # just one file
                 'METAB' = 'metab_timewise-dea-fdr')                    

list_of_dfs = c()

for(ome in names(dea_directories)){
  if(ome == "IMMUNO") next
  curr_ome_files = list.files(path = dea_directories[[ome]], 
                              pattern = dea_patterns[[ome]],
                              full.names = T)
  for (tissue in tissue_codes){
    file = curr_ome_files[grepl(tissue, curr_ome_files)]
    if(length(file) == 0) next
    if(length(file) > 1) warning(sprintf("More than one file: %s", paste(file, collapse=", ")))
    dt = fread(file, sep="\t", header=T)
    # Convert to data.frame for ease of use 
    df = as.data.frame(dt) 
    # Get tissue abbreviation
    tissue_abbreviation = MotrpacBicQC::tissue_abbr[[tissue]]
    # New name
    new_name = paste0(c(ome, gsub("-","",tissue_abbreviation), "DA"), collapse="_")
    writeLines(new_name)
    list_of_dfs = c(list_of_dfs, new_name)
    assign(new_name, df)
    # Save .rda - this keeps giving me an error early in the loop
    do.call("use_data", list(as.name(new_name), overwrite = TRUE))
  }
}

## Handle IMMUNO separately
ome = "IMMUNO"
# Read in single results file 
dt = fread(sprintf("%s/%s", dea_directories[[ome]], dea_patterns[[ome]]), sep='\t', header=T)
for (curr_tissue in unique(dt[,tissue])){
  dtsub = dt[tissue == curr_tissue]
  # Convert to data.frame for ease of use 
  df = as.data.frame(dtsub) 
  # Get tissue abbreviation
  tissue_abbreviation = MotrpacBicQC::tissue_abbr[[curr_tissue]]
  # New name
  new_name = paste0(c(ome, gsub("-","",tissue_abbreviation), "DA"), collapse="_")
  writeLines(new_name)
  list_of_dfs = c(list_of_dfs, new_name)
  assign(new_name, df)
  # Save .rda - this keeps giving me an error early in the loop
  do.call("use_data", list(as.name(new_name), overwrite = TRUE))
}

# To help with documentation:
writeLines(list_of_dfs)
library(sinew)
sinew::makeOxygen(TRNSCRPT_BLOOD_DA)
sinew::makeOxygen(PROT_CORTEX_DA)
sinew::makeOxygen(PHOSPHO_CORTEX_DA)
sinew::makeOxygen(ACETYL_HEART_DA)
sinew::makeOxygen(UBIQ_HEART_DA)
sinew::makeOxygen(METAB_PLASMA_DA)
sinew::makeOxygen(IMMUNO_PLASMA_DA)
