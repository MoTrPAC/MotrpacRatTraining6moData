#!/bin/R
# Author: Nicole Gay
# 6 July 2022

# Set things up -----------------------------------------------------------------

library(MotrpacBicQC)
library(devtools)
library(usethis)
library(data.table)
library(sinew)

setwd("/oak/stanford/groups/smontgom/nicolerg/src/MOTRPAC/MotrpacRatTraining6moData")
data_dir = "/oak/stanford/groups/smontgom/shared/motrpac/internal_releases/motrpac-data-freeze-pass/v1.1/analysis"

# Define tissue codes 
tissue_codes = unique(bic_animal_tissue_code$tissue_name_release)
tissue_codes = tissue_codes[tissue_codes != ""]

# Helper functions --------------------------------------------------------------

#' Read a single file from Google Cloud into a data.table
#'
#' @param path GCP path, i.e. starts with "gs://"
#' @param sep column separator to use with "fread"
#' @param tmpdir scratch path to download files from GCP
#' @param GSUTIL_PATH path to "gsutil" on your computer
#' @param check_first check if file exists before downloading it. read in existing file if it exists. should be set to TRUE if you are running this function in parallel
#'
#' @return A data.table
dl_read_gcp = function(path,sep='\t',tmpdir='/oak/stanford/groups/smontgom/nicolerg/tmp',GSUTIL_PATH='~/google-cloud-sdk/bin/gsutil',check_first=F){
  system(sprintf('mkdir -p %s',tmpdir))
  # download
  new_path = sprintf('%s/%s',tmpdir,basename(path))
  # only download if it doesn't exist to avoid conflicts when running this script in parallel
  if(check_first){
    if(!file.exists(new_path)){
      cmd = sprintf('%s cp %s %s', GSUTIL_PATH, path, tmpdir)
      system(cmd,ignore.stdout = T,ignore.stderr = T)
    }
  }else{
    cmd = sprintf('%s cp %s %s', GSUTIL_PATH, path, tmpdir)
    system(cmd,ignore.stdout = T,ignore.stderr = T)
  }
  # read in the data as a data.table
  if(file.exists(new_path)){
    dt = fread(new_path,sep=sep,header=T)
    return(dt)
  }
  warning(sprintf("gsutil file %s does not exist.\n",path))
  return()
}

# Vectors and lists -------------------------------------------------------------

tissue_dt = data.table(MotrpacBicQC::bic_animal_tissue_code)

TISSUE_ORDER = MotrpacBicQC::tissue_order
TISSUE_ABBREV = TISSUE_ORDER[order(TISSUE_ORDER)]

# get tissue codes corresponding to just the tissue abbreviations 
tissue_codes = unique(tissue_dt[abbreviation %in% TISSUE_ABBREV, tissue_name_release])
tissue_codes = tissue_codes[tissue_codes!=""]
tissue_codes = tissue_codes[order(tissue_codes)]

# define tissue colors
TISSUE_COLORS = MotrpacBicQC::tissue_cols[c(TISSUE_ABBREV, tissue_codes)]

# define assay colors
library(RColorBrewer)
o = c("METAB","TRNSCRPT","PROT","ACETYL","PHOSPHO","UBIQ","ATAC","METHYL","IMMUNO")
assay_cols = brewer.pal(length(o), 'Set1')
names(assay_cols) = o
assay_cols[["UBIQ"]] = "#cfc100" # make UBIQ less bright yellow
ASSAY_COLORS = assay_cols
ASSAY_COLORS = ASSAY_COLORS[order(names(ASSAY_COLORS))]

# define assay abbreviations
ASSAY_ABBREV = o[order(o)]
ASSAY_ORDER = MotrpacBicQC::assay_order[MotrpacBicQC::assay_order %in% ASSAY_ABBREV]

# get assay codes corresponding to just the assay abbreviations 
assay_codes = c("epigen-atac-seq", "epigen-rrbs", "transcript-rna-seq", "metab", "prot-pr", "prot-ac", "prot-ph", "prot-ub", "immunoassay")
assay_abbr = c("ATAC", "METHYL", "TRNSCRPT", "METAB", "PROT", "ACETYL", "PHOSPHO", "UBIQ", "IMMUNO")

## define assay code <--> abbrev vectors

ASSAY_CODE_TO_ABBREV = assay_abbr
names(ASSAY_CODE_TO_ABBREV) = assay_codes
ASSAY_CODE_TO_ABBREV = ASSAY_CODE_TO_ABBREV[order(names(ASSAY_CODE_TO_ABBREV))]

ASSAY_ABBREV_TO_CODE = assay_codes
names(ASSAY_ABBREV_TO_CODE) = assay_abbr
ASSAY_ABBREV_TO_CODE = ASSAY_ABBREV_TO_CODE[order(names(ASSAY_ABBREV_TO_CODE))]

## define tissue code <--> abbrev vectors

tissue_map = tissue_dt[abbreviation %in% TISSUE_ABBREV & tissue_name_release != '']

TISSUE_CODE_TO_ABBREV = tissue_map[,abbreviation]
names(TISSUE_CODE_TO_ABBREV) = tissue_map[,tissue_name_release]
TISSUE_CODE_TO_ABBREV = TISSUE_CODE_TO_ABBREV[order(names(TISSUE_CODE_TO_ABBREV))]

TISSUE_ABBREV_TO_CODE = tissue_map[,tissue_name_release]
names(TISSUE_ABBREV_TO_CODE) = tissue_map[,abbreviation]
TISSUE_ABBREV_TO_CODE = TISSUE_ABBREV_TO_CODE[order(names(TISSUE_ABBREV_TO_CODE))]

GROUP_COLORS = MotrpacBicQC::group_cols[c("control","1w","2w","4w","8w")]
SEX_COLORS = MotrpacBicQC::sex_cols

# write each to a separate .rda
usethis::use_data(TISSUE_ABBREV,
                  TISSUE_COLORS,
                  ASSAY_COLORS,
                  ASSAY_ABBREV,
                  ASSAY_CODE_TO_ABBREV,
                  ASSAY_ABBREV_TO_CODE,
                  TISSUE_CODE_TO_ABBREV,
                  TISSUE_ABBREV_TO_CODE,
                  GROUP_COLORS,
                  SEX_COLORS,
                  ASSAY_ORDER,
                  TISSUE_ORDER, 
                  internal=FALSE, overwrite=TRUE)


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
sinew::makeOxygen(TRNSCRPT_BLOOD_DA)
sinew::makeOxygen(PROT_CORTEX_DA)
sinew::makeOxygen(PHOSPHO_CORTEX_DA)
sinew::makeOxygen(ACETYL_HEART_DA)
sinew::makeOxygen(UBIQ_HEART_DA)
sinew::makeOxygen(METAB_PLASMA_DA)
sinew::makeOxygen(IMMUNO_PLASMA_DA)

# Phenotypic data ---------------------------------------------------------------

dmaqc_metadata = 'gs://motrpac-data-freeze-pass/pass1b-06/v1.1/results/phenotype/pass1b_6m_viallabel_data.txt'
dmaqc_dict = 'gs://motrpac-data-freeze-pass/pass1b-06/v1.1/results/phenotype/merged_dictionary.txt'

# download and format phenotypic data 
dmaqc_metadata = dl_read_gcp(dmaqc_metadata)
dmaqc_dict = dl_read_gcp(dmaqc_dict)

# make a map from old column name to new column name
code_to_string = list()
for (c in colnames(dmaqc_metadata)){
  if(c %in% dmaqc_dict[,BICUniqueID]){
    # is there an exact match?
    corresponding_code = c
  }else{
    # is there a partial match?
    corresponding_code = gsub("_.*", "", c)
  }
  string = dmaqc_dict[BICUniqueID == corresponding_code, FullName]
  if(grepl("_", c)){
    string = paste0(string, "_", gsub(".*_","",c))
  }
  code_to_string[[c]] = string
}
# rename columns
code_to_string = unlist(code_to_string, use.names=TRUE)
all(colnames(dmaqc_metadata) %in% names(code_to_string))
colnames(dmaqc_metadata) = code_to_string[colnames(dmaqc_metadata)]

# make categorical variables human-readable
for (var in unique(dmaqc_dict[`Categorical.Definitions`!="",FullName])){
  if(var == "Key.d_sacrifice") next
  
  d = dmaqc_dict[FullName == var]
  keys = unname(unlist(strsplit(d[,Categorical.Values],'\\|')))
  values = tolower(unname(unlist(strsplit(d[,Categorical.Definitions],'\\|'))))
  names(values) = keys
  values = values[values!=""]
  
  # replace column with string
  if(!var %in% colnames(dmaqc_metadata)){
    # there are multiple columns that need to be converted
    cols_to_convert = colnames(dmaqc_metadata)[grepl(sprintf("^%s_",var), colnames(dmaqc_metadata))]
    print(cols_to_convert)
    for(column in cols_to_convert){
      dmaqc_metadata[,(column) := unname(values)[match(get(column), names(values))]]
    }
  }else{
    dmaqc_metadata[,(var) := unname(values)[match(get(var), names(values))]]
  }
}

# make column names lowercase for ease of use
colnames(dmaqc_metadata) = tolower(colnames(dmaqc_metadata))

# subset to PASS1B
pass1b = dmaqc_metadata[key.protocol=="phase 1b"]

## clean up some variables 
# clean up "sacrificetime"
pass1b[,sacrificetime := sapply(key.sacrificetime, function(x) gsub(' week.*','w',x))]
# clean up 'intervention'
pass1b[,intervention := key.intervention]
pass1b[grepl('training',intervention), intervention := 'training']
# make "group" - "1w", "2w", "4w", "8w", "control"
pass1b[,group := sacrificetime]
pass1b[intervention == 'control', group := 'control']
# make tech ID a string
pass1b[,specimen.collection.bloodtechid := paste0('tech',specimen.collection.bloodtechid)]
pass1b[,specimen.collection.uterustechid := paste0('tech',specimen.collection.uterustechid)]
pass1b[,specimen.processing.techid := paste0('tech',specimen.processing.techid)]
# make "sex" column
pass1b[,sex := registration.sex]
# make viallabel char
pass1b[,viallabel := as.character(viallabel)]

# add time_to_freeze
pass1b[,time_to_freeze := calculated.variables.frozetime_after_train - calculated.variables.deathtime_after_train]

# remove all-NA cols
cols_to_remove = c()
for(c in colnames(pass1b)){
  if(length(unique(pass1b[[c]]))==1){
    print(c)
    print(unique(pass1b[[c]]))
    if(is.na(unique(pass1b[[c]]))){
      cols_to_remove = c(cols_to_remove, c)
    }
  }
}
pass1b[,(cols_to_remove) := NULL]

# convert to data.frame
pass1b = as.data.frame(pass1b)
rownames(pass1b) = pass1b$viallabel
assign("PHENO", pass1b)
usethis::use_data(PHENO, overwrite = TRUE)
sinew::makeOxygen(PHENO)

 
# Feature-to-gene map -----------------------------------------------------------

system("gsutil cp gs://motrpac-data-freeze-pass/pass1b-06/v1.1/analysis/resources/master_feature_to_gene_20211116.RData /oak/stanford/groups/smontgom/nicolerg/tmp")
load("/oak/stanford/groups/smontgom/nicolerg/tmp/master_feature_to_gene_20211116.RData")

FEATURE_TO_GENE = as.data.frame(master_feature_to_gene)
REPEATED_FEATURES = as.data.frame(repeated_feature_map)

usethis::use_data(FEATURE_TO_GENE, overwrite = TRUE)
usethis::use_data(REPEATED_FEATURES, overwrite = TRUE)

sinew::makeOxygen(FEATURE_TO_GENE)
sinew::makeOxygen(REPEATED_FEATURES)


# Rat-to-human ortholog map -----------------------------------------------------

rat_to_human = dl_read_gcp("gs://mawg-data/external-datasets/rat-id-mapping/gencode.v39.RGD.20201001.human.rat.gene.ids.txt", sep="\t")
RAT_TO_HUMAN_GENE = as.data.frame(rat_to_human)
usethis::use_data(RAT_TO_HUMAN_GENE, overwrite = TRUE)
sinew::makeOxygen(RAT_TO_HUMAN_GENE)

# RNA-seq counts  ---------------------------------------------------------------

data_dir = "/oak/stanford/groups/smontgom/shared/motrpac/internal_releases/motrpac-data-freeze-pass/v1.1/results/transcriptomics"
for(tissue_code in names(TISSUE_CODE_TO_ABBREV)){
  file = sprintf("%s/%s/transcript-rna-seq/motrpac_pass1b-06_%s_transcript-rna-seq_rsem-genes-count.txt", data_dir, tissue_code, tissue_code)
  if(file.exists(file)){
    counts = fread(file, sep="\t", header=T)
    # set row names
    counts = as.data.frame(counts)
    rownames(counts) = counts$gene_id
    counts$gene_id = NULL
    # coerce values to int
    counts_round = as.data.frame(apply(counts, c(1,2), as.integer)) 
    # change name
    new_name = sprintf("TRNSCRPT_%s_RAW_COUNTS", gsub("-","",TISSUE_CODE_TO_ABBREV[[tissue_code]]))
    writeLines(new_name)
    assign(new_name, counts_round)
    # save .rda
    do.call("use_data", list(as.name(new_name), overwrite = TRUE))
  }else{
    print(sprintf("file %s DNE", file))
  }
}

# GET data type QC metrics ------------------------------------------------------

TRNSCRPT_META = as.data.frame(dl_read_gcp("gs://motrpac-data-freeze-pass/pass1b-06/v1.1/results/transcriptomics/qa-qc/motrpac_pass1b-06_transcript-rna-seq_qa-qc-metrics.csv", sep=","))
colnames(TRNSCRPT_META) = gsub(" .*","",colnames(TRNSCRPT_META))
ATAC_META = as.data.frame(dl_read_gcp("gs://motrpac-data-freeze-pass/pass1b-06/v1.1/results/epigenomics/qa-qc/motrpac_pass1b-06_epigen-atac-seq_qa-qc-metrics.csv", sep=","))
colnames(ATAC_META) = gsub(" .*","",colnames(ATAC_META))
METHYL_META = as.data.frame(dl_read_gcp("gs://motrpac-data-freeze-pass/pass1b-06/v1.1/results/epigenomics/qa-qc/motrpac_pass1b-06_epigen-rrbs_qa-qc-metrics.csv", sep=","))
colnames(METHYL_META) = gsub(" .*","",colnames(METHYL_META))

usethis::use_data(TRNSCRPT_META, ATAC_META, METHYL_META, internal = FALSE, overwrite = TRUE)
sinew::makeOxygen(TRNSCRPT_META)
sinew::makeOxygen(ATAC_META)
sinew::makeOxygen(METHYL_META)
