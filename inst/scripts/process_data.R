#!/bin/R
# Author: Nicole Gay
# 6 July 2022

# Set things up -----------------------------------------------------------------

library(MotrpacBicQC)
library(devtools)
library(usethis)
library(data.table)
library(sinew)
load_all()

data_dir = "motrpac/internal_releases/motrpac-data-freeze-pass/v1.1"

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
dl_read_gcp = function(path,sep='\t',tmpdir='/tmp',GSUTIL_PATH='gsutil',check_first=F){
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
dea_directories = c('TRNSCRPT' = sprintf("%s/analysis/transcriptomics/transcript-rna-seq/dea", data_dir),
                    'PROT' = sprintf("%s/analysis/proteomics-untargeted/prot-pr/dea", data_dir),
                    'PHOSPHO' = sprintf("%s/analysis/proteomics-untargeted/prot-ph/dea", data_dir),
                    'ACETYL' = sprintf("%s/analysis/proteomics-untargeted/prot-ac/dea", data_dir),
                    'UBIQ' = sprintf("%s/analysis/proteomics-untargeted/prot-ub/dea", data_dir),
                    'IMMUNO' = sprintf("%s/analysis/proteomics-targeted/immunoassay-luminex/dea", data_dir),
                    'METAB' = sprintf("%s/analysis/metabolomics-named-merged/dea", data_dir),
                    'metab-metareg' = sprintf("%s/analysis/metabolomics-named-merged/dea/meta-regression", data_dir),
                    'METHYL' = sprintf("%s/analysis/epigenomics/epigen-rrbs/dea", data_dir),
                    'ATAC' = sprintf("%s/analysis/epigenomics/epigen-atac-seq/dea", data_dir))

dea_patterns = c('TRNSCRPT' = 'timewise-dea-fdr',
                 'PROT' = 'prot-pr_timewise-dea-fdr',
                 'PHOSPHO' = 'prot-ph_timewise-dea-fdr',
                 'ACETYL' = 'prot-ac_timewise-dea-fdr',
                 'UBIQ' = 'prot-ub-protein-corrected_timewise-dea-fdr',
                 'IMMUNO' = 'pass1b-06_immunoassay_timewise-dea-fdr_20211005.txt', # just one file
                 'METAB' = 'metab_timewise-dea-fdr',
                 'metab-metareg' = 'metab-meta-reg_timewise-dea-fdr',
                 'METHYL' = 'epigen-rrbs_timewise-dea-fdr',
                 'ATAC' = 'epigen-atac-seq_timewise-dea-fdr')

list_of_dfs = c()
#system("gsutil cp gs://motrpac-data-freeze-pass/pass1b-06/v1.1/analysis/resources/master_feature_to_gene_20211116.RData /tmp")
load("/tmp/master_feature_to_gene_20211116.RData")
REPEATED_FEATURES = as.data.table(repeated_feature_map)

load("data/TRAINING_REGULATED_FEATURES.rda") # from DEA tables, 5% FDR
differential_features = unique(TRAINING_REGULATED_FEATURES$feature)

rep = data.table(REPEATED_FEATURES)

for(ome in names(dea_directories)){
  if(ome == "IMMUNO") next
  if(ome == "metab-metareg"){
    assay_abbr = "METAB"
  }else{
    assay_abbr = ome
  }
  curr_ome_files = list.files(path = dea_directories[[ome]],
                              pattern = dea_patterns[[ome]],
                              full.names = T)
  for (tissue in tissue_codes){
    file = curr_ome_files[grepl(tissue, curr_ome_files)]
    if(length(file) == 0) next
    if(length(file) > 1) warning(sprintf("More than one file: %s", paste(file, collapse=", ")))
    dt = fread(file, sep="\t", header=T)

    # Get tissue abbreviation
    tissue_abbreviation = MotrpacBicQC::tissue_abbr[[tissue]]

    # add columns
    dt[,feature := sprintf("%s;%s;%s",
                           assay_abbr,
                           tissue_abbreviation,
                           feature_ID)]
    # rename some columns
    setnames(dt,
             old=c("tissue","assay","tissue_abbreviation"),
             new=c("tissue_code","assay_code","tissue"))
    dt[,assay := assay_abbr]

    if(assay_abbr=="METAB"){
      # fix repeated features
      features = dt[,feature]
      if(any(features %in% rep[,feature])){
        new_features = sprintf("METAB;%s;%s:%s", dt[,tissue], dt[,dataset], dt[,feature_ID])
        print(unique(new_features[new_features %in% rep[,new_feature]]))
        features[new_features %in% rep[,new_feature]] = new_features[new_features %in% rep[,new_feature]]
        dt[,feature := features]
      }
    }

    # reorder columns
    col_order = c("feature","assay","assay_code","tissue","tissue_code","feature_ID",
                  "dataset","site","is_targeted",
                  "sex","comparison_group",
                  "p_value","adj_p_value","logFC","logFC_se","shrunk_logFC","shrunk_logFC_se","tscore","zscore",
                  "covariates","removed_samples","numNAs",
                  "comparison_average_intensity","comparison_average_intensity_se",
                  "reference_average_intensity","reference_average_intensity_se",
                  "Chr", "Locus", "EntrezID", "Symbol", "fscore",
                  "metabolite_refmet","cv","metabolite","control_cv","mz",
                  "rt","neutral_mass","meta_reg_het_p","meta_reg_pvalue",
                  "selection_fdr")

    curr_cols = col_order[col_order %in% colnames(dt)]
    other_cols = colnames(dt)[!colnames(dt) %in% curr_cols]
    print(other_cols)
    dt = dt[,c(curr_cols, other_cols), with=F]

    # remove "feature" if it's not a differential feature
    dt[!feature %in% differential_features, feature := NA]

    # Convert to data.frame for ease of use
    df = as.data.frame(dt)

    # New name
    if(ome == "metab-metareg"){
      new_name = paste0(c(assay_abbr, gsub("-","",tissue_abbreviation), "DA_METAREG"), collapse="_")
    }else{
      new_name = paste0(c(assay_abbr, gsub("-","",tissue_abbreviation), "DA"), collapse="_")
    }

    writeLines(new_name)
    list_of_dfs = c(list_of_dfs, new_name)
    assign(new_name, df)

    if(assay_abbr %in% c("ATAC","METHYL")){
      # copy full set to GCS
      outfile = sprintf("%s/extracted_sample_level_data/%s/%s.rda",
                        data_dir,
                        assay_abbr,
                        new_name)
      do.call("save", list(as.name(new_name),
                           file = outfile,
                           compress = "bzip2",
                           compression_level = 9))
      print(outfile)
      print(head(get(new_name)))
    }else{
      print(new_name)
      print(head(get(new_name)))
      do.call("use_data", list(as.name(new_name), overwrite = TRUE))
    }
  }
}

## Fix missing refmet name for METAB
load("data/METAB_WATSC_DA_METAREG.rda")
METAB_WATSC_DA_METAREG$metabolite_refmet[METAB_WATSC_DA_METAREG$feature_ID == "TG(36:1)>TG(4:0_16:0_16:1)_and_TG(2:0_16:0_18:1)_feature4"] = "TG(36:1)_lp_d"
# save
use_data(METAB_WATSC_DA_METAREG, overwrite=T)

load("data/METAB_WATSC_DA.rda")
METAB_WATSC_DA$metabolite_refmet[METAB_WATSC_DA$feature_ID == "TG(36:1)>TG(4:0_16:0_16:1)_and_TG(2:0_16:0_18:1)_feature4"] = "TG(36:1)_lp_d"
# save
use_data(METAB_WATSC_DA, overwrite=T)

## Handle IMMUNO separately
ome = "IMMUNO"
assay_abbr = "IMMUNO"
# Read in single results file
dt = fread(sprintf("%s/%s", dea_directories[[ome]], dea_patterns[[ome]]), sep='\t', header=T)
for (curr_tissue in unique(dt[,tissue])){
  dt_sub = dt[tissue == curr_tissue]

  # Get tissue abbreviation
  tissue_abbreviation = MotrpacBicQC::tissue_abbr[[curr_tissue]]

  # add columns
  dt_sub[,feature := sprintf("IMMUNO;%s;%s",
                             tissue_abbreviation,
                             feature_ID)]
  # rename some columns
  setnames(dt_sub,
           old=c("tissue","assay","tissue_abbreviation","panel"),
           new=c("tissue_code","assay_code","tissue","dataset"))
  dt_sub[,assay := assay_abbr]

  # fix repeated features
  features = dt_sub[,feature]
  if(any(features %in% rep[,feature])){
    new_features = sprintf("IMMUNO;%s;%s:%s", dt_sub[,tissue], dt_sub[,dataset], dt_sub[,feature_ID])
    print(unique(new_features[new_features %in% rep[,new_feature]]))
    features[new_features %in% rep[,new_feature]] = new_features[new_features %in% rep[,new_feature]]
    dt_sub[,feature := features]
  }

  # reorder columns
  col_order = c("feature","assay","assay_code","tissue","tissue_code","feature_ID",
                "dataset","site","is_targeted",
                "sex","comparison_group",
                "p_value","adj_p_value","logFC","logFC_se","shrunk_logFC","shrunk_logFC_se","tscore","zscore",
                "covariates","removed_samples","numNAs",
                "comparison_average_intensity","comparison_average_intensity_se",
                "reference_average_intensity","reference_average_intensity_se",
                "Chr", "Locus", "EntrezID", "Symbol", "fscore",
                "metabolite_refmet","cv","metabolite","control_cv","mz",
                "rt","neutral_mass","meta_reg_het_p","meta_reg_pvalue",
                "selection_fdr")

  curr_cols = col_order[col_order %in% colnames(dt_sub)]
  other_cols = colnames(dt_sub)[!colnames(dt_sub) %in% curr_cols]
  print(other_cols)
  dt_sub = dt_sub[,c(curr_cols, other_cols), with=F]

  # remove non-DE feature
  dt_sub[!feature %in% differential_features, feature := NA]

  # Convert to data.frame for ease of use
  df = as.data.frame(dt_sub)

  new_name = paste0(c(assay_abbr, gsub("-","",tissue_abbreviation), "DA"), collapse="_")

  writeLines(new_name)
  list_of_dfs = c(list_of_dfs, new_name)
  assign(new_name, df)

  print(head(get(new_name)))
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
sinew::makeOxygen(METAB_PLASMA_DA_METAREG)

load(sprintf("%s/extracted_sample_level_data/ATAC/ATAC_BAT_DA.rda", data_dir))
sinew::makeOxygen(ATAC_BAT_DA)

load(sprintf("%s/extracted_sample_level_data/METHYL/METHYL_BAT_DA.rda", data_dir))
sinew::makeOxygen(METHYL_BAT_DA)

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

# add tissue column by parsing viallabel
pass1b[,tissue_code_no := sapply(viallabel, function(x){
  paste0(c("T", unname(unlist(strsplit(x, "")))[8:9]), collapse="")
})]
tissue_dt = data.table(MotrpacBicQC::bic_animal_tissue_code)
tissue_dt = tissue_dt[!is.na(abbreviation)]
pass1b[,tissue := tissue_dt[match(pass1b[,tissue_code_no], bic_tissue_code), abbreviation]]

# remove improperly calculated weight training var
pass1b[,calculated.variables.wgt_gain_after_train := NULL]

# convert to data.frame
pass1b = as.data.frame(pass1b)
rownames(pass1b) = pass1b$viallabel
assign("PHENO", pass1b)
usethis::use_data(PHENO, overwrite = TRUE)
sinew::makeOxygen(PHENO)


# Feature-to-gene map -----------------------------------------------------------

#system("gsutil cp gs://motrpac-data-freeze-pass/pass1b-06/v1.1/analysis/resources/master_feature_to_gene_20211116.RData /tmp")
load("/tmp/master_feature_to_gene_20211116.RData")

FEATURE_TO_GENE = as.data.frame(master_feature_to_gene)
REPEATED_FEATURES = repeated_feature_map

# rename some columns for consistency
head(REPEATED_FEATURES)
REPEATED_FEATURES[is.na(dataset) & !is.na(panel), dataset := panel]
REPEATED_FEATURES[,panel := NULL]
setnames(REPEATED_FEATURES, c("assay_abbr","tissue_abbreviation"), c("assay","tissue"))
REPEATED_FEATURES = as.data.frame(REPEATED_FEATURES)

usethis::use_data(FEATURE_TO_GENE, overwrite = TRUE)
usethis::use_data(REPEATED_FEATURES, overwrite = TRUE)

sinew::makeOxygen(FEATURE_TO_GENE)
sinew::makeOxygen(REPEATED_FEATURES)

# 10/10/22: Make filtered version of feature to gene map 
feat = data.table(FEATURE_TO_GENE)
data("TRAINING_REGULATED_FEATURES")
# remove non-diff epigen features
to_remove = feat[grepl("*chr.*[0-9]", feature_ID), feature_ID]
to_remove = to_remove[!to_remove %in% TRAINING_REGULATED_FEATURES$feature_ID]
feat = feat[!feature_ID %in% to_remove]
# save 
FEATURE_TO_GENE_FILT = as.data.frame(feat)
object.size(FEATURE_TO_GENE)
object.size(FEATURE_TO_GENE_FILT)
usethis::use_data(FEATURE_TO_GENE_FILT, overwrite = TRUE)
# re-compress 
tools::resaveRdaFiles(paths = 'data/FEATURE_TO_GENE_FILT.rda')


# Rat-to-human ortholog map -----------------------------------------------------

rat_to_human = dl_read_gcp("gs://mawg-data/external-datasets/rat-id-mapping/gencode.v39.RGD.20201001.human.rat.gene.ids.txt", sep="\t")
RAT_TO_HUMAN_GENE = as.data.frame(rat_to_human)
usethis::use_data(RAT_TO_HUMAN_GENE, overwrite = TRUE)
sinew::makeOxygen(RAT_TO_HUMAN_GENE)


# GET data type QC metrics ------------------------------------------------------

TRNSCRPT_META = as.data.frame(dl_read_gcp("gs://motrpac-data-freeze-pass/pass1b-06/v1.1/results/transcriptomics/qa-qc/motrpac_pass1b-06_transcript-rna-seq_qa-qc-metrics.csv", sep=","))
colnames(TRNSCRPT_META) = gsub(" .*","",colnames(TRNSCRPT_META))
ATAC_META = as.data.frame(dl_read_gcp("gs://motrpac-data-freeze-pass/pass1b-06/v1.1/results/epigenomics/qa-qc/motrpac_pass1b-06_epigen-atac-seq_qa-qc-metrics.csv", sep=","))
colnames(ATAC_META) = gsub(" .*","",colnames(ATAC_META))
METHYL_META = as.data.frame(dl_read_gcp("gs://motrpac-data-freeze-pass/pass1b-06/v1.1/results/epigenomics/qa-qc/motrpac_pass1b-06_epigen-rrbs_qa-qc-metrics.csv", sep=","))
colnames(METHYL_META) = gsub(" .*","",colnames(METHYL_META))

# refactor viallabel
TRNSCRPT_META = as.data.frame(cbind(viallabel = as.character(TRNSCRPT_META$vial_label),
                                    TRNSCRPT_META))
ATAC_META$viallabel = as.character(ATAC_META$viallabel)
METHYL_META = as.data.frame(cbind(viallabel = as.character(METHYL_META$vial_label),
                                  METHYL_META))

# remove ref stds
TRNSCRPT_META = TRNSCRPT_META[!grepl("^8",TRNSCRPT_META$viallabel),]
ATAC_META = ATAC_META[!grepl("^8",ATAC_META$viallabel),]
METHYL_META = METHYL_META[!grepl("^8",METHYL_META$viallabel),]

usethis::use_data(TRNSCRPT_META, ATAC_META, METHYL_META, internal = FALSE, overwrite = TRUE)
sinew::makeOxygen(TRNSCRPT_META)
sinew::makeOxygen(ATAC_META)
sinew::makeOxygen(METHYL_META)


# Proteomics sample-level data -------------------------------------------------

#gs://motrpac-data-freeze-pass/pass1b-06/v1.1/results/proteomics-untargeted/t58-heart/prot-pr/motrpac_pass1b-06_t58-heart_prot-pr_vial-metadata.txt
indir = sprintf("%s/results/extracted_proteomics", data_dir)
system(sprintf("gsutil -m cp gs://motrpac-data-freeze-pass/pass1b-06/v1.1/results/proteomics-untargeted/**/*vial-metadata.txt %s", indir))

for(assay_code in c('prot-pr','prot-ac','prot-ub','prot-ph')){
  dtlist = list()
  files = list.files(path=indir, pattern=assay_code, full.names=T)
  for(f in files){
    dt = fread(f, sep="\t", header=T)
    # add tissue label
    dt[,tissue := TISSUE_CODE_TO_ABBREV[[unname(unlist(strsplit(basename(f), '_')))[3]]]]
    dt[,viallabel := as.character(vial_label)]
    dt[,vial_label := NULL]
    dtlist[[f]] = dt
  }
  dt = rbindlist(dtlist)
  dt[,assay := ASSAY_CODE_TO_ABBREV[[assay_code]]]
  df = as.data.frame(dt)
  print(head(df))

  # save table
  new_name = sprintf("%s_META", ASSAY_CODE_TO_ABBREV[[assay_code]])

  assign(new_name, df)
  print(new_name)
  print(head(get(new_name)))
  do.call("use_data", list(as.name(new_name), overwrite = TRUE))
}


# RRBS feature annotation ------------------------------------------------------

indir = sprintf("%s/extracted_sample_level_data/METHYL", data_dir)
system(sprintf("gsutil -m cp gs://mawg-data/pass1b-06/epigen-rrbs-v2/data/*_epigen-rrbs_normalized-data-feature-annot.txt %s", indir))

dtlist = list()
files = list.files(path=indir, pattern="normalized-data-feature-annot", full.names=T)
for(f in files){
  dt = fread(f, sep="\t", header=T)
  # add tissue label
  dt[,tissue := TISSUE_CODE_TO_ABBREV[[unname(unlist(strsplit(basename(f), '_')))[3]]]]
  dt[,feature_ID := featureID]
  dt[,featureID := NULL]
  dtlist[[f]] = dt
  print(f)
}
dt = rbindlist(dtlist)

METHYL_FEATURE_ANNOT = as.data.frame(dt)
save(METHYL_FEATURE_ANNOT,
     file=sprintf("%s/extracted_sample_level_data/METHYL/METHYL_FEATURE_ANNOT.rda",data_dir),
     compress = "bzip2",
     compression_level = 9)

# ATAC feature annotation ------------------------------------------------------

atac_feature_annot = dl_read_gcp("gs://mawg-data/pass1b-06/epigen-atac-seq/mapping/pass1b-06_epigen-atac-seq_feature-mapping_20211110.txt")
setnames(atac_feature_annot, "assay", "assay_code")
atac_feature_annot = data.table(cbind(assay="ATAC", atac_feature_annot))
ATAC_FEATURE_ANNOT = as.data.frame(atac_feature_annot)
save(ATAC_FEATURE_ANNOT,
     file=sprintf("%s/extracted_sample_level_data/ATAC/ATAC_FEATURE_ANNOT.rda",data_dir),
     compress = "bzip2",
     compression_level = 9)
