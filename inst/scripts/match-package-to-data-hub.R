#!/bin/R
# Nicole Gay
# 1/11/23
# Updated 3/1/23

library(MotrpacRatTraining6mo) # also attaches MotrpacRatTraining6moData
library(data.table)
secret = "it's a secret"
secret2 = "it's also a secret"

#### Annotate objects in the data R package ####

objs = list_available_data(package="MotrpacRatTraining6moData")
types = c()

for(obj in objs){
  data = get(obj)
  if(is.data.frame(data)){
    type="data.frame"
  }else{
    type=typeof(data)
  }
  types = append(types, type)
}
table(types)

dt = data.table("r_object_name"=objs,
                "r_object_type"=types,
                "category"="")

# add some categories
dt[grepl("_DA$|TRAINING_REGULATED_FEATURES", r_object_name), category := "differential analysis"]
dt[grepl("_DA_METAREG", r_object_name), category := "differential analysis, meta-regression"]
dt[grepl("_NORM_DATA", r_object_name), category := "sample-level data, normalized"]
dt[grepl("_RAW_COUNTS", r_object_name), category := "sample-level data, raw"]
dt[grepl("_META$|OUTLIERS|PHENO", r_object_name), category := "sample-level metadata"]
dt[grepl("GRAPH|REPFDR", r_object_name), category := "graphical clustering"]
dt[grepl("ENRICH|PATHWAY|UNIVERSES", r_object_name), category := "pathway enrichment"]
dt[grepl("FEATURE_ANNOT|REPEATED|RAT_TO_|FEATURE_TO|FEATURE_ID", r_object_name), category := "feature mapping"]
dt[grepl("TISSUE_|ASSAY_|COLORS|ORDER", r_object_name), category := "abbreviations, colors, order"]

# add tissues
dt[,tissue := sapply(r_object_name, function(x){
  for(t in TISSUE_ABBREV){
    if(grepl(gsub("-","",t), x)){
      return(t)
    }
  }
  return("")
})]

# add omes
dt[,ome := sapply(r_object_name, function(x){
  for(t in ASSAY_ABBREV){
    if(grepl(t, x)){
      return(t)
    }
  }
  return("")
})]

dt = dt[order(category)]

# add column for data hub path
dt[,data_hub_path := ""]

#### Programmatically add some GCS URLs #### 

for(i in 1:nrow(dt)){
  
  # skip complete rows
  if(dt[i,data_hub_path] != ""){
    next
  }
  
  OBJ = dt[i, r_object_name]
  CAT = dt[i, category]
  TISSUE = dt[i, tissue] 
  OME = dt[i, ome]
  if(TISSUE=="" | OME==""){
    next
  }
  TISSUE_CODE = TISSUE_ABBREV_TO_CODE[[TISSUE]]
  ASSAY_CODE = ASSAY_ABBREV_TO_CODE[[OME]]
  
  if(CAT == "sample-level data, normalized"){
    if(OME %in% c("METAB","IMMUNO")){
      next
    }
    # normalized sample-level data
    if(OME == "ATAC"){
      command = sprintf("gsutil ls %s/analysis/**/normalized-data/**%s*%s*tsv",secret,ASSAY_CODE,TISSUE_CODE)
    }else{
      command = sprintf("gsutil ls %s/analysis/**/normalized-data/**%s*%s*txt",secret,TISSUE_CODE,ASSAY_CODE)
    }
    files = system(command, intern=TRUE)
    if(OME == "METHYL"){
      files = files[grepl("epigen-rrbs_normalized-log-M-window", files)] 
    }
  }else if(CAT == "sample-level data, raw"){
    # raw counts
    command = sprintf("gsutil ls %s/results/transcriptomics/**/motrpac_pass1b-06_%s_%s_rsem-genes-count.txt",secret,TISSUE_CODE,ASSAY_CODE)
    files = system(command, intern=TRUE)
  }else if(CAT == "differential analysis"){
    # standard differential analysis
    command = sprintf("gsutil ls %s/analysis/**/dea/*%s*%s*timewise-dea-fdr*",secret,TISSUE_CODE,ASSAY_CODE)
    files = system(command, intern=TRUE)
    if(OME=="METAB"){
      files = files[!grepl("meta-regression", files)]
    }
  }else if(CAT == "differential analysis, meta-regression"){
    # meta-regression
    command = sprintf("gsutil ls %s/analysis/metabolomics-named-merged/dea/meta-regression/*%s*timewise-dea-fdr*",secret,TISSUE_CODE)
    files = system(command, intern=TRUE)
  }
  
  if(length(files) > 1){
    if(OME %in% c("PHOSPHO","ACETYL")){
      files = files[!grepl("protein-corrected", files)]
      files = files[grepl("normalized-logratio", files)]
    }else if(OME == "UBIQ"){
      files = files[grepl("protein-corrected", files)]
      files = files[grepl("normalized-protein-corrected-logratio", files)]
    }else if(OME == "ATAC"){
      files = files[grepl("stanford", files)]
    }else if(OME == "PROT"){
      files = files[grepl("normalized-logratio", files)]
    }
  }
  
  if(length(files) == 0){
    warning(sprintf("No file: %s", OBJ))
  }else if(length(files) > 1){
    warning(sprintf("Multiple files: %s\n%s", OBJ, paste0(files, collapse="\n")))
  }else{
    message(sprintf("%s %s", OBJ, files))
    dt[i,data_hub_path := files]
  }
}

#### Write out table with annotated R package data objects #### 

write.table(dt, "~/Desktop/data-r-portal.tsv", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

#### Write a separate table of normalized data used for metabolomics

#' Read a single file from Google Cloud into a data.table
#'
#' @param path GCP path, i.e. starts with "gs://"
#' @param sep column separator to use with "fread"
#' @param tmpdir scratch path to download files from GCP
#' @param GSUTIL_PATH path to "gsutil" on your computer
#' @param check_first check if file exists before downloading it. read in existing file if it exists. should be set to TRUE if you are running this function in parallel
#'
#' @return A data.table
dl_read_gcp = function(path,sep='\t',tmpdir='/tmp',GSUTIL_PATH='gsutil',check_first=TRUE){
  system(sprintf('mkdir -p %s',tmpdir))
  # download
  new_path = sprintf('%s/%s',tmpdir,basename(path))
  # only download if it doesn't exist to avoid conflicts when running this script in parallel; clear scratch space when you're done
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

# metabolomics normalized data is a different beast
# partially copied from save_sample_level_data_metabolomics.Rmd
kw = dl_read_gcp(sprintf("%s/metabolomics/data/analysis_20210914/pass1b-06_sample-ctr-decision-table-kw-summary.txt", secret2))
tissue_codes = names(TISSUE_CODE_TO_ABBREV)
# list all normalized data
metab_norm = system(sprintf("gsutil ls %s/analysis/metabolomics-*/**/*metab*txt",secret), intern=TRUE)

get_file_for_dataset = function(d){
  possible_norm = metab_norm[grepl(d, metab_norm)]
  if(length(possible_norm) == 0){
    stop(sprintf("files not found for %s",d))
  }
  for(t in tissue_codes){
    if(grepl(t, d)){
      tissue = TISSUE_CODE_TO_ABBREV[[t]]
      tissue_code = t
      break
    }
  }
  platform = gsub("_","-",gsub(".*metab_","metab_",d))
  sample_ctr = kw[dataset == d, sample_ctr]
  data = METAB_NORM_DATA_NESTED[[platform]][[tissue]]
  if(is.null(data)){
    stop(sprintf("data not found for %s %s %s",d,platform,tissue))
  }
  if(grepl("metab-t", platform)){
    # if it's a targeted dataset
    # if >12 features, return table 3b imputed 3b.imputed_named-convert2na-log2-featurectr-featurefilt-knn
    # otherwise, return 2_named-convert2na-log2 
    if(nrow(data) > 12){
      file = possible_norm[grepl("imputed_named-convert2na-log2-featurectr-featurefilt-knn", possible_norm)]
    }else{
      file = possible_norm[grepl("2_named-convert2na-log2", possible_norm)]
    }
  }else{
    # if it's untargeted
    # if sample-centered, use 2b2 2b2_named-featurefilt-knn-samplefilt-log2-samplectr
    # otherwise use table 2 2_named-featurefilt-knn-samplefilt-log2
    if(sample_ctr==1){
      file = possible_norm[grepl("2b2_named-featurefilt-knn-samplefilt-log2-samplectr", possible_norm)]
    }else{
      file = possible_norm[grepl("2_named-featurefilt-knn-samplefilt-log2", possible_norm)]
    }
  }
  if(length(file) != 1){
    stop(sprintf("Not exactly one file: %s", paste(file, collapse=", ")))
  }
  return(list(tissue=tissue,
              platform=platform,
              dataset=d,
              file=file))
}

metab_files = list()
for(d in kw[,dataset]){
  metab_files[[d]] = get_file_for_dataset(d)
}
metab_files_dt = rbindlist(metab_files)
write.table(metab_files_dt, "~/Desktop/metab-norm-data-files.tsv", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

#### Write text files for objects that are missing from the data portal #### 

ATAC_FEATURE_ANNOT = load_atac_feature_annotation(scratchdir = "/tmp")
METHYL_FEATURE_ANNOT = load_methyl_feature_annotation(scratchdir = "/tmp")
TRNSCRPT_FEATURE_ANNOT = load_feature_annotation("TRNSCRPT", scratchdir = "/tmp")
tables = c(
  "TRAINING_REGULATED_FEATURES",
  "TRAINING_REGULATED_NORM_DATA",
  "TRAINING_REGULATED_NORM_DATA_NO_OUTLIERS",
  "GRAPH_STATES",
  "GRAPH_PW_ENRICH",
  "METAB_FEATURE_ID_MAP",
  "METAB_NORM_DATA_FLAT",
  "IMMUNO_NORM_DATA_FLAT",
  "RAT_TO_HUMAN_GENE",
  "OUTLIERS",
  "IMMUNO_META",
  "ATAC_FEATURE_ANNOT",
  "METHYL_FEATURE_ANNOT",
  "TRNSCRPT_FEATURE_ANNOT")

# start with tables because that's easy
for(table in tables){
  data = get(table)
  if(!is.data.frame(data)){
    warning(sprintf("%s is not a data frame", table))
    next
  }
  if(table == "GRAPH_PW_ENRICH"){
    data$parents = unlist(lapply(data$parents, function(x){
      paste(x, collapse=", ")
    }))
  }
  write.table(data, sprintf("~/Desktop/%s.txt", table), col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")
  message(table)
}

# handle colors, abbreviations, codes, order
library(MotrpacBicQC)

# "ASSAY_ABBREV"
# "ASSAY_ABBREV_TO_CODE"
# "ASSAY_CODE_TO_ABBREV"
# "ASSAY_COLORS"
# "ASSAY_ORDER"
assay_key = data.table(assay_codes)
# add a row for immuno 
assay_key = rbindlist(list(assay_key, 
                           data.table(omics_code="proteomics-targeted",
                                      submission_code="",
                                      assay_code="immunoassay",
                                      assay_name="Multiplexed immunoassays",
                                      cas_code="stanford")))
assay_key[,assay_abbreviation := ASSAY_CODE_TO_ABBREV[assay_code]]
assay_key[grepl("metab",omics_code), assay_abbreviation := "METAB"]
assay_key[,assay_hex_colour := ASSAY_COLORS[assay_abbreviation]]
assay_key = assay_key[!is.na(assay_abbreviation)]
# add assay order 
names(ASSAY_ORDER) = 1:length(ASSAY_ORDER)
assay_key[,assay_order := names(ASSAY_ORDER)[match(assay_abbreviation, ASSAY_ORDER)]]
write.table(assay_key, "~/Desktop/pass1b-06_assay_metadata.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")

# "TISSUE_ABBREV"
# "TISSUE_ABBREV_TO_CODE"
# "TISSUE_CODE_TO_ABBREV"
# "TISSUE_COLORS"
# "TISSUE_ORDER"
tissue_key = data.table(bic_animal_tissue_code)
# add tissue order
names(TISSUE_ORDER) = 1:length(TISSUE_ORDER)
tissue_key[,tissue_order := names(TISSUE_ORDER)[match(abbreviation, TISSUE_ORDER)]]
tissue_key = tissue_key[!is.na(tissue_order) & tissue_name_release != ""]
write.table(tissue_key, "~/Desktop/pass1b-06_tissue_metadata.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")

# handle colors separately/in addition
color_list = list()
for(obj in c("SEX_COLORS","GROUP_COLORS","ASSAY_COLORS","TISSUE_COLORS")){
  data = get(obj)
  dt = data.table(variable = tolower(gsub("_.*","",obj)),
                  value = names(data),
                  hex_colour = unname(data))
  color_list[[obj]] = dt
}
colors = rbindlist(color_list)
colors[hex_colour == "white", hex_colour := "#FFFFFF"]
write.table(colors, "~/Desktop/pass1b-06_color_codes.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")

#### Write text files for lists #### 

# "GENE_UNIVERSES"
names(GENE_UNIVERSES)
names(GENE_UNIVERSES$gene_symbol)
names(GENE_UNIVERSES$gene_symbol$IMMUNO)

# for each ID type, one column per tissue and ome 
for (id_type in names(GENE_UNIVERSES)){
  cols = list() # names: ome_tissue
  row_ome = c()
  row_tissue = c()
  longest = 0
  for(ome in names(GENE_UNIVERSES[[id_type]])){
    for(tissue in names(GENE_UNIVERSES[[id_type]][[ome]])){
      label = sprintf("%s_%s", ome, tissue)
      genes = GENE_UNIVERSES[[id_type]][[ome]][[tissue]]
      genes = genes[order(genes)]
      longest = max(longest, length(genes))
      row_ome = c(row_ome, ome)
      row_tissue = c(row_tissue, tissue)
      cols[[label]] = as.character(genes)
    }
  }
  # now extend each list to longest
  cols_filled = lapply(cols, function(x){
    c(x, rep(NA_character_, times=(longest-length(x))))
  })
  # now make it a data.table
  dt = data.table::copy(cols_filled)
  setDT(dt)   
  # add column headers
  header = data.table(V1 = row_ome, V2 = row_tissue)
  header = data.table(t(header))
  dt = rbindlist(list(header, dt), use.names=FALSE)
  
  write.table(dt, file=sprintf("~/Desktop/GENE_UNIVERSES_by_%s.txt",id_type), col.names=FALSE, row.names=FALSE, quote=FALSE, sep='\t')
}

# "PATHWAY_PARENTS"
head(names(PATHWAY_PARENTS))
dt = data.table(PATHWAY_ID = names(PATHWAY_PARENTS),
                PATHWAY_PARENTS = unname(unlist(PATHWAY_PARENTS)))
write.table(dt, file=sprintf("~/Desktop/PATHWAY_PARENTS.txt"), col.names=TRUE, row.names=FALSE, quote=FALSE, sep='\t')

# "REPFDR_INPUTS"
names(REPFDR_INPUTS)
for(f in names(REPFDR_INPUTS)){
  df = as.data.frame(REPFDR_INPUTS[[f]])
  if(!is.null(rownames(df)) & !'feature' %in% colnames(df)){
    df = cbind(feature=rownames(df), df)
  }
  write.table(df, file=sprintf("~/Desktop/REPFDR_INPUTS_%s.txt",f), col.names=TRUE, row.names=FALSE, quote=FALSE, sep='\t')
}

# "REPFDR_RES"
names(REPFDR_RES)
# repfdr_em_res
dt = data.table(cbind(data.table(feature=rownames(REPFDR_RES$repfdr_em_res$mat)),
                      REPFDR_RES$repfdr_em_res$mat))
write.table(dt, file=sprintf("~/Desktop/REPFDR_RES_repfdr_em_res_matrix.txt"), col.names=TRUE, row.names=FALSE, quote=FALSE, sep='\t')

dt = data.table(cbind(data.table(state=rownames(REPFDR_RES$repfdr_em_res$Pi)),
                      REPFDR_RES$repfdr_em_res$Pi))
write.table(dt, file=sprintf("~/Desktop/REPFDR_RES_repfdr_em_res_Pi.txt"), col.names=TRUE, row.names=FALSE, quote=FALSE, sep='\t')

# repfdr_clusters
dt = data.table(cbind(data.table(rowname=rownames(REPFDR_RES$repfdr_clusters)),
                      REPFDR_RES$repfdr_clusters))
write.table(dt, file=sprintf("~/Desktop/REPFDR_RES_repfdr_clusters.txt"), col.names=TRUE, row.names=FALSE, quote=FALSE, sep='\t')

# repfdr_clusters_str
# repfdr_clusters_pi
dt = data.table(cluster = names(REPFDR_RES$repfdr_clusters_pi),
                pi = REPFDR_RES$repfdr_clusters_pi)
write.table(dt, file=sprintf("~/Desktop/REPFDR_RES_repfdr_clusters_pi.txt"), col.names=TRUE, row.names=FALSE, quote=FALSE, sep='\t')
