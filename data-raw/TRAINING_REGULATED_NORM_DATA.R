library(MotrpacRatTraining6mo)
library(devtools)
library(data.table)
load_all() # load dev version of this library 

# compile sample-level data for all differential features 
diff_features = data.table(TRAINING_REGULATED_FEATURES)
features = diff_features[,feature]
features = unique(features)

# handle metabolomics meta-regression feature_IDs
# meta-regression "feature" have multiple corresponding feature_IDs in the sample data
# rename when filtering by feature below 
metab_map = data.table(METAB_FEATURE_ID_MAP)

# annotate features 
feature_dt = data.table::as.data.table(check_cluster_res_format(data.frame(feature=features, cluster="dummy")))
datasets = unique(feature_dt[,.(ome, tissue)])
omes = unique(datasets[,ome])
training_regulated_only = TRUE

# save versions both with and without outliers 
res = list()
res[["TRUE"]] = list()
res[["FALSE"]] = list()
for(i in 1:nrow(datasets)){
  .ome = datasets[i, ome]
  .tissue = datasets[i, tissue]
  # get sample-level data 
  for(exclude_outliers in c(TRUE, FALSE)){
    data = NULL
    if (.ome %in% c("ATAC","METHYL")){
      data = load_sample_data(.tissue, 
                              .ome, 
                              normalized=TRUE, 
                              training_regulated_only=training_regulated_only, 
                              exclude_outliers=exclude_outliers, 
                              scratchdir=scratchdir,
                              warnings=TRUE)
    }else{
      data = load_sample_data(.tissue, 
                              .ome, 
                              normalized=TRUE, 
                              training_regulated_only=FALSE, 
                              exclude_outliers=exclude_outliers, 
                              scratchdir=scratchdir,
                              warnings=TRUE)
    }
    if(is.null(data)) next
    # convert colnames to PID
    viallabel_cols = colnames(data)[grepl("^9", colnames(data))]
    if(length(viallabel_cols)>0){
      pids = viallabel_to_pid(viallabel_cols)
      stopifnot(length(pids) == length(viallabel_cols))
      stopifnot(length(pids) == length(unique(pids)))
      # rename columns 
      new_colnames = as.character(unname(pids[viallabel_cols]))
      colnames(data)[grepl("^[0-9]", colnames(data))] = new_colnames
    }
    # add new feature column 
    data = data.table::as.data.table(data)
    data[,new_feature := sprintf("%s;%s;%s", assay, tissue, feature_ID)]
    # select features 
    curr_feat = feature_dt[ome==.ome & tissue==.tissue, feature]
    data = data[feature %in% curr_feat | new_feature %in% curr_feat]
    data[,feature := ifelse(feature %in% curr_feat, feature, new_feature)]
    data[,new_feature := NULL]
    if(nrow(data)>0){
      # add to result
      res[[as.character(exclude_outliers)]][[sprintf("%s_%s",.ome,.tissue)]] = data
    }else{
      warning(sprintf("No unfiltered features for %s %s.", .ome, .tissue))
    }
  }
}

merged = list()
for(version in names(res)){
  sample_level_data = data.table::rbindlist(res[[version]], fill=TRUE)
  
  # check if features are present 
  if(!all(features %in% sample_level_data[,feature])){
    # what's missing?
    missing = features[!features %in% sample_level_data[,feature]]
    warning(sprintf("%s out of %s features were not found in the normalized sample-level data:\n%s", 
                    length(missing), 
                    length(features), 
                    paste(missing, collapse="\n")))
  }
  merged[[version]] = sample_level_data
}

# save each of these
TRAINING_REGULATED_NORM_DATA = merged[["FALSE"]]
TRAINING_REGULATED_NORM_DATA_NO_OUTLIERS = merged[["TRUE"]]

TRAINING_REGULATED_NORM_DATA = as.data.frame(TRAINING_REGULATED_NORM_DATA)
use_data(TRAINING_REGULATED_NORM_DATA, overwrite = TRUE)
tools::resaveRdaFiles(paths = 'data/TRAINING_REGULATED_NORM_DATA.rda')

TRAINING_REGULATED_NORM_DATA_NO_OUTLIERS = as.data.frame(TRAINING_REGULATED_NORM_DATA_NO_OUTLIERS)
use_data(TRAINING_REGULATED_NORM_DATA_NO_OUTLIERS, overwrite = TRUE)
tools::resaveRdaFiles(paths = 'data/TRAINING_REGULATED_NORM_DATA_NO_OUTLIERS.rda')
