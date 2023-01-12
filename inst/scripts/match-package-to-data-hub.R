#!/bin/R
# Nicole Gay
# 1/11/23

library(MotrpacRatTraining6mo)
library(MotrpacRatTraining6moData)
library(data.table)

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
                "category"="",
                "data_hub_path"="")

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

write.table(dt, "~/Desktop/data-r-portal.tsv", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
