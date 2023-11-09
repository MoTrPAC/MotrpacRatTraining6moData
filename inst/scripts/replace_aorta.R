library(data.table)
# library(devtools)
# load_all()

# copied from MotrpacRatTraining6mo
list_available_data = function(package=NULL){
  res = utils::data(package=package)
  obj = res$results[,3]
  # remove objects that can't be called directly
  obj = obj[!grepl("\\(", obj)]
  return(obj)
}

# list all available data objects
objs = list_available_data(package="MotrpacRatTraining6moData")

# identify data objects with "aorta" in them
has_aorta = c()
for (obj in objs){
  print(obj)
  # skip ATAC data, which doesn't have vena cava data
  if(grepl("ATAC", obj)){
    next()
  }
  # skip GENE_UNIVERSES, which already uses VENACV
  if(obj=="GENE_UNIVERSES"){
    next()
  }
  
  data = get(obj)
  if (is.data.frame(data)){
    if (any(data == "t65-aorta", na.rm=TRUE)){
      has_aorta = c(has_aorta, obj)
    }
    if (any(data == "aorta", na.rm=TRUE)){
      has_aorta = c(has_aorta, obj)
    }
    if (any(data == "Aorta", na.rm=TRUE)){
      has_aorta = c(has_aorta, obj)
    }
  }else if (is.vector(data)){
    if("t65-aorta" %in% data){
      has_aorta = c(has_aorta, obj)
    }
  }else{
    print(sprintf("Not sure what type %s is", obj))
  }
}

# manually edit data objects
unique(has_aorta)

data("METAB_VENACV_DA_METAREG")
METAB_VENACV_DA_METAREG$tissue_code[METAB_VENACV_DA_METAREG$tissue_code == "t65-aorta"] = "t99-vena-cava"
stopifnot(!any(METAB_VENACV_DA_METAREG == "t65-aorta", na.rm = TRUE))
usethis::use_data(METAB_VENACV_DA_METAREG, overwrite = T)

data("METAB_VENACV_DA")
METAB_VENACV_DA$tissue_code[METAB_VENACV_DA$tissue_code == "t65-aorta"] = "t99-vena-cava"
stopifnot(!any(METAB_VENACV_DA == "t65-aorta", na.rm = TRUE))
usethis::use_data(METAB_VENACV_DA, overwrite = T)

data("OUTLIERS")
OUTLIERS$tissue_code[OUTLIERS$tissue_code == "t65-aorta"] = "t99-vena-cava"
stopifnot(!any(OUTLIERS == "t65-aorta", na.rm = TRUE))
usethis::use_data(OUTLIERS, overwrite = T)

data("TISSUE_ABBREV_TO_CODE")
TISSUE_ABBREV_TO_CODE[TISSUE_ABBREV_TO_CODE == "t65-aorta"] = "t99-vena-cava"
stopifnot(!any(TISSUE_ABBREV_TO_CODE == "t65-aorta", na.rm = TRUE))
usethis::use_data(TISSUE_ABBREV_TO_CODE, overwrite = T)

data("TISSUE_CODE_TO_ABBREV")
names(TISSUE_CODE_TO_ABBREV)[names(TISSUE_CODE_TO_ABBREV) == "t65-aorta"] = "t99-vena-cava"
stopifnot(!any(TISSUE_CODE_TO_ABBREV == "t65-aorta", na.rm = TRUE))
stopifnot(!any(names(TISSUE_CODE_TO_ABBREV) == "t65-aorta", na.rm = TRUE))
usethis::use_data(TISSUE_CODE_TO_ABBREV, overwrite = T)

data("TRAINING_REGULATED_FEATURES")
TRAINING_REGULATED_FEATURES$tissue_code[TRAINING_REGULATED_FEATURES$tissue_code == "t65-aorta"] = "t99-vena-cava"
stopifnot(!any(TRAINING_REGULATED_FEATURES == "t65-aorta", na.rm = TRUE))
usethis::use_data(TRAINING_REGULATED_FEATURES, overwrite = T)

data("TRNSCRPT_VENACV_DA")
TRNSCRPT_VENACV_DA$tissue_code[TRNSCRPT_VENACV_DA$tissue_code == "t65-aorta"] = "t99-vena-cava"
stopifnot(!any(TRNSCRPT_VENACV_DA == "t65-aorta", na.rm = TRUE))
usethis::use_data(TRNSCRPT_VENACV_DA, overwrite = T)

data("PHENO")
PHENO[PHENO$tissue_code_no == "T65",]
PHENO$tissue_code_no[PHENO$tissue_code_no == "T65"] = "T99"
PHENO$tissue_description[PHENO$tissue_description == "Aorta"] = "Vena Cava"
PHENO$specimen.processing.aliquotdescription = gsub("Aorta", "Vena Cava", PHENO$specimen.processing.aliquotdescription)
PHENO$specimen.processing.sampletypedescription[PHENO$specimen.processing.sampletypedescription == "Aorta"] = "Vena Cava"
stopifnot(!any(PHENO == "t65-aorta", na.rm = TRUE))
usethis::use_data(PHENO, overwrite = T)
