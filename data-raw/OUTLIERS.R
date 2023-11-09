## code to prepare `OUTLIERS` dataset goes here
library(data.table)
library(sinew)
library(MotrpacRatTraining6moData)

OUTLIERS = fread("data-raw/outliers.tsv", sep='\t', header=T)
OUTLIERS[,viallabel := as.character(viallabel)]

# fix tissue
OUTLIERS[,tissue_code := gsub(",.*","",tissue)]
# remove MSSM outliers
OUTLIERS = OUTLIERS[!grepl("mssm", tissue)]
# remove vena cava rows
OUTLIERS = OUTLIERS[!viallabel == "all"]
OUTLIERS[,tissue_code]
OUTLIERS[,tissue := TISSUE_CODE_TO_ABBREV[tissue_code]]

# fix assay
setnames(OUTLIERS, "platform", "assay_code")
OUTLIERS[,assay := ASSAY_CODE_TO_ABBREV[assay_code]]
OUTLIERS[is.na(assay)]
OUTLIERS[,platform := NA_character_]
OUTLIERS[assay_code == "immunoassay-rat-mag27plex", platform := "rat-mag27plex"]
OUTLIERS[assay_code == "immunoassay-rat-mag27plex", assay := "IMMUNO"]
OUTLIERS[assay_code == "immunoassay-rat-mag27plex", assay_code := "immunoassay"]

# fix group
OUTLIERS[,group := gsub(",", "_", group)]

# add all 1w, 2w vena cava samples
venacv = data.table(PHENO[PHENO$specimen.processing.sampletypedescription == "Vena Cava",])
table(venacv[,specimen.processing.sampletypedescription])
venacv = venacv[sex == "female" & group %in% c("1w","2w")]
venacv[,.(viallabel, sex, group)]
out2 = data.table(viallabel = as.character(venacv[,viallabel]),
                  tissue = "VENACV",
                  tissue_code = "t99-vena-cava",
                  pid = venacv[,pid],
                  group = sprintf("female_%s", venacv[,group]),
                  reason = "BAT contamination")
OUTLIERS = rbindlist(list(OUTLIERS, out2), fill=T)
OUTLIERS[OUTLIERS==""] = NA

# sort
OUTLIERS = OUTLIERS[order(tissue, assay)]

# reorder
OUTLIERS = OUTLIERS[,.(
  viallabel, 
  tissue,
  tissue_code,
  assay,
  assay_code,
  platform, 
  pid, 
  group, 
  reason
)]

# data.frame
OUTLIERS = as.data.frame(OUTLIERS)

usethis::use_data(OUTLIERS, overwrite = TRUE)
sinew::makeOxygen(OUTLIERS)
