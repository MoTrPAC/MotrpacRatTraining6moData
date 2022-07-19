## code to prepare `OUTLIERS` dataset goes here
library(data.table)
library(sinew)
library(MotrpacRatTraining6moData)

OUTLIERS = fread("data-raw/outliers.tsv", sep='\t', header=T)
OUTLIERS[,viallabel := as.character(viallabel)]

# add all 1w, 2w vena cava samples
venacv = data.table(PHENO[PHENO$specimen.processing.sampletypedescription == "Aorta",])
table(venacv[,specimen.processing.sampletypedescription])
venacv = venacv[sex == "female" & group %in% c("1w","2w")]
venacv[,.(viallabel, sex, group)]
out2 = data.table(viallabel = as.character(venacv[,viallabel]),
                  tissue = "VENACV",
                  pid = venacv[,pid],
                  group = sprintf("female_%s", venacv[,group]),
                  reason = "BAT contamination")
OUTLIERS = rbindlist(list(OUTLIERS, out2), fill=T)
OUTLIERS[OUTLIERS==""] = NA

# sort
OUTLIERS = OUTLIERS[order(tissue, assay)]

# data.frame
OUTLIERS = as.data.frame(OUTLIERS)

usethis::use_data(OUTLIERS, overwrite = TRUE)
sinew::makeOxygen(OUTLIERS)
