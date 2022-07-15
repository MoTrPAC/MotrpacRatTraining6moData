## code to prepare `OUTLIERS` dataset goes here
library(data.table)
library(sinew)
OUTLIERS = as.data.frame(fread("data-raw/outliers.tsv", sep='\t', header=T))
usethis::use_data(OUTLIERS, overwrite = TRUE)
sinew::makeOxygen(OUTLIERS)
