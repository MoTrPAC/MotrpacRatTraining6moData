# MotrpacRatTraining6moData 1.7.0 (2022-12-27)

* Update phenotypic data (`PHENO`) to v3, which includes changes to the VO2max data as announced by the BIC. 
Retain old column names for backwards compatibility.  
* Copy all `training-dea` tables to GCP and add documentation (`TRAINING_DA`).  
* Add details of metabolomics meta-regression to `METAB_DA_METAREG`.  

# MotrpacRatTraining6moData 1.6.1 (2022-12-13)

* Fix shuffled PIDs in `METAB_NORM_DATA_FLAT` and `IMMUNO_NORM_DATA_FLAT` introduced by 
not reordering the list returned by `MotrpacRatTraining6mo::viallabel_to_pid()`. 
`METAB_NORM_DATA_NESTED` and `IMMUNO_NORM_DATA_NESTED` are unchanged. 
* Regenerate `TRAINING_REGULATED_NORM_DATA` and `TRAINING_REGULATED_NORM_DATA_NO_OUTLIERS` to account for changes above.  

# MotrpacRatTraining6moData 1.6.0 (2022-11-08)

* Add `TRAINING_REGULATED_NORM_DATA` and `TRAINING_REGULATED_NORM_DATA_NO_OUTLIERS`  
* Fix missing `feature` labels in `METAB_NORM_DATA_FLAT` 

# MotrpacRatTraining6moData 1.5.0 (2022-10-14)

* Add `RAT_TO_HUMAN_PHOSPHO` 

# MotrpacRatTraining6moData 1.4.0 (2022-10-10)

* Add `FEATURE_TO_GENE_FILT` for faster indexing 

# MotrpacRatTraining6moData 1.3.3 (2022-09-20)

First version with vignette. Additional minor improvements to documentation. 
