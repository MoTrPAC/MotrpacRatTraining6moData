#!/bin/R
library(MotrpacRatTraining6moData)
library(data.table)
library(sinew)

# for each ome, copy final training-dea tables to GCP bucket
public_bucket = "gs://motrpac-rat-training-6mo-extdata"

# get names of all files
all = system("gsutil ls gs://motrpac-data-hub/pass1b-06/analysis/**/dea/*training-dea-fdr*", intern=TRUE)
# remove protein-corrected ACETYL, PHOSPHO, non-protein-corrected UBIQ
all = all[!grepl("prot-ac-protein-corrected|prot-ph-protein-corrected|prot-ub_training-dea-fdr", all)]

# do it by ome
peek = c()
patterns = c(names(ASSAY_CODE_TO_ABBREV), "metab-meta-reg", "prot-ub-protein-corrected")
for(assay_code in patterns){
  curr = all[grepl(paste0(assay_code,"_"), all)]
  if(assay_code == "metab-meta-reg"){
    abbrev = "METAB"
  }else if(assay_code == "prot-ub-protein-corrected"){
    abbrev = "UBIQ"
  }else{
    abbrev = ASSAY_CODE_TO_ABBREV[[assay_code]]
  }
  message(assay_code)
  message(abbrev)
  writeLines(curr)
  # copy to bucket
  if(abbrev == "METAB"){
    if(assay_code == "metab"){
      command = sprintf("gsutil -m cp %s gs://motrpac-rat-training-6mo-extdata/training-da/METAB/redundant/", paste0(curr, collapse=" "))
    }else if(assay_code == "metab-meta-reg"){
      command = sprintf("gsutil -m cp %s gs://motrpac-rat-training-6mo-extdata/training-da/METAB/meta-regression/", paste0(curr, collapse=" "))
    }else{
      stop("metab assay code not recognized")
    }
  }else{
    command = sprintf("gsutil -m cp %s gs://motrpac-rat-training-6mo-extdata/training-da/%s/", paste0(curr, collapse=" "), abbrev)
  }
  system(command)
  message("---")
  peek = c(peek, curr[1])
}

# read in first few lines of each type of table to make docs 
peek = na.omit(peek)
dtlist = list()
for(f in peek){
  system(sprintf("gsutil cp %s /tmp", f))
  dt = fread(file=sprintf("/tmp/%s", basename(f)), nrows = 2, header=T)
  for(c in c("removed_samples","removed_samples_male","removed_samples_female")){
    if(c %in% colnames(dt)){
      dt[,(c) := as.character(get(c))]
    }
  }
  dtlist[[f]] = dt
}
data = rbindlist(dtlist, fill=TRUE)
makeOxygen(data)

# make sure we have all our files
written = system("gsutil ls gs://motrpac-rat-training-6mo-extdata/training-da/**", intern=TRUE)
stopifnot(length(written) == length(all))
stopifnot(all(basename(written) %in% basename(all)))
stopifnot(all(basename(all) %in% basename(written)))

# # this was run in bash to remove date suffixes.
# # otherwise it is difficult to reconstruct the URL. 
# for file in $(gsutil ls gs://motrpac-rat-training-6mo-extdata/training-da/**/*training*); do
#   new=$(echo $file | sed "s/-fdr_.*/-fdr\.txt/")
#   gsutil mv $file $new
#   echo $new
# done
