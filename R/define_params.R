#' Internal functions used to define parameters and column descriptions.
#' Please add functions in alphabetical order. 
#' 
#' To render:
#'   tissue character, tissue abbreviation, one of [MotrpacRatTraining6moData::TISSUE_ABBREV]
#'   assay character, assay abbreviation, one of [MotrpacRatTraining6moData::ASSAY_ABBREV]
#' Type:
#'   tissue @eval tissue()
#'   assay @eval assay()

assay = function(){
  "character, assay abbreviation, one of [MotrpacRatTraining6moData::ASSAY_ABBREV]"
}

assay_code = function(){
  "character, assay code used in data release. See [MotrpacBicQC::assay_codes]."
}

comparison_group = function(){
  "time point of trained animals compared to the sex-matched sedentary control animals, one of '1w', '2w', '4w', '8w'"
}

feature = function(){
  paste("character, feature identifier in the format '[MotrpacRatTrainingData::ASSAY_ABBREV];[MotrpacRatTrainingData::TISSUE_ABBREV];feature_ID'.",
        "For redundant differential features, 'feature_ID' is prepended with the specific platform to make unique identifiers.",
        "See [MotrpacRatTrainingData::REPEATED_FEATURES] for details.")
}

feature_ID = function(){
  "character, MoTrPAC feature identifier"
}

logFC = function(){
  "numeric, log fold-change where the numerator is 'comparison_group', e.g., '1w', and the denominator is the group of sex-matched sedentary control animals"
}

logFC_se = function(){
  "numeric, standard error of the log fold-change"
}

ome = function(){
  assay()
}

selection_fdr = function(){
  paste("numeric, adjusted training p-value used to select training-regulated analytes.", 
        "P-values are IHW-adjusted across all datasets within a given assay with tissue as a covariate.")
}

sex = function(){
  "character, one of 'male' or 'female'"
}

tissue = function(){
  "character, tissue abbreviation, one of [MotrpacRatTraining6moData::TISSUE_ABBREV]"
}

tissue_abbreviation = function(){
  tissue()
}

tissue_code = function(){
  "character, tissue code used in data release. See [MotrpacBicQC::bic_animal_tissue_code]."
}

viallabel = function(){
  "sample identifier"
}
