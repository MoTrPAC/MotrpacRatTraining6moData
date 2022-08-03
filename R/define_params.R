# Internal functions used to define parameters and column descriptions.
# Please add functions in alphabetical order. 
# 
# To render:
#   tissue character, tissue abbreviation, one of [MotrpacRatTraining6moData::TISSUE_ABBREV]
#   assay character, assay abbreviation, one of [MotrpacRatTraining6moData::ASSAY_ABBREV]
# Type:
#   tissue @eval tissue()
#   assay @eval assay()

adj_p_value_da = function(){
  "double, adjusted p-value from 'p_value' column. P-values are BY-adjusted across all datasets within a given assay/ome."
}

assay = function(){
  "character, assay abbreviation, one of [MotrpacRatTraining6moData::ASSAY_ABBREV]"
}

assay_abbrevation = function(){
  assay()
}

assay_code = function(){
  "character, assay code used in data release. See [MotrpacBicQC::assay_codes]."
}

comparison_average_intensity = function(){
  "double, average intensity among the replicates in 'comparison_group'"
}

comparison_average_intensity_se = function(){
  "double, standard error of 'comparison_average_intensity'"
}

comparison_group = function(){
  "character, time point of trained animals compared to the sex-matched sedentary control animals, one of '1w', '2w', '4w', '8w'"
}

covariates = function(){
  "character, comma-separated list of adjustment variables or NA"
}

dataset = function(){
  paste("character, for immunoassay and metabolomics features,",
        "this variable specifies the immune panel (rat-myokine, rat-pituitary, rat-mag27plex,",
        "rat-metabolic, ADIPONECTIN, SERPIN-E) or metabolomics platform (metab-u-ionpneg,",
        "metab-u-lrpneg, metab-u-lrppos, metab-u-hilicpos, metab-u-rpneg, metab-u-rppos,",
        "metab-t-amines, metab-t-oxylipneg, metab-t-tca, metab-t-nuc, metab-t-acoa,",
        "metab-t-ka) the feature was measured in. 'meta-reg' specifies results from",
        "the metabolomics meta-regression for repeated features"
        )
}

dataset_metab = function(){
  paste("character,",
        "metabolomics platform (metab-u-ionpneg,",
        "metab-u-lrpneg, metab-u-lrppos, metab-u-hilicpos, metab-u-rpneg, metab-u-rppos,",
        "metab-t-amines, metab-t-oxylipneg, metab-t-tca, metab-t-nuc, metab-t-acoa,",
        "metab-t-ka) the feature was measured in. 'meta-reg' specifies results from",
        "the metabolomics meta-regression for repeated features."
  )
}

feature = function(){
  paste("character, unique feature identifier in the format '[MotrpacRatTrainingData::ASSAY_ABBREV];[MotrpacRatTrainingData::TISSUE_ABBREV];feature_ID'",
        "only for training-regulated features at 5% IHW FDR.",
        "For redundant differential features, 'feature_ID' is prepended with the specific platform to make unique identifiers.",
        "See [MotrpacRatTrainingData::REPEATED_FEATURES] for details.")
}

feature_ID = function(){
  "character, MoTrPAC feature identifier"
}

logFC = function(){
  "double, log fold-change where the numerator is 'comparison_group', e.g., '1w', and the denominator is the group of sex-matched sedentary control animals"
}

logFC_se = function(){
  "double, standard error of the log fold-change"
}

meta_reg_het_p = function(){
  "double, for metabolites with multiple measurements, the meta-regression heterogeneity p-value, where a smaller p-value indicates more disagreement between platforms. One value per feature (not per training group)."
}

meta_reg_pvalue = function(){
  "double, for metabolites with multiple measurements, the meta-regression p-value. One value per feature (not per training group)."
}

numNAs = function(){
  "integer, number of missing values in the current 'comparison_group' and sex-matched sedentary control animals"
}

ome = function(){
  assay()
}

p_value_da = function(){
  "double, unadjusted p-value for the difference between 'comparison_group' and the group of sex-matched sedentary control animals"
}

reference_average_intensity = function(){
  "double, average intensity among the replicates in the group of sex-matched sedentary control animals"
}

reference_average_intensity_se = function(){
  "double, standard error of 'reference_average_intensity'"
}

removed_samples = function(){
  "character, comma-separated list of vial labels excluded from differential analysis or NA"
}

selection_fdr = function(){
  paste("double, adjusted training p-value used to select training-regulated analytes.", 
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

tscore = function(){
  "double, t statistic"
}

viallabel = function(){
  "character, sample identifier"
}

zscore = function(){
  "double, z statistic"
}
