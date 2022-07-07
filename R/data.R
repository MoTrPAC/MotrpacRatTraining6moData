#' @title Differential analysis of RNA-seq datasets 
#' @description Timewise summary statistics and training FDR from 
#'     differential analysis (DA) that tests the effect of training on the 
#'     expression of each gene within each sex. One data frame per tissue.   
#' @format A data frame with >60000 rows and 20 variables:
#' \describe{
#'   \item{\code{feature_ID}}{character, analyte identifier}
#'   \item{\code{tissue}}{character, MoTrPAC tissue release code}
#'   \item{\code{sex}}{character, one of "male" or "female"}
#'   \item{\code{comparison_group}}{character, time point compared to the sex-matched 
#'       sedentary control animals, one of "1w", "2w", "4w", "8w"}
#'   \item{\code{assay}}{character, MoTrPAC assay or "ome" code. "transcript-rna-seq" for RNA-seq datasets.}
#'   \item{\code{covariates}}{character, comma-separated list of covariates or NA}
#'   \item{\code{removed_samples}}{character, comma-separated list of vial labels excluded from differential analysis or NA}
#'   \item{\code{logFC}}{double, log fold-change where the numerator is \code{comparison_group} 
#'       and the denominator is the group of sex-matched sedentary control animals}
#'   \item{\code{logFC_se}}{double, standard error of the log fold-change}
#'   \item{\code{shrunk_logFC}}{double, log fold-change with shrinkage applied}
#'   \item{\code{shrunk_logFC_se}}{double, standard error of the shrunken log fold-change}
#'   \item{\code{zscore}}{double, z statistic}
#'   \item{\code{p_value}}{double, unadjusted p-value for the difference between 
#'       \code{comparison_group} and the group of sex-matched sedentary control animals}
#'   \item{\code{comparison_average_intensity}}{double, average intensity among the replicates in \code{comparison_group}}
#'   \item{\code{comparison_average_intensity_se}}{double, standard error of \code{comparison_average_intensity}}
#'   \item{\code{reference_average_intensity}}{double, average intensity among the 
#'       replicates in the group of sex-matched sedentary control animals}
#'   \item{\code{reference_average_intensity_se}}{double, standard error of \code{reference_average_intensity}}
#'   \item{\code{adj_p_value}}{double, adjusted p-value from \code{p_value} column. 
#'       P-values are BY-adjusted across all datasets within a given \code{assay}.}
#'   \item{\code{tissue_abbreviation}}{character, MoTrPAC tissue abbreviation used in manuscripts}
#'   \item{\code{selection_fdr}}{double, adjusted training p-value used to 
#'       select training-regulated analytes. P-values are IHW-adjusted across all 
#'       datasets within a given \code{assay} with \code{tissue} as a covariate.} 
#'}
#' @name TRNSCRPT_DA
"TRNSCRPT_BLOOD_DA"

#' @rdname TRNSCRPT_DA
"TRNSCRPT_HIPPOC_DA"

#' @rdname TRNSCRPT_DA
"TRNSCRPT_CORTEX_DA"

#' @rdname TRNSCRPT_DA
"TRNSCRPT_HYPOTH_DA"

#' @rdname TRNSCRPT_DA
"TRNSCRPT_SKMGN_DA"

#' @rdname TRNSCRPT_DA
"TRNSCRPT_SKMVL_DA"

#' @rdname TRNSCRPT_DA
"TRNSCRPT_HEART_DA"

#' @rdname TRNSCRPT_DA
"TRNSCRPT_KIDNEY_DA"

#' @rdname TRNSCRPT_DA
"TRNSCRPT_ADRNL_DA"

#' @rdname TRNSCRPT_DA
"TRNSCRPT_COLON_DA"

#' @rdname TRNSCRPT_DA
"TRNSCRPT_SPLEEN_DA"

#' @rdname TRNSCRPT_DA
"TRNSCRPT_TESTES_DA"

#' @rdname TRNSCRPT_DA
"TRNSCRPT_OVARY_DA"

#' @rdname TRNSCRPT_DA
"TRNSCRPT_VENACV_DA"

#' @rdname TRNSCRPT_DA
"TRNSCRPT_LUNG_DA"

#' @rdname TRNSCRPT_DA
"TRNSCRPT_SMLINT_DA"

#' @rdname TRNSCRPT_DA
"TRNSCRPT_LIVER_DA"

#' @rdname TRNSCRPT_DA
"TRNSCRPT_BAT_DA"

#' @rdname TRNSCRPT_DA
"TRNSCRPT_WATSC_DA"


#' @title Differential analysis of proteomics datasets 
#' @description Timewise summary statistics and training FDR from 
#'     differential analysis (DA) that tests the effect of training on each 
#'     proteomics feature (protein, phosphosite, acetylsite, ubiquitlsite) 
#'     within each sex. One data frame per data type per tissue.  
#' @format A data frame with >40000 rows and 16 variables:
#' \describe{
#'   \item{\code{feature_ID}}{character, analyte identifier}
#'   \item{\code{tissue}}{character, MoTrPAC tissue release code}
#'   \item{\code{assay}}{character, MoTrPAC assay or "ome" code. "prot-pr" for 
#'       global proteomics; "prot-ph" for phosphoproteomics; "prot-ac" for
#'       acetylomics; "prot-ub" for ubiquitylomics.}
#'   \item{\code{sex}}{character, one of "male" or "female"}
#'   \item{\code{logFC}}{double, log fold-change where the numerator is \code{comparison_group} 
#'       and the denominator is the group of sex-matched sedentary control animals}
#'   \item{\code{logFC_se}}{double, standard error of the log fold-change}
#'   \item{\code{tscore}}{double, t statistic}
#'   \item{\code{covariates}}{character, comma-separated list of covariates or NA}
#'   \item{\code{comparison_group}}{character, time point compared to the sex-matched 
#'       sedentary control animals, one of "1w", "2w", "4w", "8w"}
#'   \item{\code{p_value}}{double, unadjusted p-value for the difference between 
#'       \code{comparison_group} and the group of sex-matched sedentary control animals}
#'   \item{\code{comparison_average_intensity}}{double, average intensity among the replicates in \code{comparison_group}}
#'   \item{\code{reference_average_intensity}}{double, average intensity among the 
#'       replicates in the group of sex-matched sedentary control animals}
#'   \item{\code{numNAs}}{integer, number of missing values in the current 
#'       \code{comparison_group} and sex-matched sedentary control animals}
#'   \item{\code{adj_p_value}}{double, adjusted p-value from \code{p_value} column. 
#'       P-values are BY-adjusted across all datasets within a given \code{assay}.}
#'   \item{\code{tissue_abbreviation}}{character, MoTrPAC tissue abbreviation used in manuscripts}
#'   \item{\code{selection_fdr}}{double, adjusted training p-value used to 
#'       select training-regulated analytes. P-values are IHW-adjusted across all 
#'       datasets within a given \code{assay} with \code{tissue} as a covariate.} 
#'}
#' @name PROTEOME_DA
"PROT_CORTEX_DA"

#' @rdname PROTEOME_DA
"PROT_SKMGN_DA"

#' @rdname PROTEOME_DA
"PROT_HEART_DA"

#' @rdname PROTEOME_DA
"PROT_KIDNEY_DA"

#' @rdname PROTEOME_DA
"PROT_LUNG_DA"

#' @rdname PROTEOME_DA
"PROT_LIVER_DA"

#' @rdname PROTEOME_DA
"PROT_WATSC_DA"

#' @rdname PROTEOME_DA
"PHOSPHO_CORTEX_DA"

#' @rdname PROTEOME_DA
"PHOSPHO_SKMGN_DA"

#' @rdname PROTEOME_DA
"PHOSPHO_HEART_DA"

#' @rdname PROTEOME_DA
"PHOSPHO_KIDNEY_DA"

#' @rdname PROTEOME_DA
"PHOSPHO_LUNG_DA"

#' @rdname PROTEOME_DA
"PHOSPHO_LIVER_DA"

#' @rdname PROTEOME_DA
"PHOSPHO_WATSC_DA"

#' @rdname PROTEOME_DA
"ACETYL_HEART_DA"

#' @rdname PROTEOME_DA
"ACETYL_LIVER_DA"

#' @rdname PROTEOME_DA
"UBIQ_HEART_DA"

#' @rdname PROTEOME_DA
"UBIQ_LIVER_DA"


#' @title Differential analysis of merged metabolomics datasets 
#' @description Timewise summary statistics and training FDR from 
#'     differential analysis (DA) that tests the effect of training on each 
#'     metabolite within each sex. One data frame per tissue.  
#' @format A data frame with >800 rows and 25 variables:
#' \describe{
#'   \item{\code{feature_ID}}{character, analyte identifier}
#'   \item{\code{tissue}}{character, MoTrPAC tissue release code}
#'   \item{\code{dataset}}{character, mame of platform used by the site, e.g. "metab-u-hilicpos"}
#'   \item{\code{assay}}{character, MoTrPAC assay or "ome" code. "metab" for metabolomics datasets.}
#'   \item{\code{sex}}{character, one of "male" or "female"}
#'   \item{\code{is_targeted}}{logical, is this a targeted platform?}
#'   \item{\code{site}}{character, Chemical Analysis Site (CAS) name}
#'   \item{\code{metabolite}}{character, name of metabolite as appears in the CAS's data}
#'   \item{\code{metabolite_refmet}}{character, RefMet name of metabolite OR the site-given code for unnamed metabolites}
#'   \item{\code{logFC}}{double, log fold-change where the numerator is \code{comparison_group} 
#'       and the denominator is the group of sex-matched sedentary control animals}
#'   \item{\code{logFC_se}}{double, standard error of the log fold-change}
#'   \item{\code{tscore}}{double, t statistic}
#'   \item{\code{cv}}{double, feature coefficient of variation in the dataset}
#'   \item{\code{control_cv}}{double, feature coffeicient of variation in the control samples}
#'   \item{\code{covariates}}{character, comma-separated list of covariates or NA}
#'   \item{\code{comparison_group}}{character, time point compared to the sex-matched 
#'       sedentary control animals, one of "1w", "2w", "4w", "8w"}
#'   \item{\code{p_value}}{double, unadjusted p-value for the difference between 
#'       \code{comparison_group} and the group of sex-matched sedentary control animals}
#'   \item{\code{comparison_average_intensity}}{double, average intensity among the replicates in \code{comparison_group}}
#'   \item{\code{reference_average_intensity}}{double, average intensity among the 
#'       replicates in the group of sex-matched sedentary control animals}
#'   \item{\code{tissue_abbreviation}}{character, MoTrPAC tissue abbreviation used in manuscripts}
#'   \item{\code{mz}}{double, mass over charge}
#'   \item{\code{rt}}{double, retention time}
#'   \item{\code{neutral_mass}}{double, neutral mass}
#'   \item{\code{adj_p_value}}{double, adjusted p-value from \code{p_value} column. 
#'       P-values are BY-adjusted across all datasets within a given \code{assay}.}
#'   \item{\code{selection_fdr}}{double, adjusted training p-value used to 
#'       select training-regulated analytes. P-values are IHW-adjusted across all 
#'       datasets within a given \code{assay} with \code{tissue} as a covariate.}  
#'}
#' @name METAB_DA
"METAB_PLASMA_DA"

#' @rdname METAB_DA
"METAB_HIPPOC_DA"

#' @rdname METAB_DA
"METAB_CORTEX_DA"

#' @rdname METAB_DA
"METAB_HYPOTH_DA"

#' @rdname METAB_DA
"METAB_SKMGN_DA"

#' @rdname METAB_DA
"METAB_SKMVL_DA"

#' @rdname METAB_DA
"METAB_HEART_DA"

#' @rdname METAB_DA
"METAB_KIDNEY_DA"

#' @rdname METAB_DA
"METAB_ADRNL_DA"

#' @rdname METAB_DA
"METAB_COLON_DA"

#' @rdname METAB_DA
"METAB_SPLEEN_DA"

#' @rdname METAB_DA
"METAB_TESTES_DA"

#' @rdname METAB_DA
"METAB_OVARY_DA"

#' @rdname METAB_DA
"METAB_VENACV_DA"

#' @rdname METAB_DA
"METAB_LUNG_DA"

#' @rdname METAB_DA
"METAB_SMLINT_DA"

#' @rdname METAB_DA
"METAB_LIVER_DA"

#' @rdname METAB_DA
"METAB_BAT_DA"

#' @rdname METAB_DA
"METAB_WATSC_DA"


#' @title Differential analysis of multiplexed immunoassays 
#' @description Timewise summary statistics and training FDR from 
#'     differential analysis (DA) that tests the effect of training on each 
#'     immunoassay analyte within each sex. One data frame per tissue. 
#' @format A data frame with >180 rows and 15 variables:
#' \describe{
#'   \item{\code{feature_ID}}{character, analyte identifier}
#'   \item{\code{panel}}{character, name of LUMINEX panel, e.g., "rat-mag27plex"}
#'   \item{\code{tissue}}{character, MoTrPAC tissue release code}
#'   \item{\code{assay}}{character, MoTrPAC assay or "ome" code. "immunoassay" for multiplexed immunoassay datasets.}
#'   \item{\code{comparison_group}}{character, time point compared to the sex-matched 
#'       sedentary control animals, one of "1w", "2w", "4w", "8w"}
#'   \item{\code{logFC}}{double, log fold-change where the numerator is \code{comparison_group} 
#'       and the denominator is the group of sex-matched sedentary control animals}
#'   \item{\code{logFC_se}}{double, standard error of the log fold-change}
#'   \item{\code{p_value}}{double, unadjusted p-value for the difference between 
#'       \code{comparison_group} and the group of sex-matched sedentary control animals}
#'   \item{\code{reference_average_intensity_se}}{double, standard error of \code{reference_average_intensity}}
#'   \item{\code{comparison_average_intensity}}{double, average intensity among the replicates in \code{comparison_group}}
#'   \item{\code{sex}}{character, one of "male" or "female"}
#'   \item{\code{covariates}}{character, comma-separated list of covariates or NA}
#'   \item{\code{adj_p_value}}{double, adjusted p-value from \code{p_value} column. 
#'       P-values are BY-adjusted across all datasets within a given \code{assay}.}
#'   \item{\code{tissue_abbreviation}}{character, MoTrPAC tissue abbreviation used in manuscripts}
#'   \item{\code{selection_fdr}}{double, adjusted training p-value used to 
#'       select training-regulated analytes. P-values are IHW-adjusted across all 
#'       datasets within a given \code{assay} with \code{tissue} as a covariate.} 
#'}
#' @name IMMUNO_DA
"IMMUNO_PLASMA_DA"

#' @rdname IMMUNO_DA
"IMMUNO_HIPPOC_DA"

#' @rdname IMMUNO_DA
"IMMUNO_CORTEX_DA"

#' @rdname IMMUNO_DA
"IMMUNO_ADRNL_DA"

#' @rdname IMMUNO_DA
"IMMUNO_TESTES_DA"

#' @rdname IMMUNO_DA
"IMMUNO_OVARY_DA"

#' @rdname IMMUNO_DA
"IMMUNO_SKMGN_DA"

#' @rdname IMMUNO_DA
"IMMUNO_BAT_DA"

#' @rdname IMMUNO_DA
"IMMUNO_WATSC_DA"

#' @rdname IMMUNO_DA
"IMMUNO_COLON_DA"

#' @rdname IMMUNO_DA
"IMMUNO_SMLINT_DA"

#' @rdname IMMUNO_DA
"IMMUNO_LIVER_DA"

#' @rdname IMMUNO_DA
"IMMUNO_SKMVL_DA"

#' @rdname IMMUNO_DA
"IMMUNO_HEART_DA"

#' @rdname IMMUNO_DA
"IMMUNO_KIDNEY_DA"

#' @rdname IMMUNO_DA
"IMMUNO_SPLEEN_DA"

#' @rdname IMMUNO_DA
"IMMUNO_LUNG_DA"
