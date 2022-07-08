#' @title Tissue abbreviations 
#' @description Tissue abbreviations used in tables and figures 
#' @format Unnamed character vector
"TISSUE_ABBREV"

#' @title Tissue colors 
#' @description Tissue colors used in figures. Values are hex codes and names are tissue abbreviations (\code{TISSUE_ABBREV}) OR tissue codes. 
#' @format Named character vector 
"TISSUE_COLORS"

#' @title Assay colors 
#' @description Assay colors used in figures. Values are hex codes and names are assay abbreviations (\code{ASSAY_ABBREV}). 
#' @format Named character vector 
"ASSAY_COLORS"

#' @title Assay or "ome" abbreviations 
#' @description Assay/ome abbreviations used in tables and figures 
#' @format Unnamed character vector
"ASSAY_ABBREV"

#' @title Assay code-to-abbreviation mapping 
#' @description Values are abbreviations (\code{ASSAY_ABBREV}) and names are codes. 
#'     See [MotrpacBicQC::assay_codes] for more details.  
#' @format Named character vector 
"ASSAY_CODE_TO_ABBREV"

#' @title Assay abbreviation-to-code mapping 
#' @description Values are codes and names are abbreviations (\code{ASSAY_ABBREV}).  
#'     See [MotrpacBicQC::assay_codes] for more details.
#' @format Named character vector 
"ASSAY_ABBREV_TO_CODE"

#' @title Tissue code-to-abbreviation mapping 
#' @description Values are abbreviations (\code{TISSUE_ABBREV}) and names are codes. 
#'     See [MotrpacBicQC::bic_animal_tissue_code] for more details.  
#' @format Named character vector 
"TISSUE_CODE_TO_ABBREV"

#' @title Tissue abbreviation-to-code mapping 
#' @description Values are codes and names are abbreviations (\code{TISSUE_ABBREV}).  
#'     See [MotrpacBicQC::bic_animal_tissue_code] for more details.
#' @format Named character vector 
"TISSUE_ABBREV_TO_CODE"

#' @title Intervention group colors 
#' @description Intervention group colors used in figures. Values are hex codes and names are intervention groups.  
#' @format Named character vector
"GROUP_COLORS"

#' @title Sex colors 
#' @description Sex colors used in figures. Values are hex codes and names include multiple ways to indicate male or female.  
#' @format Named character vector
"SEX_COLORS"

#' @title Assay order 
#' @description Biological order of \code{ASSAY_ABBREV} used in figures  
#' @format Unnamed character vector
"ASSAY_ORDER"

#' @title Tissue order 
#' @description Biological order of \code{TISSUE_ABBREV} used in figures  
#' @format Unnamed character vector
"TISSUE_ORDER"


#' @name FEATURE_TO_GENE
#' @title Feature-to-gene map
#' @description This feature-to-gene map associates every feature tested in 
#'     differential analysis with a gene and includes all current gene identifiers available in RGD as of 11/12/2020. 
#' @format A data frame with 4044034 rows and 9 variables:
#' \describe{
#'   \item{\code{entrez_gene}}{double, Entrez gene ID}
#'   \item{\code{feature_ID}}{character, MoTrPAC feature identifier}
#'   \item{\code{rgd_gene}}{integer, RGD gene ID}
#'   \item{\code{gene_symbol}}{character, official gene symbol}
#'   \item{\code{old_gene_symbol}}{character, semicolon-separated list of deprecated or alias gene symbols}
#'   \item{\code{ensembl_gene}}{character, Ensembl gene ID}
#'   \item{\code{relationship_to_gene}}{double, for ATAC and METHYL features only. Distance 
#'       from the closest edge of the feature to the start or end of the closest gene, whichever is closer.
#'       A value of 0 means there is non-zero overlap between the feature and the gene.
#'       A negative value means the feature is upstream of "geneStart".
#'       A a positive value means the feature is downstream of "geneEnd".
#'       Note that "geneStart" and "geneEnd" are strand-agnostic, i.e. "geneStart" 
#'       is always less than "geneEnd", even if the gene is on the negative strand ("geneStrand" == 2).}
#'   \item{\code{custom_annotation}}{character, a version of the \code{ChIPseeker} annotations with many corrections. Values include:
#'       "Distal Intergenic", "Promoter (<=1kb)", "Exon", "Promoter (1-2kb)", "Downstream (<5kb)", "Upstream (<5kb)", "5' UTR", "Intron", "3' UTR", "Overlaps Gene"}
#'   \item{\code{kegg_id}}{character, KEGG ID for METAB features only. See [MotrpacBicQC::get_and_validate_mdd()] for more details.} 
#'}
#' @details All proteomics feature IDs (RefSeq accessions) were mapped to gene 
#'     symbols and Entrez IDs using NCBI’s "gene2refseq" mapping files 
#'     (<https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2refseq.gz>, downloaded on 2020/12/18). 
#'     Epigenomics features were mapped to the nearest gene using the R function [ChIPseeker::annotatePeak()] 
#'     with Ensembl gene annotation (Rattus norvegicus release 96). 
#'     Gene symbols, Entrez IDs, Ensembl IDs, and RGD IDs were mapped to each 
#'     other using RGD’s rat gene annotation 
#'     (<https://download.rgd.mcw.edu/data_release/RAT/GENES_RAT.txt>, generated on 2021/11/12).
#' 
#'     For fast(er) indexing, convert this object to a \code{data.table} and use 
#'     [data.table::setkey()] to set the key to the column you are matching. 
#'     This dramatically improves performance. 
#'     
"FEATURE_TO_GENE"


#' @name REPEATED_FEATURES
#' @title Repeated feature info
#' @description Documentation of 91 IMMUNO and METAB features that were measured on multiple platforms 
#'     and have multiple differential analysis results. New unique identifiers are defined. These new 
#'     unique identifiers were used in the graphical analysis. 
#' @format A data frame with 825 rows and 10 variables:
#' \describe{
#'   \item{\code{feature_ID}}{character, MoTrPAC feature identifier}
#'   \item{\code{assay_abbr}}{character, assay abbreviation}
#'   \item{\code{tissue_abbreviation}}{character, tissue abbreviation}
#'   \item{\code{dataset}}{character, platform for metabolomics features only. 
#'       "meta-reg" indicates the feature came from the merged meta-regression results.
#'       See more details in [MotrpacBicQC::assay_codes].}
#'   \item{\code{panel}}{character, name of LUMINEX panel for immunoassay features only}
#'   \item{\code{selection_fdr}}{double, adjusted training p-value used to 
#'       select training-regulated analytes. P-values are IHW-adjusted across all 
#'       datasets within a given \code{assay} with \code{tissue} as a covariate.}
#'   \item{\code{metabolite_refmet}}{character, RefMet name of metabolite OR the site-given code for unnamed metabolites}
#'   \item{\code{feature}}{character, duplicated \code{feature} in the format \code{[ASSAY_ABBREV];[TISSUE_ABBREV];[feature_ID]}}
#'   \item{\code{new_feature_ID}}{character, new unique \code{feature_ID} with \code{panel} or \code{dataset} prepended}
#'   \item{\code{new_feature}}{character, new unique \code{feature} in the format \code{[ASSAY_ABBREV];[TISSUE_ABBREV];[new_feature_ID]}} 
#'}
#' @details 91 IMMUNO and METAB features were measured on multiple platforms and have multiple differential analysis results. 
#'     In order to distinguish between the different sets of results, these 
#'     feature_IDs were prepended with \code{panel} (for IMMUNO) or \code{dataset} (for METAB). 
#'     For example, "BDNF" is separated into "rat-myokine:BDNF" and "rat-pituitary:BDNF". 
#'     For completeness, both the modified and unmodified feature_IDs are included in both 
#'     the feature-to-gene map and universe lists. Note that the graphical analysis
#'     uses the modified feature_IDs; differential analysis results use the unmodified feature_IDs.
"REPEATED_FEATURES"


#' @title Phenotypic data 
#' @description Phenotypic data for samples from the MoTrPAC endurance exercise 
#'     training study in 6-month-old rats. One row per sample (\code{viallabel}).
#' @format A data frame with 5955 rows and 509 variables:
#' \describe{
#'   \item{\code{pid}}{integer, unique, randomly generated 8 digit numeric identifier used in linkage to phenotypic data}
#'   \item{\code{bid}}{integer, unique, randomly generated 8 digit numeric identifier used in linkage to phenotypic data}
#'   \item{\code{labelid}}{double, unique 11 digit identifier for specimen label ID, originating at the collection site, 
#'       that provides a link to specimen processing and used for shipments to the biorepository (same as 
#'       \code{vialLabel} only in instances where aliquots were not further processed at the biorepository)}
#'   \item{\code{viallabel}}{character, unique 11 digit numeric identifier of a sample vial made up of \code{bid}, 
#'       \code{key.sacrificetime}, \code{specimen.processing.sampletypedescription}, and aliquot number}
#'   \item{\code{key.protocol}}{character, protocol in which the procedures for sample collection were outlined}
#'   \item{\code{key.agegroup}}{character, age group for the animal}
#'   \item{\code{key.d_arrive}}{character, date animal arrived at site (same as \code{registration.d_arrive})}
#'   \item{\code{key.d_sacrifice}}{character, date animal was sacrificed (same as \code{specimen.collection.d_visit})}
#'   \item{\code{key.intervention}}{character, exercise group defined by those who complete exercise training}
#'   \item{\code{key.sacrificetime}}{character, weeks since the beginning of the intervention}
#'   \item{\code{key.anirandgroup}}{character, description of the randomization group}
#'   \item{\code{key.siteid}}{character, site where the samples were collected}
#'   \item{\code{registration.d_visit}}{character, date registration was completed}
#'   \item{\code{registration.days_visit}}{integer, day registration was completed relative to \code{registration.d_arrive} (count)}
#'   \item{\code{registration.staffid}}{integer, registration staff ID}
#'   \item{\code{registration.siteid}}{character, registration site}
#'   \item{\code{registration.formname}}{character, registration form name}
#'   \item{\code{registration.versionnbr}}{integer, registration form version number}
#'   \item{\code{registration.ratid}}{character, site rat ID}
#'   \item{\code{registration.d_arrive}}{character, date the rat arrived at the site}
#'   \item{\code{registration.days_arrive}}{integer, day of arrival relative to \code{registration.d_arrive} (count)}
#'   \item{\code{registration.batchnumber}}{integer, registration batch number}
#'   \item{\code{registration.d_reverselight}}{character, date the rat was housed in reverse light}
#'   \item{\code{registration.days_reverselight}}{integer, day the rat was housed in reverse light relative to \code{registration.d_arrive} (count)}
#'   \item{\code{registration.d_birth}}{character, rat date of birth (+/- two weeks)}
#'   \item{\code{registration.days_birth}}{integer, rat day of birth relative to \code{registration.d_arrive} (count)}
#'   \item{\code{registration.sex}}{character, sex, one of "male", "female"}
#'   \item{\code{registration.weight}}{double, rat's body weight}
#'   \item{\code{registration.cagenumber}}{integer, cage number in which the rat was housed}
#'   \item{\code{familiarization.d_visit}}{character, date familiarization was completed}
#'   \item{\code{familiarization.days_visit}}{integer, day familiarization was completed relative to \code{registration.d_arrive} (count)}
#'   \item{\code{familiarization.staffid}}{integer, familiarization staff ID}
#'   \item{\code{familiarization.siteid}}{character, familiarization site}
#'   \item{\code{familiarization.formname}}{character, familiarization form name}
#'   \item{\code{familiarization.versionnbr}}{integer, familiarization form version number}
#'   \item{\code{familiarization.d_treadmillbegin}}{character, date the rat began treadmill familiarization}
#'   \item{\code{familiarization.days_treadmillbegin}}{integer, day the rat began treadmill familiarization relative to \code{registration.d_arrive} (count)}
#'   \item{\code{familiarization.d_treadmillcomplete}}{character, date the rat completed treadmill familiarization}
#'   \item{\code{familiarization.days_treadmillcomplete}}{integer, day the rat completed treadmill familiarization relative to \code{registration.d_arrive} (count)}
#'   \item{\code{familiarization.activity_score}}{double, rat's activity score}
#'   \item{\code{familiarization.compliant}}{character, was this rat compliant?}
#'   \item{\code{familiarization.weight}}{double, rat's body weight}
#'   \item{\code{training.d_visit}}{character, date intervention started}
#'   \item{\code{training.days_visit}}{integer, day intervention started relative to \code{registration.d_arrive}}
#'   \item{\code{training.staffid}}{integer, training staff ID}
#'   \item{\code{training.siteid}}{character, training site}
#'   \item{\code{training.formname}}{character, training form name}
#'   \item{\code{training.versionnbr}}{integer, training form version number}
#'   \item{\code{training.day*date}}{character, date of given training day}
#'   \item{\code{training.day*_days}}{integer, day of given training day relative to \code{registration.d_arrive} (count)}
#'   \item{\code{training.day*time}}{character, time of given training day}
#'   \item{\code{training.day*_treadmillspeed}}{double, treadmill speed on given training day}
#'   \item{\code{training.day*_treadmillincline}}{double, treadmill incline on given training day}
#'   \item{\code{training.day*_timeontreadmill}}{integer, time on treadmill on given training day}
#'   \item{\code{training.day*_weight}}{double, weight on given training day}
#'   \item{\code{training.day*_posttrainlact}}{double, post-training lactate on given training day}
#'   \item{\code{training.day*_score}}{integer, score on given training day}
#'   \item{\code{training.day*_comments}}{character, comments on given training day}
#'   \item{\code{vo2.max.test.d_visit_1}}{character, date first VO2 max test was completed}
#'   \item{\code{vo2.max.test.d_visit_2}}{character, date second VO2 max test was completed}
#'   \item{\code{vo2.max.test.days_visit_1}}{integer, day first VO2 max test was completed relative to \code{registration.d_arrive} (count)}
#'   \item{\code{vo2.max.test.days_visit_2}}{integer, day second VO2 max test was completed relative to \code{registration.d_arrive} (count)}
#'   \item{\code{vo2.max.test.staffid_1}}{integer, first VO2 max test staff ID}
#'   \item{\code{vo2.max.test.staffid_2}}{integer, second VO2 max test staff ID}
#'   \item{\code{vo2.max.test.siteid_1}}{character, first VO2 max test site}
#'   \item{\code{vo2.max.test.siteid_2}}{character, second VO2 max test site}
#'   \item{\code{vo2.max.test.visitcode_1}}{character, visit code for first VO2 max test}
#'   \item{\code{vo2.max.test.visitcode_2}}{character, visit code for second VO2 max test}
#'   \item{\code{vo2.max.test.formname_1}}{character, form name for first VO2 max test}
#'   \item{\code{vo2.max.test.formname_2}}{character, form name for second VO2 max test}
#'   \item{\code{vo2.max.test.versionnbr_1}}{integer, form version number for first VO2 max test}
#'   \item{\code{vo2.max.test.versionnbr_2}}{integer, form version number for second VO2 max test}
#'   \item{\code{vo2.max.test.d_vo2_1}}{character, date of first VO2 max test}
#'   \item{\code{vo2.max.test.d_vo2_2}}{character, date of second VO2 max test}
#'   \item{\code{vo2.max.test.days_vo2_1}}{integer, date of first VO2 max test relative to \code{registration.d_arrive} (count)}
#'   \item{\code{vo2.max.test.days_vo2_2}}{integer, date of second VO2 max test relative to \code{registration.d_arrive} (count)}
#'   \item{\code{vo2.max.test.t_vo2_1}}{character, start time of first VO2 max test?}
#'   \item{\code{vo2.max.test.t_vo2_2}}{character, start time of second VO2 max test?}
#'   \item{\code{vo2.max.test.blactate_begin_1}}{double, blood lactate at the beginning of exercise for first VO2 max test}
#'   \item{\code{vo2.max.test.blactate_begin_2}}{double, blood lactate at the beginning of exercise for second VO2 max test}
#'   \item{\code{vo2.max.test.vo2_max_1}}{double, VO2 at max for first VO2 max test}
#'   \item{\code{vo2.max.test.vo2_max_2}}{double, VO2 at max for second VO2 max test}
#'   \item{\code{vo2.max.test.vco2_max_1}}{double, VCO2 at max for first VO2 max test}
#'   \item{\code{vo2.max.test.vco2_max_2}}{double, VCO2 at max for second VO2 max test}
#'   \item{\code{vo2.max.test.rer_max_1}}{double, RER at max for first VO2 max test}
#'   \item{\code{vo2.max.test.rer_max_2}}{double, RER at max for second VO2 max test}
#'   \item{\code{vo2.max.test.speed_max_1}}{double, max speed of running for first VO2 max test}
#'   \item{\code{vo2.max.test.speed_max_2}}{double, max speed of running for second VO2 max test}
#'   \item{\code{vo2.max.test.t_complete_1}}{character, time of completion for first VO2 max test}
#'   \item{\code{vo2.max.test.t_complete_2}}{character, time of completion for second VO2 max test}
#'   \item{\code{vo2.max.test.blactate_end_1}}{double, blood lactate at end of exercise for first VO2 max test}
#'   \item{\code{vo2.max.test.blactate_end_2}}{double, blood lactate at end of exercise for second VO2 max test}
#'   \item{\code{vo2.max.test.vo2_comments_1}}{character, comments about first VO2 max test}
#'   \item{\code{vo2.max.test.vo2_comments_2}}{character, comments about second VO2 max test}
#'   \item{\code{nmr.testing.d_visit_1}}{character, date first NMR test was completed}
#'   \item{\code{nmr.testing.d_visit_2}}{character, date second NMR test was completed}
#'   \item{\code{nmr.testing.days_visit_1}}{integer, day first NMR test was completed relative to \code{registration.d_arrive} (count)}
#'   \item{\code{nmr.testing.days_visit_2}}{integer, day second NMR test was completed relative to \code{registration.d_arrive} (count)}
#'   \item{\code{nmr.testing.staffid_1}}{integer, first NMR test staff ID}
#'   \item{\code{nmr.testing.staffid_2}}{integer, second NMR test staff ID}
#'   \item{\code{nmr.testing.siteid_1}}{character, first NMR test site}
#'   \item{\code{nmr.testing.siteid_2}}{character, second NMR test site}
#'   \item{\code{nmr.testing.visitcode_1}}{character, visit code for first NMR test}
#'   \item{\code{nmr.testing.visitcode_2}}{character, visit code for second NMR test}
#'   \item{\code{nmr.testing.formname_1}}{character, form name for first NMR test}
#'   \item{\code{nmr.testing.formname_2}}{character, form name for second NMR test}
#'   \item{\code{nmr.testing.versionnbr_1}}{integer, form version number for first NMR test}
#'   \item{\code{nmr.testing.versionnbr_2}}{integer, form version number for second NMR test}
#'   \item{\code{nmr.testing.d_nmr_1}}{character, date data was collected for the first NMR test}
#'   \item{\code{nmr.testing.d_nmr_2}}{character, date data was collected for the second NMR test}
#'   \item{\code{nmr.testing.days_nmr_1}}{integer, day data was collected for the first NMR test relative to \code{registration.d_arrive} (count)}
#'   \item{\code{nmr.testing.days_nmr_2}}{integer, day data was collected for the second NMR test relative to \code{registration.d_arrive} (count)}
#'   \item{\code{nmr.testing.nmr_weight_1}}{double, rat's body weight for first NMR test}
#'   \item{\code{nmr.testing.nmr_weight_2}}{double, rat's body weight for second NMR test}
#'   \item{\code{nmr.testing.nmr_fat_1}}{double, rat's percent body fat for first NMR test}
#'   \item{\code{nmr.testing.nmr_fat_2}}{double, rat's percent body fat for second NMR test}
#'   \item{\code{nmr.testing.nmr_lean_1}}{double, rat's percent lean for first NMR test}
#'   \item{\code{nmr.testing.nmr_lean_2}}{double, rat's percent lean for second NMR test}
#'   \item{\code{nmr.testing.nmr_fluid_1}}{double, rat's percent fluid for first NMR test}
#'   \item{\code{nmr.testing.nmr_fluid_2}}{double, rat's percent fluid for second NMR test}
#'   \item{\code{nmr.testing.nmr_comments_2}}{character, comments about second NMR test}
#'   \item{\code{terminal.weight.bw}}{double, terminal body weight (g)}
#'   \item{\code{terminal.weight.lg}}{integer, terminal lateral gastrocnemius weight (mg)}
#'   \item{\code{terminal.weight.mg}}{integer, terminal medial gastrocnemius weight (mg)}
#'   \item{\code{terminal.weight.pl}}{integer, terminal plantaris weight (mg)}
#'   \item{\code{terminal.weight.sol}}{integer, terminal soleus weight (mg)}
#'   \item{\code{specimen.collection.d_visit}}{character, date of specimen collection}
#'   \item{\code{specimen.collection.days_visit}}{integer, day of specimen collection relative to \code{registration.d_arrive} (count)}
#'   \item{\code{specimen.collection.staffid}}{integer, specimen collection staff ID}
#'   \item{\code{specimen.collection.siteid}}{character, specimen collection site}
#'   \item{\code{specimen.collection.visitcode}}{character, specimen collection site visit code}
#'   \item{\code{specimen.collection.formname}}{character, specimen collection form name}
#'   \item{\code{specimen.collection.versionnbr}}{integer, specimen collection form version number}
#'   \item{\code{specimen.collection.anesthesiaid}}{integer, anesthesia administrator ID}
#'   \item{\code{specimen.collection.t_anesthesia}}{character, anesthesia administration time}
#'   \item{\code{specimen.collection.bloodtype}}{character, blood type}
#'   \item{\code{specimen.collection.bloodtube}}{character, blood tube}
#'   \item{\code{specimen.collection.bloodcomplete}}{character, blood collection complete}
#'   \item{\code{specimen.collection.bloodtechid}}{character, blood collection technician ID}
#'   \item{\code{specimen.collection.t_bloodstart}}{character, blood collection start time}
#'   \item{\code{specimen.collection.t_bloodstop}}{character, blood collection stop time}
#'   \item{\code{specimen.collection.t_edtafill}}{character, EDTA fill time}
#'   \item{\code{specimen.collection.uterustype}}{character, uterus type}
#'   \item{\code{specimen.collection.uterustechid}}{character, uterus technician ID}
#'   \item{\code{specimen.collection.uteruscomplete}}{character, uterus complete}
#'   \item{\code{specimen.collection.t_death}}{character, time of death}
#'   \item{\code{specimen.collection.deathtype}}{character, type of death}
#'   \item{\code{specimen.processing.versionnbr}}{integer, specimen processing form version number}
#'   \item{\code{specimen.processing.formname}}{character, specimen processing form name}
#'   \item{\code{specimen.processing.visitcode}}{character, specimen processing visit code}
#'   \item{\code{specimen.processing.siteid}}{character, specimen processing site}
#'   \item{\code{specimen.processing.sampletypedescription}}{character, specimen tissue type}
#'   \item{\code{specimen.processing.samplenumber}}{integer, specimen sample number}
#'   \item{\code{specimen.processing.timepoint}}{integer, specimen time point}
#'   \item{\code{specimen.processing.aliquotdescription}}{character, specimen aliquot tube}
#'   \item{\code{specimen.processing.volume}}{character, specimen aliquot volume}
#'   \item{\code{specimen.processing.hemolyzed}}{character, specimen hemolyzed?}
#'   \item{\code{specimen.processing.t_collection}}{character, specimen collection time}
#'   \item{\code{specimen.processing.t_edtaspin}}{character, EDTA spin time}
#'   \item{\code{specimen.processing.t_freeze}}{character, specimen freeze time}
#'   \item{\code{specimen.processing.techid}}{character, specimen processing tech ID}
#'   \item{\code{calculated.variables.wgt_gain_after_train}}{double, weight gain after training: 
#'       \code{max(nmr.testing.nmr_weight_2, nmr.testing.nmr_weight_1) - registration.weight}}
#'   \item{\code{calculated.variables.pct_body_fat_change}}{double, percent body fat change after training: 
#'       \code{nmr.testing.nmr_fat_2 - nmr.testing.nmr_fat_1}}
#'   \item{\code{calculated.variables.pct_body_lean_change}}{double, percent body lean change after training:
#'       \code{nmr.testing.nmr_lean_2 - nmr.testing.nmr_lean_1}}
#'   \item{\code{calculated.variables.pct_body_fluid_change}}{double, percent body fluid change after training:
#'       \code{nmr.testing.nmr_fluid_2 - nmr.testing.nmr_fluid_1}}
#'   \item{\code{calculated.variables.lactate_change_dueto_train}}{double, lactate change due to training:
#'       \code{vo2.max.test.blactate_end_2 - vo2.max.test.blactate_begin_1}}
#'   \item{\code{calculated.variables.vo2_max_change}}{double, change in VO2 max after training:
#'       \code{vo2.max.test.vo2_max_2 - vo2.max.test.vo2_max_1}}
#'   \item{\code{calculated.variables.coll_time_train}}{integer, sample collection time after training (seconds): 
#'       specimen collection date and time - last time trained}
#'   \item{\code{calculated.variables.deathtime_after_train}}{integer, time of death after training (seconds):
#'       death date and time - last time trained}
#'   \item{\code{calculated.variables.frozetime_after_train}}{integer, time of sample freezing after training (seconds): 
#'       specimen collection date and freeze time - last time trained}
#'   \item{\code{biclabeldata.shiptositeid}}{character, Chemical Analysis Site (CAS), i.e. the omics site)}
#'   \item{\code{sacrificetime}}{character, sacrifice time point in weeks, one of "1w", "2w", "4w", "8w"}
#'   \item{\code{intervention}}{character, intervention, one of "training", "control"}
#'   \item{\code{group}}{character, intervention group, one of "1w", "2w", "4w", "8w", "control"}
#'   \item{\code{sex}}{character, sex, one of "male", "female"}
#'   \item{\code{time_to_freeze}}{integer, sample time to freeze, \code{calculated.variables.frozetime_after_train - calculated.variables.deathtime_after_train}} 
#'}
"PHENO"


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
