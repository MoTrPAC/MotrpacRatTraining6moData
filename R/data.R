## Colors and abbreviations ####

#' @title Tissue abbreviations 
#' @description Tissue abbreviations used in tables and figures 
#' @format Unnamed character vector
"TISSUE_ABBREV"

#' @title Tissue colors 
#' @description Tissue colors used in figures. Values are hex codes and names are tissue abbreviations [TISSUE_ABBREV] OR tissue codes. 
#' @format Named character vector 
"TISSUE_COLORS"

#' @title Assay colors 
#' @description Assay colors used in figures. Values are hex codes and names are assay abbreviations [ASSAY_ABBREV]. 
#' @format Named character vector 
"ASSAY_COLORS"

#' @title Assay or "ome" abbreviations 
#' @description Assay/ome abbreviations used in tables and figures 
#' @format Unnamed character vector
"ASSAY_ABBREV"

#' @title Assay code-to-abbreviation mapping 
#' @description Values are abbreviations [ASSAY_ABBREV] and names are codes. 
#'     See [MotrpacBicQC::assay_codes] for more details.  
#' @format Named character vector 
"ASSAY_CODE_TO_ABBREV"

#' @title Assay abbreviation-to-code mapping 
#' @description Values are codes and names are abbreviations [ASSAY_ABBREV].  
#'     See [MotrpacBicQC::assay_codes] for more details.
#' @format Named character vector 
"ASSAY_ABBREV_TO_CODE"

#' @title Tissue code-to-abbreviation mapping 
#' @description Values are abbreviations [TISSUE_ABBREV] and names are codes. 
#'     See [MotrpacBicQC::bic_animal_tissue_code] for more details.  
#' @format Named character vector 
"TISSUE_CODE_TO_ABBREV"

#' @title Tissue abbreviation-to-code mapping 
#' @description Values are codes and names are abbreviations [TISSUE_ABBREV].  
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
#' @description Biological order of [ASSAY_ABBREV] used in figures  
#' @format Unnamed character vector
"ASSAY_ORDER"

#' @title Tissue order 
#' @description Biological order of [TISSUE_ABBREV] used in figures  
#' @format Unnamed character vector
"TISSUE_ORDER"


## Gene mappings ####

#' @name FEATURE_TO_GENE
#' @title Feature-to-gene map
#' @description This feature-to-gene map associates every feature tested in 
#'     differential analysis with a gene and includes all current gene identifiers available in RGD as of 11/12/2020. 
#' @format A data frame with 4044034 rows and 9 variables:
#' \describe{
#'   \item{\code{entrez_gene}}{double, Entrez gene ID}
#'   \item{\code{feature_ID}}{`r feature_ID()`}
#'   \item{\code{rgd_gene}}{integer, RGD gene ID}
#'   \item{\code{gene_symbol}}{character, official gene symbol}
#'   \item{\code{old_gene_symbol}}{character, semicolon-separated list of deprecated or alias gene symbols}
#'   \item{\code{ensembl_gene}}{character, Ensembl gene ID}
#'   \item{\code{relationship_to_gene}}{double, for ATAC and METHYL features only. Distance 
#'     from the closest edge of the feature to the start or end of the closest gene, whichever is closer.
#'     A value of 0 means there is non-zero overlap between the feature and the gene.
#'     A negative value means the feature is upstream of "geneStart".
#'     A a positive value means the feature is downstream of "geneEnd".
#'     Note that "geneStart" and "geneEnd" are strand-agnostic, i.e. "geneStart" 
#'     is always less than "geneEnd", even if the gene is on the negative strand ("geneStrand" == 2).}
#'   \item{\code{custom_annotation}}{character, a version of the \code{ChIPseeker} annotations with many corrections. Values include:
#'     "Distal Intergenic", "Promoter (<=1kb)", "Exon", "Promoter (1-2kb)", "Downstream (<5kb)", "Upstream (<5kb)", "5' UTR", "Intron", "3' UTR", "Overlaps Gene"}
#'   \item{\code{kegg_id}}{character, KEGG ID for METAB features only. See [MotrpacBicQC::metabolomics_data_dictionary] for more details.} 
#' }
#' @details All proteomics feature IDs (RefSeq accessions) were mapped to gene 
#'     symbols and Entrez IDs using NCBI’s "gene2refseq" mapping files 
#'     (<https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2refseq.gz>, downloaded on 2020/12/18). 
#'     Epigenomics features were mapped to the nearest gene using the R function [ChIPseeker::annotatePeak()] 
#'     with Ensembl gene annotation (Rattus norvegicus release 96). 
#'     Gene symbols, Entrez IDs, Ensembl IDs, and RGD IDs were mapped to each 
#'     other using RGD’s rat gene annotation 
#'     (<https://download.rgd.mcw.edu/data_release/RAT/GENES_RAT.txt>, generated on 2021/11/12).
#' 
#'     For fast(er) indexing, convert this object to a [data.table::data.table()] and use 
#'     [data.table::setkey()] to set the key to the column you are matching. 
#'     This dramatically improves performance. 
#'     
"FEATURE_TO_GENE"


#' @title Rat-to-human gene ortholog map
#' @description Mapping between various rat and human gene identifiers for orthologs as reported by [RGD](https://rgd.mcw.edu/) (2020-01-10). 
#' @format A data frame with 21461 rows and 12 variables:
#' \describe{
#'   \item{\code{RAT_SYMBOL}}{character, official rat gene symbol}
#'   \item{\code{NAME}}{character, gene name}
#'   \item{\code{RAT_GENE_RGD_ID}}{integer, rat RGD gene ID}
#'   \item{\code{RAT_NCBI_GENE_ID}}{integer, rat NCBI/Entrez gene ID}
#'   \item{\code{RAT_ENSEMBL_ID}}{character, rat Ensembl gene ID}
#'   \item{\code{RAT_UNIPROT_ID}}{character, semicolon-separated list of rat UniProt gene IDs}
#'   \item{\code{RAT_OLD_SYMBOL}}{character, semicolon-separated list of rat deprecated or alias gene symbols}
#'   \item{\code{HUMAN_ORTHOLOG_SYMBOL}}{character, human official gene symbol}
#'   \item{\code{HUMAN_ORTHOLOG_RGD}}{character, human RGD gene ID}
#'   \item{\code{HUMAN_ORTHOLOG_NCBI_GENE_ID}}{character, human NCBI/Entrez gene ID}
#'   \item{\code{HUMAN_ORTHOLOG_ENSEMBL_ID}}{character, human Entrez gene ID}
#'   \item{\code{HUMAN_ORTHOLOG_SOURCE}}{character, human ortholog source} 
#' }
#' @details This map was compiled from several external sources. 
#'     GENCODE metadata and annotation files were used to map between human Ensembl transcript IDs, Entrez IDs, 
#'     GENCODE IDs, and Ensembl gene IDs [(Frankish et al., 2021)](https://pubmed.ncbi.nlm.nih.gov/33270111/). RGD files were used to map between human and 
#'     rat gene symbols as well as between various rat gene identifiers [(Smith et al., 2020)](https://pubmed.ncbi.nlm.nih.gov/31713623/). 
#' 
#'     There is one row per unique combination of \code{RAT_SYMBOL}, \code{RAT_GENE_RGD_ID}, \code{RAT_NCBI_GENE_ID}, \code{RAT_ENSEMBL_ID},
#'     \code{RAT_UNIPROT_ID}, \code{HUMAN_ORTHOLOG_SYMBOL}, \code{HUMAN_ORTHOLOG_RGD}, \code{HUMAN_ORTHOLOG_NCBI_GENE_ID}, and \code{HUMAN_ORTHOLOG_ENSEMBL_ID}, 
#'     so some genes correspond to multiple rows in the table. For example, there are 24 rows corresponding to \code{RAT_SYMBOL == "Tnf"} 
#'     because 3 \code{RAT_ENSEMBL_ID}s correspond to Tnf, and 8 \code{HUMAN_ORTHOLOG_ENSEMBL_ID}s correspond to TNF.  
#'     
#' @source
#'     <https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.metadata.EntrezGene.gz>
#'     <https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.metadata.HGNC.gz>
#'     <https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.chr_patch_hapl_scaff.basic.annotation.gtf.gz>
#'     <https://download.rgd.mcw.edu/pub/data_release/RGD_ORTHOLOGS.txt> (2020-01-10)
#'     <https://download.rgd.mcw.edu/data_release/RAT/GENES_RAT.txt> (2020-11-15)
"RAT_TO_HUMAN_GENE"


#' @name REPEATED_FEATURES
#' @title Repeated feature info
#' @description Documentation of 91 IMMUNO and METAB features that were measured on multiple platforms 
#'     and have multiple differential analysis results. New unique identifiers are defined. These new 
#'     unique identifiers were used in the graphical analysis. 
#' @format A data frame with 825 rows and 10 variables:
#' \describe{
#'   \item{\code{feature_ID}}{`r feature_ID()`}
#'   \item{\code{assay}}{`r assay()`}
#'   \item{\code{tissue}}{`r tissue()`}
#'   \item{\code{dataset}}{`r dataset()`}
#'   \item{\code{selection_fdr}}{`r selection_fdr()`}
#'   \item{\code{metabolite_refmet}}{character, RefMet name of metabolite}
#'   \item{\code{feature}}{character, duplicated \code{feature} in the format \code{\link{ASSAY_ABBREV};\link{TISSUE_ABBREV};[feature_ID]}}
#'   \item{\code{new_feature_ID}}{character, new unique \code{feature_ID} with \code{dataset} prepended}
#'   \item{\code{new_feature}}{character, new unique \code{feature} in the format \code{\link{ASSAY_ABBREV};\link{TISSUE_ABBREV};[new_feature_ID]}} 
#' }
#' @details 91 IMMUNO and METAB features were measured on multiple platforms and have multiple differential analysis results. 
#'     In order to distinguish between the different sets of results, these 
#'     feature_IDs were prepended with \code{dataset}. 
#'     For example, "BDNF" is separated into "rat-myokine:BDNF" and "rat-pituitary:BDNF". 
#'     For completeness, both the modified and unmodified feature_IDs are included in both 
#'     the feature-to-gene map and universe lists. Note that the graphical analysis
#'     uses the modified feature_IDs; differential analysis results use the unmodified feature_IDs.
"REPEATED_FEATURES"


#' @title Metabolite feature IDs and metadata
#' @description Mapping between various metabolite feature identifiers used
#'   at different stages of data processing. Also includes some feature metadata
#' @format A data frame with 14420 rows and 13 variables:
#' \describe{
#'   \item{\code{tissue}}{`r tissue()`}
#'   \item{\code{dataset}}{`r dataset_metab()`}
#'   \item{\code{metabolite_name}}{character, name of metabolite as appears in the Chemical Analysis Site's data}
#'   \item{\code{metabolite_refmet}}{character, RefMet name of metabolite}
#'   \item{\code{feature_ID_sample_data}}{character, \code{feature_ID} used in the sample-level data, e.g. [METAB_NORM_DATA_FLAT]}
#'   \item{\code{feature_ID_da}}{character, \code{feature_ID} used in the standard
#'     differential analysis results tables, e.g. [METAB_BAT_DA]}
#'   \item{\code{feature_ID_metareg}}{character, \code{feature_ID} used in the meta-regression
#'     results tables, e.g. [METAB_BAT_DA_METAREG]}
#'   \item{\code{dataset_metareg}}{character, version of the differential analysis results
#'     retained after meta-regression, either \code{dataset} or 'meta-reg'}
#'   \item{\code{feature}}{`r feature()`}
#'   \item{\code{rt}}{double, retention time}
#'   \item{\code{mz}}{double, mass over charge}
#'   \item{\code{neutral_mass}}{double, neutral mass}
#'   \item{\code{formula}}{character, chemical formula} 
#' }
#' @details 
#'   \code{metabolite_name} is always equal to \code{feature_ID_sample_data} and \code{feature_ID_da}.
#'   
#'   \code{feature_ID_metareg} is different than \code{metabolite_name} when the combined 
#'   results from multiple measurements of the same metabolite are reported in the meta-regression differential 
#'   analysis results, in which case \code{feature_ID_metareg} is \code{metabolite_refmet}. 
#'   112 rows (57 unique metabolites in 9 tissues) do not have corresponding meta-regression results, 
#'   in which case \code{feature_ID_metareg} is \code{NA}. 
#'   
#'   \code{feature} is only non-NA for training-regulated features at 5% IHW FDR (see [TRAINING_REGULATED_FEATURES]).
#'   Each non-NA value of \code{feature} is unique. In the case of repeated measurements, this means \code{dataset}
#'   was added to \code{metabolite_name}. See [REPEATED_FEATURES] for more details.
#'   
#'   Note that RefMet IDs are frequently updated. This table provides the version of the RefMet IDs
#'   at the time that these data were generated. 
#'   
"METAB_FEATURE_ID_MAP"


#' @title RRBS feature annotation
#' @format A data frame with 7585076 rows and 12 variables:
#' \describe{
#'   \item{\code{Chr}}{integer, chromosome}
#'   \item{\code{Locus}}{character, base pair range of feature}
#'   \item{\code{EntrezID}}{character, Entrez ID of closest gene}
#'   \item{\code{Symbol}}{character, gene symbol of closest gene}
#'   \item{\code{Strand}}{character, strand}
#'   \item{\code{Width}}{integer, width of the genomic locus represented by the feature}
#'   \item{\code{NumSites}}{integer, the number of sites merged to generate the feature, sites are merged by their correlation pattern in the data using an unsupervised analysis}
#'   \item{\code{Sites}}{character, a comma separated string with the specific sites that were merged to generate the feature}
#'   \item{\code{LocStart}}{integer, feature start in base pairs}
#'   \item{\code{LocEnd}}{integer, feature end in base pairs}
#'   \item{\code{tissue}}{`r tissue()`}
#'   \item{\code{feature_ID}}{`r feature_ID()`}
#' }
#' @details 
#'   METHYL feature annotation is only available via download from Google Cloud Storage:
#'   <https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/METHYL_FEATURE_ANNOT.rda>.
#'   You can use [MotrpacRatTraining6mo::load_methyl_feature_annotation()] to download and return this file. 
#' 
#'   Only CpG sites with methylation coverage of >=10 in all samples were included for downstream analysis,
#'   and normalization was performed separately in each tissue. Individual CpG sites were divided into 500 base-pair
#'   windows and were clustered using the Markov Clustering algorithm via the MCL R package (Jager, 2015). To apply MCL,
#'   for each 500 base-pair window an undirected graph was constructed, linking individual sites if their correlation
#'   was >=0.7. MCL was chosen for this task as it: (1) determines the number of clusters internally, (2) identifies
#'   homogeneous clusters, and (3) keeps single sites that are not correlated with either sites as singletons (clusters of size one).
#'   
#'   Given these sites, this table was generated using [MotrpacRatTraining6mo::get_peak_annotations()]. 
#'   \code{relationship_to_gene} is the shortest distance between the feature and the start or end of the closest gene. 
#'   It is 0 if the feature has any overlap with the gene. 
#'   \code{custom_annotation} fixes many issues with the \code{ChIPseeker} annotation (v1.22.1). 
#'   
#' @name METHYL_FEATURE_ANNOT
NULL


#' @title ATAC-seq feature annotation
#' @format A data frame with 1209773 rows and 16 variables:
#' \describe{
#'   \item{\code{assay}}{`r assay()`}
#'   \item{\code{assay_code}}{`r assay_code()`}
#'   \item{\code{feature_ID}}{`r feature_ID()`}
#'   \item{\code{chrom}}{character, chromosome: 1-20, X, or Y}
#'   \item{\code{start}}{double, base pair of feature start}
#'   \item{\code{end}}{double, base bair of feature end}
#'   \item{\code{width}}{integer, width of feature in base pairs}
#'   \item{\code{chipseeker_annotation}}{character, annotation from [ChIPseeker::annotatePeak()]}
#'   \item{\code{custom_annotation}}{character, a version of the ChIPseeker annotations with many corrections. Values include:
#'     "Distal Intergenic", "Promoter (<=1kb)", "Exon", "Promoter (1-2kb)", "Downstream (<5kb)", "Upstream (<5kb)", "5' UTR",
#'     "Intron", "3' UTR", "Overlaps Gene", where "Overlaps Gene" means the feature has a non-zero overlap with
#'     either the start or end of the gene but was not otherwise asssigned an annotation.}
#'   \item{\code{distanceToTSS}}{double, minimum distance from one end of the feature to the transcription start site.}
#'   \item{\code{relationship_to_gene}}{double, distance from the closest edge 
#'     of the feature to the start or end of the closest gene, whichever is closer.
#'     A value of 0 means there is non-zero overlap between the feature and the gene.
#'     A negative value means the feature is upstream of \code{geneStart}.
#'     A a positive value means the feature is downstream of \code{geneEnd}.
#'     Note that \code{geneStart} and \code{geneEnd} are strand-agnostic, i.e. \code{geneStart} 
#'     is always less than \code{geneEnd}, even if the gene is on the negative strand (\code{geneStrand == 2}).}
#'   \item{\code{ensembl_gene}}{character, Ensembl gene ID from release 96 of the Rattus norvegicus gene annotation}
#'   \item{\code{geneStart}}{integer, base pair start of gene; strand-agnostic, meaning \code{geneStart} is always less than \code{geneEnd}}
#'   \item{\code{geneEnd}}{integer, base pair end of gene; strand-agnostic, meaning \code{geneStart} is always less than \code{geneEnd}}
#'   \item{\code{geneLength}}{integer, length of gene in base pairs}
#'   \item{\code{geneStrand}}{integer, 1 (forward strand) or 2 (reverse strand)} 
#' }
#' @details 
#'   ATAC feature annotation is only available via download from Google Cloud Storage:
#'   <https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/ATAC_FEATURE_ANNOT.rda>.
#'   You can use [MotrpacRatTraining6mo::load_atac_feature_annotation()] to download and return this file. 
#' 
#'   This table was generated using [MotrpacRatTraining6mo::get_peak_annotations()]. 
#'   \code{relationship_to_gene} is the shortest distance between the feature and the start or end of the closest gene. 
#'   It is 0 if the feature has any overlap with the gene. 
#'   \code{custom_annotation} fixes many issues with the \code{ChIPseeker} annotation (v1.22.1). 
#'   
#' @name ATAC_FEATURE_ANNOT
NULL


## Phenotypic data ####

# Helper function to repeat variables for 40 days in PHENO docs
doctext = function(day){
  str = paste0(c("#'   \\item{\\code{training.day*date}}{character, date of given training day}",
                 "#'   \\item{\\code{training.day*_days}}{integer, day of given training day relative to \\code{registration.d_arrive} (count)}",
                 "#'   \\item{\\code{training.day*time}}{character, time of given training day}",
                 "#'   \\item{\\code{training.day*_treadmillspeed}}{double, treadmill speed on given training day}",
                 "#'   \\item{\\code{training.day*_treadmillincline}}{double, treadmill incline on given training day}",
                 "#'   \\item{\\code{training.day*_timeontreadmill}}{integer, time on treadmill on given training day}",
                 "#'   \\item{\\code{training.day*_weight}}{double, weight on given training day}",
                 "#'   \\item{\\code{training.day*_posttrainlact}}{double, post-training lactate on given training day}",
                 "#'   \\item{\\code{training.day*_score}}{integer, score on given training day}",
                 "#'   \\item{\\code{training.day*_comments}}{character, comments on given training day}"), collapse="\n")
  cat(gsub("\\*",day,str),"\n")
}
# for (i in 1:40){
#   doctext(i)
# }

#' @title Phenotypic data 
#' @description Phenotypic data for samples from the MoTrPAC endurance exercise 
#'     training study in 6-month-old rats. One row per sample (\code{viallabel}).
#' @format A data frame with 5955 rows and 510 variables:
#' \describe{
#'   \item{\code{pid}}{integer, unique, randomly generated 8 digit numeric identifier used in linkage to phenotypic data}
#'   \item{\code{bid}}{integer, unique, randomly generated 8 digit numeric identifier used in linkage to phenotypic data}
#'   \item{\code{labelid}}{double, unique 11 digit identifier for specimen label ID, originating at the collection site, 
#'       that provides a link to specimen processing and used for shipments to the biorepository (same as 
#'       \code{viallabel} only in instances where aliquots were not further processed at the biorepository)}
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
#'   \item{\code{training.day1date}}{character, date of given training day}
#'   \item{\code{training.day1_days}}{integer, day of given training day relative to \code{registration.d_arrive} (count)}
#'   \item{\code{training.day1time}}{character, time of given training day}
#'   \item{\code{training.day1_treadmillspeed}}{double, treadmill speed on given training day}
#'   \item{\code{training.day1_treadmillincline}}{double, treadmill incline on given training day}
#'   \item{\code{training.day1_timeontreadmill}}{integer, time on treadmill on given training day}
#'   \item{\code{training.day1_weight}}{double, weight on given training day}
#'   \item{\code{training.day1_posttrainlact}}{double, post-training lactate on given training day}
#'   \item{\code{training.day1_score}}{integer, score on given training day}
#'   \item{\code{training.day1_comments}}{character, comments on given training day} 
#'   \item{\code{training.day2date}}{character, date of given training day}
#'   \item{\code{training.day2_days}}{integer, day of given training day relative to \code{registration.d_arrive} (count)}
#'   \item{\code{training.day2time}}{character, time of given training day}
#'   \item{\code{training.day2_treadmillspeed}}{double, treadmill speed on given training day}
#'   \item{\code{training.day2_treadmillincline}}{double, treadmill incline on given training day}
#'   \item{\code{training.day2_timeontreadmill}}{integer, time on treadmill on given training day}
#'   \item{\code{training.day2_score}}{integer, score on given training day}
#'   \item{\code{training.day2_comments}}{character, comments on given training day} 
#'   \item{\code{training.day3date}}{character, date of given training day}
#'   \item{\code{training.day3_days}}{integer, day of given training day relative to \code{registration.d_arrive} (count)}
#'   \item{\code{training.day3time}}{character, time of given training day}
#'   \item{\code{training.day3_treadmillspeed}}{double, treadmill speed on given training day}
#'   \item{\code{training.day3_treadmillincline}}{double, treadmill incline on given training day}
#'   \item{\code{training.day3_timeontreadmill}}{integer, time on treadmill on given training day}
#'   \item{\code{training.day3_score}}{integer, score on given training day}
#'   \item{\code{training.day3_comments}}{character, comments on given training day} 
#'   \item{\code{training.day4date}}{character, date of given training day}
#'   \item{\code{training.day4_days}}{integer, day of given training day relative to \code{registration.d_arrive} (count)}
#'   \item{\code{training.day4time}}{character, time of given training day}
#'   \item{\code{training.day4_treadmillspeed}}{double, treadmill speed on given training day}
#'   \item{\code{training.day4_treadmillincline}}{double, treadmill incline on given training day}
#'   \item{\code{training.day4_timeontreadmill}}{integer, time on treadmill on given training day}
#'   \item{\code{training.day4_score}}{integer, score on given training day}
#'   \item{\code{training.day4_comments}}{character, comments on given training day} 
#'   \item{\code{training.day5date}}{character, date of given training day}
#'   \item{\code{training.day5_days}}{integer, day of given training day relative to \code{registration.d_arrive} (count)}
#'   \item{\code{training.day5time}}{character, time of given training day}
#'   \item{\code{training.day5_treadmillspeed}}{double, treadmill speed on given training day}
#'   \item{\code{training.day5_treadmillincline}}{double, treadmill incline on given training day}
#'   \item{\code{training.day5_timeontreadmill}}{integer, time on treadmill on given training day}
#'   \item{\code{training.day5_posttrainlact}}{double, post-training lactate on given training day}
#'   \item{\code{training.day5_score}}{integer, score on given training day}
#'   \item{\code{training.day5_comments}}{character, comments on given training day} 
#'   \item{\code{training.day6date}}{character, date of given training day}
#'   \item{\code{training.day6_days}}{integer, day of given training day relative to \code{registration.d_arrive} (count)}
#'   \item{\code{training.day6time}}{character, time of given training day}
#'   \item{\code{training.day6_treadmillspeed}}{double, treadmill speed on given training day}
#'   \item{\code{training.day6_treadmillincline}}{double, treadmill incline on given training day}
#'   \item{\code{training.day6_timeontreadmill}}{integer, time on treadmill on given training day}
#'   \item{\code{training.day6_weight}}{double, weight on given training day}
#'   \item{\code{training.day6_posttrainlact}}{double, post-training lactate on given training day}
#'   \item{\code{training.day6_score}}{integer, score on given training day}
#'   \item{\code{training.day6_comments}}{character, comments on given training day} 
#'   \item{\code{training.day7date}}{character, date of given training day}
#'   \item{\code{training.day7_days}}{integer, day of given training day relative to \code{registration.d_arrive} (count)}
#'   \item{\code{training.day7time}}{character, time of given training day}
#'   \item{\code{training.day7_treadmillspeed}}{double, treadmill speed on given training day}
#'   \item{\code{training.day7_treadmillincline}}{double, treadmill incline on given training day}
#'   \item{\code{training.day7_timeontreadmill}}{integer, time on treadmill on given training day}
#'   \item{\code{training.day7_posttrainlact}}{double, post-training lactate on given training day}
#'   \item{\code{training.day7_score}}{integer, score on given training day}
#'   \item{\code{training.day7_comments}}{character, comments on given training day} 
#'   \item{\code{training.day8date}}{character, date of given training day}
#'   \item{\code{training.day8_days}}{integer, day of given training day relative to \code{registration.d_arrive} (count)}
#'   \item{\code{training.day8time}}{character, time of given training day}
#'   \item{\code{training.day8_treadmillspeed}}{double, treadmill speed on given training day}
#'   \item{\code{training.day8_treadmillincline}}{double, treadmill incline on given training day}
#'   \item{\code{training.day8_timeontreadmill}}{integer, time on treadmill on given training day}
#'   \item{\code{training.day8_posttrainlact}}{double, post-training lactate on given training day}
#'   \item{\code{training.day8_score}}{integer, score on given training day}
#'   \item{\code{training.day8_comments}}{character, comments on given training day} 
#'   \item{\code{training.day9date}}{character, date of given training day}
#'   \item{\code{training.day9_days}}{integer, day of given training day relative to \code{registration.d_arrive} (count)}
#'   \item{\code{training.day9time}}{character, time of given training day}
#'   \item{\code{training.day9_treadmillspeed}}{double, treadmill speed on given training day}
#'   \item{\code{training.day9_treadmillincline}}{double, treadmill incline on given training day}
#'   \item{\code{training.day9_timeontreadmill}}{integer, time on treadmill on given training day}
#'   \item{\code{training.day9_posttrainlact}}{double, post-training lactate on given training day}
#'   \item{\code{training.day9_score}}{integer, score on given training day}
#'   \item{\code{training.day9_comments}}{character, comments on given training day} 
#'   \item{\code{training.day10date}}{character, date of given training day}
#'   \item{\code{training.day10_days}}{integer, day of given training day relative to \code{registration.d_arrive} (count)}
#'   \item{\code{training.day10time}}{character, time of given training day}
#'   \item{\code{training.day10_treadmillspeed}}{double, treadmill speed on given training day}
#'   \item{\code{training.day10_treadmillincline}}{double, treadmill incline on given training day}
#'   \item{\code{training.day10_timeontreadmill}}{integer, time on treadmill on given training day}
#'   \item{\code{training.day10_posttrainlact}}{double, post-training lactate on given training day}
#'   \item{\code{training.day10_score}}{integer, score on given training day}
#'   \item{\code{training.day10_comments}}{character, comments on given training day} 
#'   \item{\code{training.day11date}}{character, date of given training day}
#'   \item{\code{training.day11_days}}{integer, day of given training day relative to \code{registration.d_arrive} (count)}
#'   \item{\code{training.day11time}}{character, time of given training day}
#'   \item{\code{training.day11_treadmillspeed}}{double, treadmill speed on given training day}
#'   \item{\code{training.day11_treadmillincline}}{double, treadmill incline on given training day}
#'   \item{\code{training.day11_timeontreadmill}}{integer, time on treadmill on given training day}
#'   \item{\code{training.day11_weight}}{double, weight on given training day}
#'   \item{\code{training.day11_posttrainlact}}{double, post-training lactate on given training day}
#'   \item{\code{training.day11_score}}{integer, score on given training day}
#'   \item{\code{training.day11_comments}}{character, comments on given training day} 
#'   \item{\code{training.day12date}}{character, date of given training day}
#'   \item{\code{training.day12_days}}{integer, day of given training day relative to \code{registration.d_arrive} (count)}
#'   \item{\code{training.day12time}}{character, time of given training day}
#'   \item{\code{training.day12_treadmillspeed}}{double, treadmill speed on given training day}
#'   \item{\code{training.day12_treadmillincline}}{double, treadmill incline on given training day}
#'   \item{\code{training.day12_timeontreadmill}}{integer, time on treadmill on given training day}
#'   \item{\code{training.day12_score}}{integer, score on given training day}
#'   \item{\code{training.day12_comments}}{character, comments on given training day} 
#'   \item{\code{training.day13date}}{character, date of given training day}
#'   \item{\code{training.day13_days}}{integer, day of given training day relative to \code{registration.d_arrive} (count)}
#'   \item{\code{training.day13time}}{character, time of given training day}
#'   \item{\code{training.day13_treadmillspeed}}{double, treadmill speed on given training day}
#'   \item{\code{training.day13_treadmillincline}}{double, treadmill incline on given training day}
#'   \item{\code{training.day13_timeontreadmill}}{integer, time on treadmill on given training day}
#'   \item{\code{training.day13_score}}{integer, score on given training day}
#'   \item{\code{training.day13_comments}}{character, comments on given training day} 
#'   \item{\code{training.day14date}}{character, date of given training day}
#'   \item{\code{training.day14_days}}{integer, day of given training day relative to \code{registration.d_arrive} (count)}
#'   \item{\code{training.day14time}}{character, time of given training day}
#'   \item{\code{training.day14_treadmillspeed}}{double, treadmill speed on given training day}
#'   \item{\code{training.day14_treadmillincline}}{double, treadmill incline on given training day}
#'   \item{\code{training.day14_timeontreadmill}}{integer, time on treadmill on given training day}
#'   \item{\code{training.day14_score}}{integer, score on given training day}
#'   \item{\code{training.day14_comments}}{character, comments on given training day} 
#'   \item{\code{training.day15date}}{character, date of given training day}
#'   \item{\code{training.day15_days}}{integer, day of given training day relative to \code{registration.d_arrive} (count)}
#'   \item{\code{training.day15time}}{character, time of given training day}
#'   \item{\code{training.day15_treadmillspeed}}{double, treadmill speed on given training day}
#'   \item{\code{training.day15_treadmillincline}}{double, treadmill incline on given training day}
#'   \item{\code{training.day15_timeontreadmill}}{integer, time on treadmill on given training day}
#'   \item{\code{training.day15_posttrainlact}}{double, post-training lactate on given training day}
#'   \item{\code{training.day15_score}}{integer, score on given training day}
#'   \item{\code{training.day15_comments}}{character, comments on given training day} 
#'   \item{\code{training.day16date}}{character, date of given training day}
#'   \item{\code{training.day16_days}}{integer, day of given training day relative to \code{registration.d_arrive} (count)}
#'   \item{\code{training.day16time}}{character, time of given training day}
#'   \item{\code{training.day16_treadmillspeed}}{double, treadmill speed on given training day}
#'   \item{\code{training.day16_treadmillincline}}{double, treadmill incline on given training day}
#'   \item{\code{training.day16_timeontreadmill}}{integer, time on treadmill on given training day}
#'   \item{\code{training.day16_weight}}{double, weight on given training day}
#'   \item{\code{training.day16_posttrainlact}}{double, post-training lactate on given training day}
#'   \item{\code{training.day16_score}}{integer, score on given training day}
#'   \item{\code{training.day16_comments}}{character, comments on given training day} 
#'   \item{\code{training.day17date}}{character, date of given training day}
#'   \item{\code{training.day17_days}}{integer, day of given training day relative to \code{registration.d_arrive} (count)}
#'   \item{\code{training.day17time}}{character, time of given training day}
#'   \item{\code{training.day17_treadmillspeed}}{double, treadmill speed on given training day}
#'   \item{\code{training.day17_treadmillincline}}{double, treadmill incline on given training day}
#'   \item{\code{training.day17_timeontreadmill}}{integer, time on treadmill on given training day}
#'   \item{\code{training.day17_score}}{integer, score on given training day}
#'   \item{\code{training.day17_comments}}{character, comments on given training day} 
#'   \item{\code{training.day18date}}{character, date of given training day}
#'   \item{\code{training.day18_days}}{integer, day of given training day relative to \code{registration.d_arrive} (count)}
#'   \item{\code{training.day18time}}{character, time of given training day}
#'   \item{\code{training.day18_treadmillspeed}}{double, treadmill speed on given training day}
#'   \item{\code{training.day18_treadmillincline}}{double, treadmill incline on given training day}
#'   \item{\code{training.day18_timeontreadmill}}{integer, time on treadmill on given training day}
#'   \item{\code{training.day18_score}}{integer, score on given training day}
#'   \item{\code{training.day18_comments}}{character, comments on given training day} 
#'   \item{\code{training.day19date}}{character, date of given training day}
#'   \item{\code{training.day19_days}}{integer, day of given training day relative to \code{registration.d_arrive} (count)}
#'   \item{\code{training.day19time}}{character, time of given training day}
#'   \item{\code{training.day19_treadmillspeed}}{double, treadmill speed on given training day}
#'   \item{\code{training.day19_treadmillincline}}{double, treadmill incline on given training day}
#'   \item{\code{training.day19_timeontreadmill}}{integer, time on treadmill on given training day}
#'   \item{\code{training.day19_score}}{integer, score on given training day}
#'   \item{\code{training.day19_comments}}{character, comments on given training day} 
#'   \item{\code{training.day20date}}{character, date of given training day}
#'   \item{\code{training.day20_days}}{integer, day of given training day relative to \code{registration.d_arrive} (count)}
#'   \item{\code{training.day20time}}{character, time of given training day}
#'   \item{\code{training.day20_treadmillspeed}}{double, treadmill speed on given training day}
#'   \item{\code{training.day20_treadmillincline}}{double, treadmill incline on given training day}
#'   \item{\code{training.day20_timeontreadmill}}{integer, time on treadmill on given training day}
#'   \item{\code{training.day20_posttrainlact}}{double, post-training lactate on given training day}
#'   \item{\code{training.day20_score}}{integer, score on given training day}
#'   \item{\code{training.day20_comments}}{character, comments on given training day} 
#'   \item{\code{training.day21date}}{character, date of given training day}
#'   \item{\code{training.day21_days}}{integer, day of given training day relative to \code{registration.d_arrive} (count)}
#'   \item{\code{training.day21time}}{character, time of given training day}
#'   \item{\code{training.day21_treadmillspeed}}{double, treadmill speed on given training day}
#'   \item{\code{training.day21_treadmillincline}}{double, treadmill incline on given training day}
#'   \item{\code{training.day21_timeontreadmill}}{integer, time on treadmill on given training day}
#'   \item{\code{training.day21_weight}}{double, weight on given training day}
#'   \item{\code{training.day21_posttrainlact}}{double, post-training lactate on given training day}
#'   \item{\code{training.day21_score}}{integer, score on given training day}
#'   \item{\code{training.day21_comments}}{character, comments on given training day} 
#'   \item{\code{training.day22date}}{character, date of given training day}
#'   \item{\code{training.day22_days}}{integer, day of given training day relative to \code{registration.d_arrive} (count)}
#'   \item{\code{training.day22time}}{character, time of given training day}
#'   \item{\code{training.day22_treadmillspeed}}{double, treadmill speed on given training day}
#'   \item{\code{training.day22_treadmillincline}}{double, treadmill incline on given training day}
#'   \item{\code{training.day22_timeontreadmill}}{integer, time on treadmill on given training day}
#'   \item{\code{training.day22_score}}{integer, score on given training day}
#'   \item{\code{training.day22_comments}}{character, comments on given training day} 
#'   \item{\code{training.day23date}}{character, date of given training day}
#'   \item{\code{training.day23_days}}{integer, day of given training day relative to \code{registration.d_arrive} (count)}
#'   \item{\code{training.day23time}}{character, time of given training day}
#'   \item{\code{training.day23_treadmillspeed}}{double, treadmill speed on given training day}
#'   \item{\code{training.day23_treadmillincline}}{double, treadmill incline on given training day}
#'   \item{\code{training.day23_timeontreadmill}}{integer, time on treadmill on given training day}
#'   \item{\code{training.day23_score}}{integer, score on given training day}
#'   \item{\code{training.day23_comments}}{character, comments on given training day} 
#'   \item{\code{training.day24date}}{character, date of given training day}
#'   \item{\code{training.day24_days}}{integer, day of given training day relative to \code{registration.d_arrive} (count)}
#'   \item{\code{training.day24time}}{character, time of given training day}
#'   \item{\code{training.day24_treadmillspeed}}{double, treadmill speed on given training day}
#'   \item{\code{training.day24_treadmillincline}}{double, treadmill incline on given training day}
#'   \item{\code{training.day24_timeontreadmill}}{integer, time on treadmill on given training day}
#'   \item{\code{training.day24_score}}{integer, score on given training day}
#'   \item{\code{training.day24_comments}}{character, comments on given training day} 
#'   \item{\code{training.day25date}}{character, date of given training day}
#'   \item{\code{training.day25_days}}{integer, day of given training day relative to \code{registration.d_arrive} (count)}
#'   \item{\code{training.day25time}}{character, time of given training day}
#'   \item{\code{training.day25_treadmillspeed}}{double, treadmill speed on given training day}
#'   \item{\code{training.day25_treadmillincline}}{double, treadmill incline on given training day}
#'   \item{\code{training.day25_timeontreadmill}}{integer, time on treadmill on given training day}
#'   \item{\code{training.day25_posttrainlact}}{double, post-training lactate on given training day}
#'   \item{\code{training.day25_score}}{integer, score on given training day}
#'   \item{\code{training.day25_comments}}{character, comments on given training day} 
#'   \item{\code{training.day26date}}{character, date of given training day}
#'   \item{\code{training.day26_days}}{integer, day of given training day relative to \code{registration.d_arrive} (count)}
#'   \item{\code{training.day26time}}{character, time of given training day}
#'   \item{\code{training.day26_treadmillspeed}}{double, treadmill speed on given training day}
#'   \item{\code{training.day26_treadmillincline}}{double, treadmill incline on given training day}
#'   \item{\code{training.day26_timeontreadmill}}{integer, time on treadmill on given training day}
#'   \item{\code{training.day26_weight}}{double, weight on given training day}
#'   \item{\code{training.day26_posttrainlact}}{double, post-training lactate on given training day}
#'   \item{\code{training.day26_score}}{integer, score on given training day}
#'   \item{\code{training.day26_comments}}{character, comments on given training day} 
#'   \item{\code{training.day27date}}{character, date of given training day}
#'   \item{\code{training.day27_days}}{integer, day of given training day relative to \code{registration.d_arrive} (count)}
#'   \item{\code{training.day27time}}{character, time of given training day}
#'   \item{\code{training.day27_treadmillspeed}}{double, treadmill speed on given training day}
#'   \item{\code{training.day27_treadmillincline}}{double, treadmill incline on given training day}
#'   \item{\code{training.day27_timeontreadmill}}{integer, time on treadmill on given training day}
#'   \item{\code{training.day27_score}}{integer, score on given training day}
#'   \item{\code{training.day27_comments}}{character, comments on given training day} 
#'   \item{\code{training.day28date}}{character, date of given training day}
#'   \item{\code{training.day28_days}}{integer, day of given training day relative to \code{registration.d_arrive} (count)}
#'   \item{\code{training.day28time}}{character, time of given training day}
#'   \item{\code{training.day28_treadmillspeed}}{double, treadmill speed on given training day}
#'   \item{\code{training.day28_treadmillincline}}{double, treadmill incline on given training day}
#'   \item{\code{training.day28_timeontreadmill}}{integer, time on treadmill on given training day}
#'   \item{\code{training.day28_score}}{integer, score on given training day}
#'   \item{\code{training.day28_comments}}{character, comments on given training day} 
#'   \item{\code{training.day29date}}{character, date of given training day}
#'   \item{\code{training.day29_days}}{integer, day of given training day relative to \code{registration.d_arrive} (count)}
#'   \item{\code{training.day29time}}{character, time of given training day}
#'   \item{\code{training.day29_treadmillspeed}}{double, treadmill speed on given training day}
#'   \item{\code{training.day29_treadmillincline}}{double, treadmill incline on given training day}
#'   \item{\code{training.day29_timeontreadmill}}{integer, time on treadmill on given training day}
#'   \item{\code{training.day29_score}}{integer, score on given training day}
#'   \item{\code{training.day29_comments}}{character, comments on given training day} 
#'   \item{\code{training.day30date}}{character, date of given training day}
#'   \item{\code{training.day30_days}}{integer, day of given training day relative to \code{registration.d_arrive} (count)}
#'   \item{\code{training.day30time}}{character, time of given training day}
#'   \item{\code{training.day30_treadmillspeed}}{double, treadmill speed on given training day}
#'   \item{\code{training.day30_treadmillincline}}{double, treadmill incline on given training day}
#'   \item{\code{training.day30_timeontreadmill}}{integer, time on treadmill on given training day}
#'   \item{\code{training.day30_posttrainlact}}{double, post-training lactate on given training day}
#'   \item{\code{training.day30_score}}{integer, score on given training day}
#'   \item{\code{training.day30_comments}}{character, comments on given training day} 
#'   \item{\code{training.day31date}}{character, date of given training day}
#'   \item{\code{training.day31_days}}{integer, day of given training day relative to \code{registration.d_arrive} (count)}
#'   \item{\code{training.day31time}}{character, time of given training day}
#'   \item{\code{training.day31_treadmillspeed}}{double, treadmill speed on given training day}
#'   \item{\code{training.day31_treadmillincline}}{double, treadmill incline on given training day}
#'   \item{\code{training.day31_timeontreadmill}}{integer, time on treadmill on given training day}
#'   \item{\code{training.day31_weight}}{double, weight on given training day}
#'   \item{\code{training.day31_posttrainlact}}{double, post-training lactate on given training day}
#'   \item{\code{training.day31_score}}{integer, score on given training day}
#'   \item{\code{training.day31_comments}}{character, comments on given training day} 
#'   \item{\code{training.day32date}}{character, date of given training day}
#'   \item{\code{training.day32_days}}{integer, day of given training day relative to \code{registration.d_arrive} (count)}
#'   \item{\code{training.day32time}}{character, time of given training day}
#'   \item{\code{training.day32_treadmillspeed}}{double, treadmill speed on given training day}
#'   \item{\code{training.day32_treadmillincline}}{double, treadmill incline on given training day}
#'   \item{\code{training.day32_timeontreadmill}}{integer, time on treadmill on given training day}
#'   \item{\code{training.day32_score}}{integer, score on given training day}
#'   \item{\code{training.day32_comments}}{character, comments on given training day} 
#'   \item{\code{training.day33date}}{character, date of given training day}
#'   \item{\code{training.day33_days}}{integer, day of given training day relative to \code{registration.d_arrive} (count)}
#'   \item{\code{training.day33time}}{character, time of given training day}
#'   \item{\code{training.day33_treadmillspeed}}{double, treadmill speed on given training day}
#'   \item{\code{training.day33_treadmillincline}}{double, treadmill incline on given training day}
#'   \item{\code{training.day33_timeontreadmill}}{integer, time on treadmill on given training day}
#'   \item{\code{training.day33_score}}{integer, score on given training day}
#'   \item{\code{training.day33_comments}}{character, comments on given training day} 
#'   \item{\code{training.day34date}}{character, date of given training day}
#'   \item{\code{training.day34_days}}{integer, day of given training day relative to \code{registration.d_arrive} (count)}
#'   \item{\code{training.day34time}}{character, time of given training day}
#'   \item{\code{training.day34_treadmillspeed}}{double, treadmill speed on given training day}
#'   \item{\code{training.day34_treadmillincline}}{double, treadmill incline on given training day}
#'   \item{\code{training.day34_timeontreadmill}}{integer, time on treadmill on given training day}
#'   \item{\code{training.day34_score}}{integer, score on given training day}
#'   \item{\code{training.day34_comments}}{character, comments on given training day} 
#'   \item{\code{training.day35date}}{character, date of given training day}
#'   \item{\code{training.day35_days}}{integer, day of given training day relative to \code{registration.d_arrive} (count)}
#'   \item{\code{training.day35time}}{character, time of given training day}
#'   \item{\code{training.day35_treadmillspeed}}{double, treadmill speed on given training day}
#'   \item{\code{training.day35_treadmillincline}}{double, treadmill incline on given training day}
#'   \item{\code{training.day35_timeontreadmill}}{integer, time on treadmill on given training day}
#'   \item{\code{training.day35_posttrainlact}}{double, post-training lactate on given training day}
#'   \item{\code{training.day35_score}}{integer, score on given training day}
#'   \item{\code{training.day35_comments}}{character, comments on given training day} 
#'   \item{\code{training.day36date}}{character, date of given training day}
#'   \item{\code{training.day36_days}}{integer, day of given training day relative to \code{registration.d_arrive} (count)}
#'   \item{\code{training.day36time}}{character, time of given training day}
#'   \item{\code{training.day36_weight}}{double, weight on given training day}
#'   \item{\code{training.day36_comments}}{character, comments on given training day} 
#'   \item{\code{training.day37date}}{character, date of given training day}
#'   \item{\code{training.day37_days}}{integer, day of given training day relative to \code{registration.d_arrive} (count)}
#'   \item{\code{training.day37time}}{character, time of given training day}
#'   \item{\code{training.day37_treadmillspeed}}{double, treadmill speed on given training day}
#'   \item{\code{training.day37_treadmillincline}}{double, treadmill incline on given training day}
#'   \item{\code{training.day37_timeontreadmill}}{integer, time on treadmill on given training day}
#'   \item{\code{training.day37_score}}{integer, score on given training day}
#'   \item{\code{training.day37_comments}}{character, comments on given training day} 
#'   \item{\code{training.day38date}}{character, date of given training day}
#'   \item{\code{training.day38_days}}{integer, day of given training day relative to \code{registration.d_arrive} (count)}
#'   \item{\code{training.day38time}}{character, time of given training day}
#'   \item{\code{training.day38_treadmillspeed}}{double, treadmill speed on given training day}
#'   \item{\code{training.day38_treadmillincline}}{double, treadmill incline on given training day}
#'   \item{\code{training.day38_timeontreadmill}}{integer, time on treadmill on given training day}
#'   \item{\code{training.day38_score}}{integer, score on given training day}
#'   \item{\code{training.day38_comments}}{character, comments on given training day} 
#'   \item{\code{training.day39date}}{character, date of given training day}
#'   \item{\code{training.day39_days}}{integer, day of given training day relative to \code{registration.d_arrive} (count)}
#'   \item{\code{training.day39time}}{character, time of given training day}
#'   \item{\code{training.day39_treadmillspeed}}{double, treadmill speed on given training day}
#'   \item{\code{training.day39_treadmillincline}}{double, treadmill incline on given training day}
#'   \item{\code{training.day39_timeontreadmill}}{integer, time on treadmill on given training day}
#'   \item{\code{training.day39_posttrainlact}}{double, post-training lactate on given training day}
#'   \item{\code{training.day39_score}}{integer, score on given training day}
#'   \item{\code{training.day39_comments}}{character, comments on given training day} 
#'   \item{\code{training.day40date}}{character, date of given training day}
#'   \item{\code{training.day40_days}}{integer, day of given training day relative to \code{registration.d_arrive} (count)}
#'   \item{\code{training.day40time}}{character, time of given training day}
#'   \item{\code{training.day40_treadmillspeed}}{double, treadmill speed on given training day}
#'   \item{\code{training.day40_treadmillincline}}{double, treadmill incline on given training day}
#'   \item{\code{training.day40_timeontreadmill}}{integer, time on treadmill on given training day}
#'   \item{\code{training.day40_posttrainlact}}{double, post-training lactate on given training day}
#'   \item{\code{training.day40_score}}{integer, score on given training day}
#'   \item{\code{training.day40_comments}}{character, comments on given training day} 
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
#'   \item{\code{sex}}{`r sex()`}
#'   \item{\code{time_to_freeze}}{integer, sample time to freeze, \code{calculated.variables.frozetime_after_train - calculated.variables.deathtime_after_train}} 
#'   \item{\code{tissue_code_no}}{character, BIC tissue code. See [MotrpacBicQC::bic_animal_tissue_code].}
#'   \item{\code{tissue}}{`r tissue()`} 
#' }
"PHENO"


## RNA-seq sample-level data ####

#' @title RNA-seq raw counts
#' @description RNA-seq raw counts as quantified by RSEM
#' @format A data frame with 32883 rows (\code{feature_ID}) and samples in columns (\code{viallabel})
#' @details STAR (v2.7.0d) was used to index and align reads to release 96 
#'     of the Ensembl Rattus norvegicus (rn6) genome and gene annotation. 
#'     RSEM (v1.3.1) was used to quantify transcriptome-coordinate-sorted 
#'     alignments using a forward probability of 0.5 to indicate a non-strand-specific protocol. 
#' @name TRNSCRPT_RAW_COUNTS
"TRNSCRPT_ADRNL_RAW_COUNTS"

#' @rdname TRNSCRPT_RAW_COUNTS
"TRNSCRPT_BAT_RAW_COUNTS"

#' @rdname TRNSCRPT_RAW_COUNTS
"TRNSCRPT_BLOOD_RAW_COUNTS"

#' @rdname TRNSCRPT_RAW_COUNTS
"TRNSCRPT_COLON_RAW_COUNTS"

#' @rdname TRNSCRPT_RAW_COUNTS
"TRNSCRPT_CORTEX_RAW_COUNTS"

#' @rdname TRNSCRPT_RAW_COUNTS
"TRNSCRPT_HEART_RAW_COUNTS"

#' @rdname TRNSCRPT_RAW_COUNTS
"TRNSCRPT_HIPPOC_RAW_COUNTS"

#' @rdname TRNSCRPT_RAW_COUNTS
"TRNSCRPT_HYPOTH_RAW_COUNTS"

#' @rdname TRNSCRPT_RAW_COUNTS
"TRNSCRPT_KIDNEY_RAW_COUNTS"

#' @rdname TRNSCRPT_RAW_COUNTS
"TRNSCRPT_LIVER_RAW_COUNTS"

#' @rdname TRNSCRPT_RAW_COUNTS
"TRNSCRPT_LUNG_RAW_COUNTS"

#' @rdname TRNSCRPT_RAW_COUNTS
"TRNSCRPT_OVARY_RAW_COUNTS"

#' @rdname TRNSCRPT_RAW_COUNTS
"TRNSCRPT_SKMGN_RAW_COUNTS"

#' @rdname TRNSCRPT_RAW_COUNTS
"TRNSCRPT_SKMVL_RAW_COUNTS"

#' @rdname TRNSCRPT_RAW_COUNTS
"TRNSCRPT_SMLINT_RAW_COUNTS"

#' @rdname TRNSCRPT_RAW_COUNTS
"TRNSCRPT_SPLEEN_RAW_COUNTS"

#' @rdname TRNSCRPT_RAW_COUNTS
"TRNSCRPT_TESTES_RAW_COUNTS"

#' @rdname TRNSCRPT_RAW_COUNTS
"TRNSCRPT_VENACV_RAW_COUNTS"

#' @rdname TRNSCRPT_RAW_COUNTS
"TRNSCRPT_WATSC_RAW_COUNTS"


#' @title Normalized RNA-seq data
#' @description Normalized sample-level RNA-seq (TRNSCRPT) data used for visualization
#' @format A data frame with genes in rows (\code{feature_ID}) and samples in columns (\code{viallabel}) 
#' @details Filtering of lowly expressed genes and normalization were performed 
#'   separately in each tissue. RSEM gene counts were used to remove lowly expressed genes, 
#'   defined as having 0.5 or fewer counts per million in all but one sample. 
#'   To generate normalized sample-level data, filtered gene counts were 
#'   TMM-normalized using [edgeR::calcNormFactors()], followed by conversion to log counts per million with [edgeR::cpm()].
#' @source \code{gs://motrpac-data-freeze-pass/pass1b-06/v1.1/analysis/transcriptomics/transcript-rna-seq/normalized-data/*normalized-log-cpm*}
#' @name TRNSCRPT_NORM_DATA
"TRNSCRPT_BLOOD_NORM_DATA"

#' @rdname TRNSCRPT_NORM_DATA
"TRNSCRPT_HIPPOC_NORM_DATA"

#' @rdname TRNSCRPT_NORM_DATA
"TRNSCRPT_CORTEX_NORM_DATA"

#' @rdname TRNSCRPT_NORM_DATA
"TRNSCRPT_HYPOTH_NORM_DATA"

#' @rdname TRNSCRPT_NORM_DATA
"TRNSCRPT_SKMGN_NORM_DATA"

#' @rdname TRNSCRPT_NORM_DATA
"TRNSCRPT_SKMVL_NORM_DATA"

#' @rdname TRNSCRPT_NORM_DATA
"TRNSCRPT_HEART_NORM_DATA"

#' @rdname TRNSCRPT_NORM_DATA
"TRNSCRPT_KIDNEY_NORM_DATA"

#' @rdname TRNSCRPT_NORM_DATA
"TRNSCRPT_ADRNL_NORM_DATA"

#' @rdname TRNSCRPT_NORM_DATA
"TRNSCRPT_COLON_NORM_DATA"

#' @rdname TRNSCRPT_NORM_DATA
"TRNSCRPT_SPLEEN_NORM_DATA"

#' @rdname TRNSCRPT_NORM_DATA
"TRNSCRPT_TESTES_NORM_DATA"

#' @rdname TRNSCRPT_NORM_DATA
"TRNSCRPT_OVARY_NORM_DATA"

#' @rdname TRNSCRPT_NORM_DATA
"TRNSCRPT_VENACV_NORM_DATA"

#' @rdname TRNSCRPT_NORM_DATA
"TRNSCRPT_LUNG_NORM_DATA"

#' @rdname TRNSCRPT_NORM_DATA
"TRNSCRPT_SMLINT_NORM_DATA"

#' @rdname TRNSCRPT_NORM_DATA
"TRNSCRPT_LIVER_NORM_DATA"

#' @rdname TRNSCRPT_NORM_DATA
"TRNSCRPT_BAT_NORM_DATA"

#' @rdname TRNSCRPT_NORM_DATA
"TRNSCRPT_WATSC_NORM_DATA"


## ATAC-seq sample-level data ####

#' @title ATAC-seq raw counts
#' @format A data frame with chromatin accessibility peaks in rows and samples (vial labels) in columns
#' @details 
#'   Unfiltered ATAC sample-level data are only available via download from Google Cloud Storage. 
#'   For example, <https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/ATAC_BAT_RAW_COUNTS.rda> 
#'   is the file for brown adipose tissue (BAT) data. You can change the name of the file to specify other tissues including:
#'   HEART, HIPPOC, KIDNEY, LIVER, LUNG, SKMGN (gastrocnemius skeletal muscle), and WATSC (subcutaneous white adipose tissue).
#'   You can also use [MotrpacRatTraining6mo::load_sample_data()] or [MotrpacRatTraining6mo::get_rdata_from_url()] 
#'   to download raw and normalized sample-level data for ATAC and METHYL.
#'   For more details about these files see the readme of this repository at 
#'   <https://github.com/MoTrPAC/MotrpacRatTraining6mo/blob/main/README.md>. 
#' 
#'   Data was processed with the [ENCODE ATAC-seq pipeline (v1.7.0)](https://github.com/ENCODE-DCC/atac-seq-pipeline).
#'   Samples from a single sex and training time point, e.g., males trained for 2 weeks, were analyzed together as biological
#'   replicates in a single workflow. Briefly, adapters were trimmed with cutadapt v2.5 (Martin, 2011) and aligned to release 96
#'   of the Ensembl Rattus norvegicus (rn6) genome (Dobin et al., 2013) with Bowtie 2 v2.3.4.3 (Langmead and Salzberg, 2012).
#'   Duplicate reads and reads mapping to the mitochondrial chromosome were removed. Signal files and peak calls were generated
#'   using MACS2 v2.2.4 (Gaspar, 2018), both from reads from each sample and pooled reads from all biological replicates.
#'   Pooled peaks were compared with the peaks called for each replicate individually using Irreproducibility Discovery Rate (Li et al., 2011)
#'   and thresholded to generate an optimal set of peaks.
#'
#'   The cloud implementation of the ENCODE ATAC-seq pipeline and source code for the post-processing steps are available at 
#'   <https://github.com/MoTrPAC/motrpac-atac-seq-pipeline>.
#'   Optimal peaks (overlap.optimal_peak.narrowPeak.bed.gz) from all workflows were concatenated, trimmed to 200 base pairs around the summit,
#'   and sorted and merged with bedtools v2.29.0 (Quinlan and Hall, 2010) to generate a master peak list. This peak list was intersected with
#'   the filtered alignments from each sample using bedtools coverage with options \code{-nonamecheck} and \code{-counts} to generate a peak by
#'   sample matrix of raw counts.
#'   
#' @name ATAC_RAW_COUNTS
NULL


#' @title Normalized ATAC-seq data
#' @description Normalized sample-level ATAC-seq (ATAC) data used for visualization and differential analysis
#' @format A data frame with peaks in rows (\code{feature_ID}) and samples in columns (\code{viallabel})
#' @details 
#'   Unfiltered ATAC sample-level data are only available via download from Google Cloud Storage. 
#'   For example, <https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/ATAC_BAT_NORM_DATA.rda> 
#'   is the file for brown adipose tissue (BAT) data. You can change the name of the file to specify other tissues including:
#'   HEART, HIPPOC, KIDNEY, LIVER, LUNG, SKMGN (gastrocnemius skeletal muscle), and WATSC (subcutaneous white adipose tissue).
#'   You can also use [MotrpacRatTraining6mo::load_sample_data()] or [MotrpacRatTraining6mo::get_rdata_from_url()] 
#'   to download raw and normalized sample-level data for ATAC and METHYL.
#'   For more details about these files see the readme of this repository at 
#'   <https://github.com/MoTrPAC/MotrpacRatTraining6mo/blob/main/README.md>. 
#' 
#'   Data was processed with the [ENCODE ATAC-seq pipeline (v1.7.0)](https://github.com/ENCODE-DCC/atac-seq-pipeline).
#'   Samples from a single sex and training time point, e.g., males trained for 2 weeks, were analyzed together as biological
#'   replicates in a single workflow. Briefly, adapters were trimmed with cutadapt v2.5 (Martin, 2011) and aligned to release 96
#'   of the Ensembl Rattus norvegicus (rn6) genome (Dobin et al., 2013) with Bowtie 2 v2.3.4.3 (Langmead and Salzberg, 2012).
#'   Duplicate reads and reads mapping to the mitochondrial chromosome were removed. Signal files and peak calls were generated
#'   using MACS2 v2.2.4 (Gaspar, 2018), both from reads from each sample and pooled reads from all biological replicates.
#'   Pooled peaks were compared with the peaks called for each replicate individually using Irreproducibility Discovery Rate (Li et al., 2011)
#'   and thresholded to generate an optimal set of peaks.
#'
#'   The cloud implementation of the ENCODE ATAC-seq pipeline and source code for the post-processing steps are available at 
#'   <https://github.com/MoTrPAC/motrpac-atac-seq-pipeline>.
#'   Optimal peaks (overlap.optimal_peak.narrowPeak.bed.gz) from all workflows were concatenated, trimmed to 200 base pairs around the summit,
#'   and sorted and merged with bedtools v2.29.0 (Quinlan and Hall, 2010) to generate a master peak list. This peak list was intersected with
#'   the filtered alignments from each sample using bedtools coverage with options \code{-nonamecheck} and \code{-counts} to generate a peak by
#'   sample matrix of raw counts.
#'
#'   The remaining steps were applied separately on raw counts from each tissue. Peaks from non-autosomal chromosomes were removed,
#'   as well as peaks that did not have at least 10 read counts in four samples. Filtered raw counts were then quantile-normalized with
#'   limma-voom (Law et al., 2014). 
#'   
#'   For the subset of normalized data corresponding to training-regulated features at 5% IHW FDR, see [ATAC_NORM_DATA_05FDR].
#'   
#' @source \code{gs://motrpac-data-freeze-pass/pass1b-06/v1.1/analysis/epigenomics/epigen-atac-seq/normalized-data/*quant-norm*}
#' @name ATAC_NORM_DATA
NULL


#' @title Normalized ATAC-seq data for training-regulated features
#' @description Normalized sample-level ATAC-seq (ATAC) data used for visualization and differential analysis.
#'   Only for training-regulated features at 5% IHW FDR. For sample-level data for \emph{all} features, 
#'   see [ATAC_NORM_DATA] and [ATAC_RAW_COUNTS]. 
#' @format A data frame with peaks in rows (\code{feature_ID}) and samples in columns (\code{viallabel})
#' @details Data was processed with the [ENCODE ATAC-seq pipeline (v1.7.0)](https://github.com/ENCODE-DCC/atac-seq-pipeline).
#'   Samples from a single sex and training time point, e.g., males trained for 2 weeks, were analyzed together as biological
#'   replicates in a single workflow. Briefly, adapters were trimmed with cutadapt v2.5 (Martin, 2011) and aligned to release 96
#'   of the Ensembl Rattus norvegicus (rn6) genome (Dobin et al., 2013) with Bowtie 2 v2.3.4.3 (Langmead and Salzberg, 2012).
#'   Duplicate reads and reads mapping to the mitochondrial chromosome were removed. Signal files and peak calls were generated
#'   using MACS2 v2.2.4 (Gaspar, 2018), both from reads from each sample and pooled reads from all biological replicates.
#'   Pooled peaks were compared with the peaks called for each replicate individually using Irreproducibility Discovery Rate (Li et al., 2011)
#'   and thresholded to generate an optimal set of peaks.
#'
#'   The cloud implementation of the ENCODE ATAC-seq pipeline and source code for the post-processing steps are available at 
#'   <https://github.com/MoTrPAC/motrpac-atac-seq-pipeline>.
#'   Optimal peaks (overlap.optimal_peak.narrowPeak.bed.gz) from all workflows were concatenated, trimmed to 200 base pairs around the summit,
#'   and sorted and merged with bedtools v2.29.0 (Quinlan and Hall, 2010) to generate a master peak list. This peak list was intersected with
#'   the filtered alignments from each sample using bedtools coverage with options \code{-nonamecheck} and \code{-counts} to generate a peak by
#'   sample matrix of raw counts.
#'
#'   The remaining steps were applied separately on raw counts from each tissue. Peaks from non-autosomal chromosomes were removed,
#'   as well as peaks that did not have at least 10 read counts in four samples. Filtered raw counts were then quantile-normalized with
#'   limma-voom (Law et al., 2014). 
#'   
#'   After performing differential analysis, training-regulated features were selected at 5% IHW FDR. The normalized data were
#'   filtered down to these features and provided here. 
#'   
#'   For the full set of normalized sample-level data, see [ATAC_NORM_DATA] and [ATAC_RAW_COUNTS]. 
#'
#' @source \code{gs://motrpac-data-freeze-pass/pass1b-06/v1.1/analysis/epigenomics/epigen-atac-seq/normalized-data/*quant-norm*}
#' @name ATAC_NORM_DATA_05FDR
"ATAC_HIPPOC_NORM_DATA_05FDR"

#' @rdname ATAC_NORM_DATA_05FDR
"ATAC_SKMGN_NORM_DATA_05FDR"

#' @rdname ATAC_NORM_DATA_05FDR
"ATAC_HEART_NORM_DATA_05FDR"

#' @rdname ATAC_NORM_DATA_05FDR
"ATAC_KIDNEY_NORM_DATA_05FDR"

#' @rdname ATAC_NORM_DATA_05FDR
"ATAC_LUNG_NORM_DATA_05FDR"

#' @rdname ATAC_NORM_DATA_05FDR
"ATAC_LIVER_NORM_DATA_05FDR"

#' @rdname ATAC_NORM_DATA_05FDR
"ATAC_BAT_NORM_DATA_05FDR"

#' @rdname ATAC_NORM_DATA_05FDR
"ATAC_WATSC_NORM_DATA_05FDR"


## RRBS sample-level data ####

#' @title RRBS raw data
#' @description RRBS raw read counts; created by loading all bismark files into a single data object.
#' @format A [edgeR::DGEList-class] object
#' @details 
#'   Raw METHYL data are only available via download from Google Cloud Storage. 
#'   For example, <https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/METHYL_BAT_RAW_DATA.rda> 
#'   is the file for brown adipose tissue (BAT) data. You can change the name of the file to specify other tissues including:
#'   HEART, HIPPOC, KIDNEY, LIVER, LUNG, SKMGN (gastrocnemius skeletal muscle), and WATSC (subcutaneous white adipose tissue).
#'   You can also use [MotrpacRatTraining6mo::get_rdata_from_url()] 
#'   to download and return raw METHYL data.
#'   For more details about these files see the readme of this repository at 
#'   <https://github.com/MoTrPAC/MotrpacRatTraining6mo/blob/main/README.md>. 
#'   
#'   Unlike the [METHYL_RAW_COUNTS] data, these objects were not filtered to remove low-count features.
#'   
#'   Reads were demultiplexed with bcl2fastq (version 2.20) using options 
#'   \code{--use-bases-mask Y*,I8Y*,I*,Y* --mask-short-adapter-reads 0 --minimum-trimmed-read-length 0} 
#'   (Illumina, San Diego, CA, USA), and UMIs in the index FASTQ files were attached to the read FASTQ files. 
#'   The regular 5' and 3' adapters were trimmed with TrimGalore (v1.18), and the diversity adapter 
#'   that is about 0 to 3 bases of RDD (R={A or G} and D={A, G, or T}) that is added before YGG 
#'   (Y={C or T} depending on the methylation) from the YGG MspI cut signature was trimmed with 
#'   the NuGEN script "trimRRBSdiversityAdaptCustomers.py" (<https://github.com/nugentechnologies/NuMetRRBS>).  
#'   FastQC (v0.11.8) was used to generate pre-alignment QC metrics3. Bismark (v0.20.0) was used to index 
#'   and align reads to release 96 of the Ensembl Rattus norvegicus (rn6) genome and gene annotation. 
#'   As the lambda genome was spiked into each sample to determine the bisulfite conversion efficiency, 
#'   the lambda genome (GenBank: J02459.1) was also indexed. Default parameters were used for Bismark’s 
#'   bismark_genome_preparation in the alignment step. Bismark output BAM files were first formatted 
#'   using a custom script; Bismark’s \code{deduplicate_bismark -p --barcode} was used 
#'   to remove PCR duplicates from the bam files; and Bismark’s \code{bismark_methylation_extractor 
#'   --comprehensive --bedgraph} was used to quantify methylated and unmethylated coverages for all the 
#'   CpG sites. Bowtie 2 (v2.3.4.3) was used to index and align reads to globin, rRNA, and phix sequences 
#'   in order to quantify the percent of reads that mapped to these contaminants and spike-ins5. 
#'   SAMtools (v1.3.1) was used to compute mapping percentages to different chromosomes6. 
#'   UMIs were used to accurately quantify PCR duplicates with NuGEN’s "nodup.py" script 
#'   (<https://github.com/tecangenomics/nudup>). QC metrics from every stage of the quantification pipeline 
#'   were compiled, in part with multiQC (v1.6)7. The openWDL-based implementation of the RRBS pipeline on 
#'   Google Cloud Platform is available on GitHub (<https://github.com/MoTrPAC/motrpac-rrbs-pipeline>).
#' @name METHYL_RAW_DATA
NULL

#' @title RRBS raw counts
#' @description RRBS read counts data after filtering for CpG sites with methylation coverage of >=10 in all samples; used as input for the site clustering pipeline.
#' @format A data frame of sites as rows. Each sample has two columns: one for the methylated counts ("Me"), and nother for the unmethylated counts ("Un").
#' @details 
#'   Unfiltered METHYL sample-level data are only available via download from Google Cloud Storage. 
#'   For example, <https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/METHYL_BAT_RAW_COUNTS.rda> 
#'   is the file for brown adipose tissue (BAT) data. You can change the name of the file to specify other tissues including:
#'   HEART, HIPPOC, KIDNEY, LIVER, LUNG, SKMGN (gastrocnemius skeletal muscle), and WATSC (subcutaneous white adipose tissue).
#'   You can also use [MotrpacRatTraining6mo::load_sample_data()] or [MotrpacRatTraining6mo::get_rdata_from_url()] 
#'   to download raw and normalized sample-level data for ATAC and METHYL.
#'   For more details about these files see the readme of this repository at 
#'   <https://github.com/MoTrPAC/MotrpacRatTraining6mo/blob/main/README.md>. 
#'   
#'   Reads were demultiplexed with bcl2fastq (version 2.20) using options 
#'   \code{--use-bases-mask Y*,I8Y*,I*,Y* --mask-short-adapter-reads 0 --minimum-trimmed-read-length 0} 
#'   (Illumina, San Diego, CA, USA), and UMIs in the index FASTQ files were attached to the read FASTQ files. 
#'   The regular 5' and 3' adapters were trimmed with TrimGalore (v1.18), and the diversity adapter 
#'   that is about 0 to 3 bases of RDD (R={A or G} and D={A, G, or T}) that is added before YGG 
#'   (Y={C or T} depending on the methylation) from the YGG MspI cut signature was trimmed with 
#'   the NuGEN script "trimRRBSdiversityAdaptCustomers.py" (<https://github.com/nugentechnologies/NuMetRRBS>).  
#'   FastQC (v0.11.8) was used to generate pre-alignment QC metrics3. Bismark (v0.20.0) was used to index 
#'   and align reads to release 96 of the Ensembl Rattus norvegicus (rn6) genome and gene annotation. 
#'   As the lambda genome was spiked into each sample to determine the bisulfite conversion efficiency, 
#'   the lambda genome (GenBank: J02459.1) was also indexed. Default parameters were used for Bismark’s 
#'   bismark_genome_preparation in the alignment step. Bismark output BAM files were first formatted 
#'   using a custom script; Bismark’s \code{deduplicate_bismark -p --barcode} was used 
#'   to remove PCR duplicates from the bam files; and Bismark’s \code{bismark_methylation_extractor 
#'   --comprehensive --bedgraph} was used to quantify methylated and unmethylated coverages for all the 
#'   CpG sites. Bowtie 2 (v2.3.4.3) was used to index and align reads to globin, rRNA, and phix sequences 
#'   in order to quantify the percent of reads that mapped to these contaminants and spike-ins5. 
#'   SAMtools (v1.3.1) was used to compute mapping percentages to different chromosomes6. 
#'   UMIs were used to accurately quantify PCR duplicates with NuGEN’s "nodup.py" script 
#'   (<https://github.com/tecangenomics/nudup>). QC metrics from every stage of the quantification pipeline 
#'   were compiled, in part with multiQC (v1.6)7. The openWDL-based implementation of the RRBS pipeline on 
#'   Google Cloud Platform is available on GitHub (<https://github.com/MoTrPAC/motrpac-rrbs-pipeline>).
#'   
#' @name METHYL_RAW_COUNTS
NULL


#' @title Normalized DNA methylation data
#' @description Normalized DNA methylation (METHYL) data used for visualization.
#'   Only for training-regulated features at 5% IHW FDR. For sample-level data for \emph{all} features, 
#'   see [METHYL_RAW_DATA], [METHYL_RAW_COUNTS], and [METHYL_NORM_DATA].
#' @format A data frame with CpG sites in rows (\code{feature_ID}) and samples in columns (\code{viallabel})
#' @details 
#'   Unfiltered METHYL sample-level data are only available via download from Google Cloud Storage. 
#'   For example, <https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/METHYL_BAT_NORM_DATA.rda> 
#'   is the file for brown adipose tissue (BAT) data. You can change the name of the file to specify other tissues including:
#'   HEART, HIPPOC, KIDNEY, LIVER, LUNG, SKMGN (gastrocnemius skeletal muscle), and WATSC (subcutaneous white adipose tissue).
#'   You can also use [MotrpacRatTraining6mo::load_sample_data()] or [MotrpacRatTraining6mo::get_rdata_from_url()] 
#'   to download raw and normalized sample-level data for ATAC and METHYL.
#'   For more details about these files see the readme of this repository at 
#'   <https://github.com/MoTrPAC/MotrpacRatTraining6mo/blob/main/README.md>. 
#' 
#'   Only CpG sites with methylation coverage of >=10 in all samples were included for downstream analysis,
#'   and normalization was performed separately in each tissue. Individual CpG sites were divided into 500 base-pair
#'   windows and were clustered using the Markov Clustering algorithm via the MCL R package (Jager, 2015). To apply MCL,
#'   for each 500 base-pair window an undirected graph was constructed, linking individual sites if their correlation
#'   was >=0.7. MCL was chosen for this task as it: (1) determines the number of clusters internally, (2) identifies
#'   homogeneous clusters, and (3) keeps single sites that are not correlated with either sites as singletons (clusters of size one).
#'   The resulting sites/clusters were used as input for normalization and differential analysis with edgeR (Robinson et al., 2010).
#'   To generate this normalized sample-level data, the methylation coverages of filtered sites/clusters were first log2-transformed,
#'   and normalization was performed using [preprocessCore::normalize.quantiles.robust()] (Bolstad, 2021).
#'   
#'   After performing differential analysis, training-regulated features were selected at 5% IHW FDR. The normalized data were
#'   filtered down to these features and provided here. 
#'   
#'   For the subset of normalized data corresponding to training-regulated features at 5% IHW FDR, see [METHYL_NORM_DATA_05FDR].
#'   
#' @source \code{gs://motrpac-data-freeze-pass/pass1b-06/v1.1/analysis/epigenomics/epigen-rrbs/normalized-data/*normalized-log-M-window.txt}
#' @name METHYL_NORM_DATA
NULL


#' @title Normalized DNA methylation data for training-regulated features
#' @description Normalized DNA methylation (METHYL) data used for visualization
#' @format A data frame with CpG sites in rows (\code{feature_ID}) and samples in columns (\code{viallabel})
#' @details Only CpG sites with methylation coverage of >=10 in all samples were included for downstream analysis,
#'   and normalization was performed separately in each tissue. Individual CpG sites were divided into 500 base-pair
#'   windows and were clustered using the Markov Clustering algorithm via the MCL R package (Jager, 2015). To apply MCL,
#'   for each 500 base-pair window an undirected graph was constructed, linking individual sites if their correlation
#'   was >=0.7. MCL was chosen for this task as it: (1) determines the number of clusters internally, (2) identifies
#'   homogeneous clusters, and (3) keeps single sites that are not correlated with either sites as singletons (clusters of size one).
#'   The resulting sites/clusters were used as input for normalization and differential analysis with edgeR (Robinson et al., 2010).
#'   To generate this normalized sample-level data, the methylation coverages of filtered sites/clusters were first log2-transformed,
#'   and normalization was performed using [preprocessCore::normalize.quantiles.robust()] (Bolstad, 2021).
#'   
#'   For the full set of normalized sample-level data, see [METHYL_NORM_DATA]. 
#'   
#' @source \code{gs://motrpac-data-freeze-pass/pass1b-06/v1.1/analysis/epigenomics/epigen-rrbs/normalized-data/*normalized-log-M-window.txt}
#' @name METHYL_NORM_DATA_05FDR
"METHYL_HIPPOC_NORM_DATA_05FDR"

#' @rdname METHYL_NORM_DATA_05FDR
"METHYL_SKMGN_NORM_DATA_05FDR"

#' @rdname METHYL_NORM_DATA_05FDR
"METHYL_HEART_NORM_DATA_05FDR"

#' @rdname METHYL_NORM_DATA_05FDR
"METHYL_KIDNEY_NORM_DATA_05FDR"

#' @rdname METHYL_NORM_DATA_05FDR
"METHYL_LUNG_NORM_DATA_05FDR"

#' @rdname METHYL_NORM_DATA_05FDR
"METHYL_LIVER_NORM_DATA_05FDR"

#' @rdname METHYL_NORM_DATA_05FDR
"METHYL_BAT_NORM_DATA_05FDR"

#' @rdname METHYL_NORM_DATA_05FDR
"METHYL_WATSC_NORM_DATA_05FDR"


## Proteomics sample-level data ####

#' @title Normalized protein expression data
#' @description Median-MAD normalized protein expression (PROT) data used for visualization and differential analysis
#' @format A data frame with protein IDs in rows (\code{feature_ID}) and samples in columns (\code{viallabel}) 
#' @details Log2 TMT ratios to the common reference were used as quantitative values for all proteomics features 
#'   (proteins, phosphosites, acetylsite, and ubiquitylsites). Proteomics features not fully quantified in at least 
#'   two plexes within a tissue and non-rat contaminants were removed. Log2 TMT ratios were sample-normalized by 
#'   median-centering and mean absolute deviation scaling. Plex batch effects were removed using linear models 
#'   implemented by the [limma::removeBatchEffect()].  
#' @source \code{gs://motrpac-data-freeze-pass/pass1b-06/v1.1/analysis/proteomics-untargeted/prot-pr/normalized-data/*med-mad-normalized-logratio.txt}
#' @name PROT_NORM_DATA
"PROT_CORTEX_NORM_DATA"

#' @rdname PROT_NORM_DATA
"PROT_SKMGN_NORM_DATA"

#' @rdname PROT_NORM_DATA
"PROT_HEART_NORM_DATA"

#' @rdname PROT_NORM_DATA
"PROT_KIDNEY_NORM_DATA"

#' @rdname PROT_NORM_DATA
"PROT_LUNG_NORM_DATA"

#' @rdname PROT_NORM_DATA
"PROT_LIVER_NORM_DATA"

#' @rdname PROT_NORM_DATA
"PROT_WATSC_NORM_DATA"


#' @title Normalized protein acetylation data
#' @description Median-MAD normalized protein acetylation (ACETYL) data used for visualization and differential analysis
#' @format A data frame with acetylsites in rows (\code{feature_ID}) and samples in columns (\code{viallabel}) 
#' @details Log2 TMT ratios to the common reference were used as quantitative values for all proteomics features 
#'   (proteins, phosphosites, acetylsite, and ubiquitylsites). Proteomics features not fully quantified in at least 
#'   two plexes within a tissue and non-rat contaminants were removed. Log2 TMT ratios were sample-normalized by 
#'   median-centering and mean absolute deviation scaling. Plex batch effects were removed using linear models 
#'   implemented by the [limma::removeBatchEffect()].  
#' @source \code{gs://motrpac-data-freeze-pass/pass1b-06/v1.1/analysis/proteomics-untargeted/prot-ac/normalized-data/*med-mad-normalized-logratio.txt}
#' @name ACETYL_NORM_DATA
"ACETYL_HEART_NORM_DATA"

#' @rdname ACETYL_NORM_DATA
"ACETYL_LIVER_NORM_DATA"


#' @title Normalized protein phosphorylation data
#' @description Median-MAD normalized protein phosphorylation (PHOSPHO) data used for visualization and differential analysis
#' @format A data frame with phosphosites in rows (\code{feature_ID}) and samples in columns (\code{viallabel}) 
#' @details Log2 TMT ratios to the common reference were used as quantitative values for all proteomics features 
#'   (proteins, phosphosites, acetylsite, and ubiquitylsites). Proteomics features not fully quantified in at least 
#'   two plexes within a tissue and non-rat contaminants were removed. Log2 TMT ratios were sample-normalized by 
#'   median-centering and mean absolute deviation scaling. Plex batch effects were removed using linear models 
#'   implemented by the [limma::removeBatchEffect()].  
#' @source \code{gs://motrpac-data-freeze-pass/pass1b-06/v1.1/analysis/proteomics-untargeted/prot-ph/normalized-data/*med-mad-normalized-logratio.txt}
#' @name PHOSPHO_NORM_DATA
"PHOSPHO_CORTEX_NORM_DATA"

#' @rdname PHOSPHO_NORM_DATA
"PHOSPHO_SKMGN_NORM_DATA"

#' @rdname PHOSPHO_NORM_DATA
"PHOSPHO_HEART_NORM_DATA"

#' @rdname PHOSPHO_NORM_DATA
"PHOSPHO_KIDNEY_NORM_DATA"

#' @rdname PHOSPHO_NORM_DATA
"PHOSPHO_LUNG_NORM_DATA"

#' @rdname PHOSPHO_NORM_DATA
"PHOSPHO_LIVER_NORM_DATA"

#' @rdname PHOSPHO_NORM_DATA
"PHOSPHO_WATSC_NORM_DATA"


#' @title Normalized protein ubiquitynation data
#' @description Median-MAD normalized protein ubiquitynation (UBIQ) data used for visualization and differential analysis
#' @format A data frame with phosphosites in rows (\code{feature_ID}) and samples in columns (\code{viallabel}) 
#' @details Log2 TMT ratios to the common reference were used as quantitative values for all proteomics features 
#'   (proteins, phosphosites, acetylsite, and ubiquitylsites). Proteomics features not fully quantified in at least 
#'   two plexes within a tissue and non-rat contaminants were removed. Log2 TMT ratios were sample-normalized by 
#'   median-centering and mean absolute deviation scaling. Plex batch effects were removed using linear models 
#'   implemented by the [limma::removeBatchEffect()]. The ubiquitynation datasets were corrected for changes in protein abundances 
#'   by fitting a global linear model between the ubiquitylsite and the cognate protein and extracting the residuals. 
#' @source \code{gs://motrpac-data-freeze-pass/pass1b-06/v1.1/analysis/proteomics-untargeted/prot-ph/normalized-data/*med-mad-normalized-protein-corrected-logratio.txt}
#' @name UBIQ_NORM_DATA
"UBIQ_HEART_NORM_DATA"

#' @rdname UBIQ_NORM_DATA
"UBIQ_LIVER_NORM_DATA"


## Immunoassay sample-level data ####

#' @title Processed immunoassay data used for differential analysis
#' @description Normalized, imputed, and filtered multiplexed immunoassay data used for differential analysis  
#' @format A nested list of data frames
#' @details 
#'   IMMUNO sample-level data is in a different format than sample-level data for other assays/omes. 
#'   Extract data from a panel and tissue using \code{\link{IMMUNO_NORM_DATA_NESTED}[[dataset]][[tissue]]}, where \code{dataset} 
#'   is one of "ADIPONECTIN", "SERPIN-E", "rat-mag27plex", "rat-metabolic", "rat-myokine", "rat-pituitary", and 
#'   \code{tissue} is a tissue abbreviation (see [TISSUE_ABBREV]). Samples (vial labels) are in rows, and analytes are
#'   in columns. Column names, barring the first "viallabel" column, correspond to \code{feature_ID}s. 
#'   
#'   Raw mean fluorescent intensities (MFIs) were log2-transformed, and measurements corresponding to wells with less than 20 beads were removed. 
#'   For remaining measurements in each panel and tissue, samples with more than 50\% missing values (i.e., due to low bead count) were removed; 
#'   features with at least two missing values for a single experimental group (e.g., males trained for 2 weeks) were removed 
#'   (this affected four colon analytes and one spleen analyte, all in the rat-mag27plex panel). Remaining missing values were 
#'   imputed with k-nearest neighbors (k=5 features). Within each panel, tissue, and analyte, we calculated the mean and standard deviation 
#'   and removed outlying measurements more than 4 standard deviations away from the mean. This is the version of the data provided by this object
#'   and used for differential analysis. 
#'   
#'   SERPIN-E and ADIPONECTIN were both originally in a panel called rat-adipokine. They 
#'   were split into their own "panels" because this panel was run with two dilutions, 
#'   and one was optimal for SERPIN-E while the other was optimal for ADIPONECTIN.
#'   Hence, in some places \code{dataset} refers to "rat-adipokine" while in other places, 
#'   like here, \code{dataset} refers to "SERPIN-E" or "ADIPONECTIN", among others. 
#' @source \code{gs://mawg-data/pass1b-06/immunoassay/data/release/pass1b-06*_mfi-log2-filt-imputed-na-outliers.txt} 
"IMMUNO_NORM_DATA_NESTED"


#' @title Combined immunoassay data used for visualization
#' @description Normalized, imputed, and filtered multiplexed immunoassay data used for visualization.
#'   Data are equivalent to the data provided in [IMMUNO_NORM_DATA_NESTED]. [IMMUNO_NORM_DATA_NESTED] is compatible with the 
#'   differential analysis functions while this format is compatible with visualization functions. 
#' @format A data frame with analytes in rows and participant IDs (PIDs) in columns
#' @source \code{gs://mawg-data/pass1b-06/immunoassay/data/release/pass1b-06*_mfi-log2-filt-imputed-na-outliers.txt} 
"IMMUNO_NORM_DATA_FLAT"


## Metabolomics sample-level data ####

#' @title Processed metabolomics data used for differential analysis
#' @description Combined sample-level data organized by metabolomics platforms and tissue used for differential analysis
#' @format A nested list of data frames
#' @details 
#'   METAB sample-level data is in a different format than sample-level data for other assays/omes. 
#'   Extract data from a panel and tissue using \code{METAB_NORM_DATA_NESTED[[platform]][[tissue]]}, where \code{platform} 
#'   is one of "metab-t-amines", "metab-t-acoa", "metab-t-nuc", "metab-t-oxylipneg", "metab-t-ka", "metab-t-etamidpos", "metab-t-tca", 
#'   "metab-u-lrppos", "metab-u-lrpneg", "metab-u-hilicpos", "metab-u-rppos", "metab-u-rpneg", "metab-u-ionpneg" (see [MotrpacBicQC::assay_codes] for details) and 
#'   \code{tissue} is a tissue abbreviation (see [TISSUE_ABBREV]). Samples (vial labels) are in columns, and metabolites (feature IDs) are
#'   in rows. Feature IDs in the row names in these tables correspond to \code{\link{METAB_FEATURE_ID_MAP}$feature_ID_sample_data}.
#'   
#'   Not all data sets were processed similarly. The processing for each tissue/platform data follows this criteria:  
#'   \itemize{
#'     \item Metabolomics targeted 't' platform with >12 features: KNN-imputed and log2-transformed
#'     \item Metabolomics targeted 't' platform with <= 12 features: Log2-transformed (no imputation) 
#'     \item Metabolomics untargeted 'u' platform: sample-centered, KNN-imputed, log2-transformed
#'   }
#'   
"METAB_NORM_DATA_NESTED"


#' @title Combined metabolomics data used for visualization
#' @description Combined sample-level metabolomics data used for visualization.
#'   Data are equivalent to the data provided in [METAB_NORM_DATA_NESTED]. [METAB_NORM_DATA_NESTED] is compatible with the 
#'   differential analysis functions while this format is compatible with visualization functions. 
#' @format A data frame with metabolites in rows and participant IDs (PIDs) in columns 
"METAB_NORM_DATA_FLAT"


## Differential analysis ####

### Ome-specific metadata ####

#' @title RNA-seq metadata and QC  
#' @description RNA-seq experimental and quantification QC metrics for transcriptomic (TRNSCRPT) data
#' @format A data frame with 935 rows and 82 variables:
#' \describe{
#'   \item{\code{viallabel}}{`r viallabel()`}
#'   \item{\code{vial_label}}{double, sample identifier}
#'   \item{\code{2D_barcode}}{double, sample barcode}
#'   \item{\code{Species}}{character, species}
#'   \item{\code{BID}}{integer, biospecimen ID}
#'   \item{\code{PID}}{double, participant ID, one per animal}
#'   \item{\code{Tissue}}{character, tissue description}
#'   \item{\code{Sample_category}}{character, study sample ("study") or reference standard ("ref)}
#'   \item{\code{GET_site}}{character, which Genomics, Epigenomics, and Transcriptomics (GET) 
#'     site performed the assay, "Stanford" or "MSSM" (Icahn School of Medicine at Mount Sinai)}
#'   \item{\code{RNA_extr_plate_ID}}{character, RNA extraction plate ID}
#'   \item{\code{RNA_extr_date}}{character, RNA extraction date}
#'   \item{\code{RNA_extr_conc}}{double, RNA concentration (ng/uL)}
#'   \item{\code{RIN}}{double, RNA Integrity Number}
#'   \item{\code{r_260_280}}{double, 260/280 ratio}
#'   \item{\code{r_260_230}}{double 260/230 ratio}
#'   \item{\code{Lib_prep_date}}{character, library preparation date in MM/DD/YYYY format}
#'   \item{\code{Lib_RNA_conc}}{double, RNA concentration used for library prep (ng/uL)}
#'   \item{\code{Lib_RNA_vol}}{integer, RNA volume used for library prep (uL)}
#'   \item{\code{Lib_robot}}{character, robot used for library prep}
#'   \item{\code{Lib_vendor}}{character, library prep vendor}
#'   \item{\code{Lib_type}}{character, library prep type}
#'   \item{\code{Lib_kit_id}}{character, library prep kit ID}
#'   \item{\code{Lib_batch_ID}}{character, library prep batch ID that distinguished different sample processing batches}
#'   \item{\code{Lib_barcode_well}}{character, well}
#'   \item{\code{Lib_index_1}}{character, i7 index}
#'   \item{\code{Lib_index_2}}{character, i5 index}
#'   \item{\code{Lib_adapter_1}}{character, Truseq I7 index with 16bp index}
#'   \item{\code{Lib_adapter_2}}{character, Truseq I5 index with 8bp index}
#'   \item{\code{Lib_UMI_cycle_num}}{integer, number of bases of UMI}
#'   \item{\code{Lib_adapter_size}}{integer, total size of the two adapters}
#'   \item{\code{Lib_frag_size}}{integer, average library fragment size (bp)}
#'   \item{\code{Lib_DNA_conc}}{double, DNA concentration of original stock of the library (ng/uL)}
#'   \item{\code{Lib_molarity}}{double, library molarity (nM)}
#'   \item{\code{Seq_platform}}{character, sequencing platform}
#'   \item{\code{Seq_date}}{integer, sequencing date, YYMMDD format}
#'   \item{\code{Seq_machine_ID}}{character, serial number of the sequencer}
#'   \item{\code{Seq_flowcell_ID}}{character, flow cell ID}
#'   \item{\code{Seq_flowcell_run}}{integer, flow cell run}
#'   \item{\code{Seq_flowcell_lane}}{character, flow cell lane}
#'   \item{\code{Seq_flowcell_type}}{character, flow cell type, e.g., S4}
#'   \item{\code{Seq_length}}{integer, read length}
#'   \item{\code{Seq_end_type}}{integer, 1=single-end, 2=paired-end}
#'   \item{\code{Phase}}{character, study phase, "PASS1B-06"}
#'   \item{\code{Seq_batch}}{character, unique identifier for sequencing batch}
#'   \item{\code{reads_raw}}{double, number of read pairs in the raw FASTQ}
#'   \item{\code{pct_adapter_detected}}{double, percent of reads with adapter detected}
#'   \item{\code{pct_trimmed}}{double, percent of reads that were trimmed}
#'   \item{\code{pct_trimmed_bases}}{double, percent of bases that were trimmed}
#'   \item{\code{reads}}{double, number of read pairs in the trimmed FASTQ files}
#'   \item{\code{pct_GC}}{double, percent GC content in trimmed FASTQ files}
#'   \item{\code{pct_dup_sequence}}{double, percent of duplicated sequences in trimmed FASTQ files}
#'   \item{\code{pct_rRNA}}{double, percent of rRNA reads in trimmed FASTQ files}
#'   \item{\code{pct_globin}}{double, percent of globin reads in trimmed FASTQ files}
#'   \item{\code{pct_phix}}{double, percent of phix reads in trimmed FASTQ files}
#'   \item{\code{pct_picard_dup}}{double, PCR duplication assessed by Picard’s tool MarkDuplicate}
#'   \item{\code{pct_umi_dup}}{double, PCR duplication rate assessed using UMIs (Unique Molecular Identifiers)}
#'   \item{\code{avg_input_read_length}}{double, average input read length}
#'   \item{\code{uniquely_mapped}}{double, number of uniquely mapped reads}
#'   \item{\code{pct_uniquely_mapped}}{double, percent of uniquely mapped reads}
#'   \item{\code{avg_mapped_read_length}}{double, average input mapped length}
#'   \item{\code{num_splices}}{double, number of splices}
#'   \item{\code{num_annotated_splices}}{double, number of annotated splices}
#'   \item{\code{num_GTAG_splices}}{double, number of GT/AG and CT/AC splices}
#'   \item{\code{num_GCAG_splices}}{double, number of GC/AG and CT/GC splices}
#'   \item{\code{num_ATAC_splices}}{double, number of AT/AC and GT/TA splices}
#'   \item{\code{num_noncanonical_splices}}{double, number of non-canonical splices}
#'   \item{\code{pct_multimapped}}{double, percent of reads that multimapped}
#'   \item{\code{pct_multimapped_toomany}}{double, percent of reads that multimapped too many times}
#'   \item{\code{pct_unmapped_mismatches}}{double, percent of unmapped reads due to mismatches}
#'   \item{\code{pct_unmapped_tooshort}}{double, percent of unmapped reads due to shortness}
#'   \item{\code{pct_unmapped_other}}{double, percent of unmapped reads for other reason}
#'   \item{\code{pct_chimeric}}{double, percent chimeric reads}
#'   \item{\code{pct_chrX}}{double, percent of reads mapped to chromosome X}
#'   \item{\code{pct_chrY}}{double, percent of reads mapped to chromosome Y}
#'   \item{\code{pct_chrM}}{double, percent of reads mapped to the mitochondrial genome}
#'   \item{\code{pct_chrAuto}}{double, percent of reads mapped to autosomal chromosomes}
#'   \item{\code{pct_contig}}{double, percent of reads mapped to contigs}
#'   \item{\code{pct_coding}}{double, percent of bases mapped to coding}
#'   \item{\code{pct_utr}}{double, percent of bases mapped to untranslated region}
#'   \item{\code{pct_intronic}}{double, percent of bases mapped to introns}
#'   \item{\code{pct_intergenic}}{double, percent of bases mapped to intergenic}
#'   \item{\code{pct_mrna}}{double, percent of bases mapped to mRNA}
#'   \item{\code{median_5_3_bias}}{double, median 5' to 3' bias} 
#' }
#' @source <gs://motrpac-data-freeze-pass/pass1b-06/v1.1/results/transcriptomics/qa-qc/motrpac_pass1b-06_transcript-rna-seq_qa-qc-metrics.csv>
"TRNSCRPT_META"


#' @title ATAC-seq metadata and QC 
#' @description ATAC-seq experimental and quantification QC metrics for chromatin accessibility (ATAC) data 
#' @format A data frame with 520 rows and 106 variables:
#' \describe{
#'   \item{\code{viallabel}}{`r viallabel()`}
#'   \item{\code{general.description}}{character, pipeline workflow description}
#'   \item{\code{replicate}}{character, replicate in the pipeline workflow}
#'   \item{\code{general.date}}{double, date the pipeline workflow was run}
#'   \item{\code{general.title}}{character, pipeline workflow title}
#'   \item{\code{general.pipeline_ver}}{character, ENCODE ATAC-seq pipeline version}
#'   \item{\code{general.pipeline_type}}{character, "atac"}
#'   \item{\code{general.genome}}{character, reference genome}
#'   \item{\code{general.aligner}}{character, read aligner}
#'   \item{\code{general.peak_caller}}{character, peak caller}
#'   \item{\code{general.seq_endedness.paired_end}}{logical, are the reads paired-ended}
#'   \item{\code{replication.num_peaks.num_peaks}}{integer, number of replication peaks (max 300000)}
#'   \item{\code{peak_stat.peak_region_size.min_size}}{integer, minimum peak width}
#'   \item{\code{peak_stat.peak_region_size.25_pct}}{integer, 25th percentile of peak width}
#'   \item{\code{peak_stat.peak_region_size.50_pct}}{integer, 50th percentile of peak width}
#'   \item{\code{peak_stat.peak_region_size.75_pct}}{integer, 75th percentile of peak width}
#'   \item{\code{peak_stat.peak_region_size.max_size}}{integer, max peak width}
#'   \item{\code{peak_stat.peak_region_size.mean}}{double, mean peak width}
#'   \item{\code{peak_enrich.frac_reads_in_peaks.macs2.frip}}{double, fraction of reads in MACS2 peaks}
#'   \item{\code{align.samstat.total_reads}}{integer, total number of alignments, including multimappers}
#'   \item{\code{align.samstat.mapped_reads}}{integer, total number of mapped reads}
#'   \item{\code{align.samstat.pct_mapped_reads}}{double, percent of reads that mapped}
#'   \item{\code{align.samstat.paired_reads}}{integer, number of paired reads}
#'   \item{\code{align.samstat.read1}}{integer, number of read 1 reads}
#'   \item{\code{align.samstat.read2}}{integer, number of read 2 reads}
#'   \item{\code{align.samstat.properly_paired_reads}}{integer, number of properly paired reads}
#'   \item{\code{align.samstat.pct_properly_paired_reads}}{double, percent of reads that were properly paired}
#'   \item{\code{align.samstat.with_itself}}{integer, number of reads paired with its pair}
#'   \item{\code{align.samstat.singletons}}{integer, number of singleton reads}
#'   \item{\code{align.samstat.pct_singletons}}{double, percent of reads that were singleton}
#'   \item{\code{align.samstat.diff_chroms}}{integer}
#'   \item{\code{align.dup.paired_reads}}{integer, number of paired reads for duplication step}
#'   \item{\code{align.dup.paired_duplicate_reads}}{integer, number of duplicate paired reads}
#'   \item{\code{align.dup.paired_optical_duplicate_reads}}{integer, number of optical duplicate paired reads}
#'   \item{\code{align.dup.pct_duplicate_reads}}{double, percent of reads that are duplicate}
#'   \item{\code{align.frac_mito.non_mito_reads}}{integer, percent of reads that align to non-mitochondrial DNA}
#'   \item{\code{align.frac_mito.mito_reads}}{integer, number of reads that align to mitochondrial DNA}
#'   \item{\code{align.frac_mito.frac_mito_reads}}{double, fraction of reads that align to mitochondrial DNA}
#'   \item{\code{align.nodup_samstat.total_reads}}{integer, number of reads after applying all filters}
#'   \item{\code{align.nodup_samstat.mapped_reads}}{integer, number of mapped reads after applying all filters}
#'   \item{\code{align.nodup_samstat.paired_reads}}{integer, number of paired reads after applying all filters}
#'   \item{\code{align.nodup_samstat.read1}}{integer, number of read 1 reads after applying all filters}
#'   \item{\code{align.nodup_samstat.read2}}{integer, number of read 2 reads after applying all filters}
#'   \item{\code{align.nodup_samstat.properly_paired_reads}}{integer, number of properly paired reads after applying all filters}
#'   \item{\code{align.nodup_samstat.with_itself}}{integer, number of reads paired with its pair after applying all filters}
#'   \item{\code{align.frag_len_stat.frac_reads_in_nfr}}{double, fraction of reads in nucelosome-free region. Should be a value greater than 0.4.}
#'   \item{\code{align.frag_len_stat.frac_reads_in_nfr_qc_pass}}{logical, does \code{align.frag_len_stat.frac_reads_in_nfr} pass the cutoff?}
#'   \item{\code{align.frag_len_stat.frac_reads_in_nfr_qc_reason}}{character, reason for \code{align.frag_len_stat.frac_reads_in_nfr_qc_pass}}
#'   \item{\code{align.frag_len_stat.nfr_over_mono_nuc_reads}}{double, reads in nucleosome-free-region 
#'     versus reads in mononucleosomal peak. Should be a value greater than 2.5.}
#'   \item{\code{align.frag_len_stat.nfr_over_mono_nuc_reads_qc_pass}}{logical, does \code{align.frag_len_stat.frac_reads_in_nfr_qc_pass} pass the cutoff?}
#'   \item{\code{align.frag_len_stat.nfr_over_mono_nuc_reads_qc_reason}}{character, reason for \code{align.frag_len_stat.nfr_over_mono_nuc_reads_qc_pass}}
#'   \item{\code{align.frag_len_stat.nfr_peak_exists}}{logical, does a nucleosome-free peak exist?}
#'   \item{\code{align.frag_len_stat.mono_nuc_peak_exists}}{logical, does a mono-nucleosomal peak exist?}
#'   \item{\code{align.frag_len_stat.di_nuc_peak_exists}}{logical, does a di-nucleosomal peak exist?}
#'   \item{\code{lib_complexity.lib_complexity.total_fragments}}{integer, total number of fragments}
#'   \item{\code{lib_complexity.lib_complexity.distinct_fragments}}{integer, number of distinct fragments}
#'   \item{\code{lib_complexity.lib_complexity.positions_with_one_read}}{integer, number of positions with one read}
#'   \item{\code{lib_complexity.lib_complexity.NRF}}{double, non-reduandant fraction. Measure of library complexity. Ideally >0.9}
#'   \item{\code{lib_complexity.lib_complexity.PBC1}}{double, PCR bottlenecking coefficient 1. Measure of library complexity. Ideally >0.9}
#'   \item{\code{lib_complexity.lib_complexity.PBC2}}{double PCR bottlenecking coefficient 2. Measure of library complexity. Ideally >3}
#'   \item{\code{align_enrich.tss_enrich.tss_enrich}}{double, transcription start site enrichment of peaks}
#'   \item{\code{2D_barcode}}{double, sample barcode}
#'   \item{\code{Tissue}}{character, tissue description}
#'   \item{\code{Species}}{character, species}
#'   \item{\code{Sample_category}}{character, study sample ("study") or reference standard ("ref)}
#'   \item{\code{GET_site}}{character, which Genomics, Epigenomics, and Transcriptomics (GET) 
#'     site performed the assay, "Stanford" or "MSSM" (Icahn School of Medicine at Mount Sinai)}
#'   \item{\code{Sample_batch}}{integer, numeric batch number for batch in which this sample was manually processed}
#'   \item{\code{Lib_adapter_1}}{character, Adapter sequence for read 1}
#'   \item{\code{Lib_adapter_2}}{character, Adapter sequence for read 2}
#'   \item{\code{Lib_index_1}}{character, i7 index}
#'   \item{\code{Lib_index_2}}{character, i5 index}
#'   \item{\code{Nuclei_extr_date}}{character, nuclei extraction date}
#'   \item{\code{Nuclei_extr_count}}{integer, nuclei count}
#'   \item{\code{Nuclei_tagmentation}}{integer, number of nuclei used in each tagmentation reaction}
#'   \item{\code{Tagmentation_date}}{character, tagmentation date, MM/DD/YYYY format}
#'   \item{\code{Tagmentation_enzyme_cat}}{integer, catalog number of tagmentation enzyme TDE1 (Tn5)}
#'   \item{\code{Tagmentation_enzyme_lot}}{integer, lot number of tagmentation enzyme TDE1 (Tn5)}
#'   \item{\code{Tagmentation_buffer_cat}}{integer, catalog number of tagmentation buffer}
#'   \item{\code{Tagmentation_buffer_lot}}{integer, lot number of tagmentation buffer}
#'   \item{\code{Tagmentation_reaction_vol}}{integer, volume of tagmentation (uL)}
#'   \item{\code{Tagmentation_purification_kit}}{character, purification kit}
#'   \item{\code{Tagmentation_purified_DNA_vol}}{double, volume of purified DNA (uL)}
#'   \item{\code{PCR_date}}{character, PCR date, MM/DD/YYYY format}
#'   \item{\code{PCR_cycle_nr}}{integer, number of PCR cycles}
#'   \item{\code{PCR_purification_beads_ul}}{character, volume of SPRIselect beads for lower size selection}
#'   \item{\code{Lib_DNA_conc}}{double, DNA concentration for the library (ng/uL)}
#'   \item{\code{Lib_DNA_molarity}}{double, DNA molarity of library (nM)}
#'   \item{\code{Lib_frag_size}}{integer, average library fragment size}
#'   \item{\code{Lib_BA_quality}}{integer, visual inspection of the library quality with the Bioanalyzer track (1=good, 0=bad)}
#'   \item{\code{Seq_DNA_molarity}}{double, DNA molarity for sequencing (nM)}
#'   \item{\code{Seq_platform}}{character, sequencing platform}
#'   \item{\code{Seq_date}}{integer, sequencing date, YYMMDD format}
#'   \item{\code{Seq_machine_ID}}{character, serial number of the sequencer}
#'   \item{\code{Seq_flowcell_ID}}{character, flow cell ID}
#'   \item{\code{Seq_flowcell_run}}{integer, flow cell run}
#'   \item{\code{Seq_flowcell_lane}}{character, flow cell lane}
#'   \item{\code{Seq_flowcell_type}}{character, flow cell type, e.g., S4}
#'   \item{\code{Seq_length}}{double, read length}
#'   \item{\code{Seq_end_type}}{integer, 1=single-end, 2=paired-end}
#'   \item{\code{total_primary_alignments}}{integer, number of primary alignments}
#'   \item{\code{pct_chrX}}{double, number of reads mapped to chromosome X}
#'   \item{\code{pct_chrY}}{double, number of reads mapped to chromosome Y}
#'   \item{\code{pct_chrM}}{double, number of reads mapped to chromosome M}
#'   \item{\code{pct_auto}}{double, number of reads mapped to autosomal chromosomes}
#'   \item{\code{pct_contig}}{double, number of reads mapped to contigs}
#'   \item{\code{Seq_batch}}{character, unique identifier for sequencing batch} 
#' }
#' @details The [ENCODE ATAC-seq pipeline v1.7.0](https://github.com/ENCODE-DCC/atac-seq-pipeline/tree/v1.7.0) was used to quantify ATAC-seq data. 
#'   Columns with a period are QC metrics from this pipeline. Note that the ENCODE pipeline reports alignments per paired-end read, 
#'   so \code{align.samstat.total_reads} reports the number of paired-end reads that align, which corresponds to twice the number of sequenced fragments.    
#' 
#' @source <gs://motrpac-data-freeze-pass/pass1b-06/v1.1/results/epigenomics/qa-qc/motrpac_pass1b-06_epigen-atac-seq_qa-qc-metrics.csv>
"ATAC_META"


#' @title RRBS metadata and QC
#' @description RRBS experimental and quantification QC metrics for DNA methylation (METHYL) data 
#' @format A data frame with 416 rows and 75 variables:
#' \describe{
#'   \item{\code{viallabel}}{`r viallabel()`}
#'   \item{\code{vial_label}}{double, sample identifier}
#'   \item{\code{2D_barcode}}{double, sample barcode}
#'   \item{\code{Species}}{character, species}
#'   \item{\code{BID}}{integer, biospecimen ID}
#'   \item{\code{PID}}{double, participant ID, one per animal}
#'   \item{\code{Tissue}}{character, tissue description}
#'   \item{\code{Sample_category}}{character, study sample ("study") or reference standard ("ref)}
#'   \item{\code{GET_site}}{character, which Genomics, Epigenomics, and Transcriptomics (GET) 
#'     site performed the assay, "Stanford" or "MSSM" (Icahn School of Medicine at Mount Sinai)}
#'   \item{\code{DNA_extr_plate_ID}}{integer, DNA extraction plate ID}
#'   \item{\code{DNA_extr_date}}{character, DNA extraction date}
#'   \item{\code{DNA_extr_protocol}}{character, DNA extraction protocol}
#'   \item{\code{DNA_extr_robot}}{character, robot used for DNA extraction}
#'   \item{\code{DNA_conc}}{double, DNA concentration (ng/uL)}
#'   \item{\code{A280/260}}{double, 280/260 ratio}
#'   \item{\code{A260/230}}{double, 260/230 ratio}
#'   \item{\code{Lib_prep_date}}{character, library preparation date}
#'   \item{\code{Lib_DNA_mass}}{double, DNA mass used for library prep (ng)}
#'   \item{\code{Lib_DNA_vol}}{double, volume of the library (uL)}
#'   \item{\code{lambda_DNA_mass}}{double, spiked-in Lambda DNA mass (ng)}
#'   \item{\code{Lib_robot}}{character, robot used for library prep}
#'   \item{\code{Lib_kit_vendor}}{character, library prep vendor}
#'   \item{\code{Lib_kit_type}}{character, library prep kit}
#'   \item{\code{Lib_kit_ID}}{character, library kit ID}
#'   \item{\code{Lib_batch_ID}}{character, library prep batch ID}
#'   \item{\code{Lib_index_1}}{character, i7 index}
#'   \item{\code{Lib_index_2}}{logical, i5 index}
#'   \item{\code{Lib_adapter_1}}{character, Trueseq i7 index with 16 bp index}
#'   \item{\code{Lib_adapter_2}}{character, include customized Metseq primer}
#'   \item{\code{Lib_UMI_cycle_num}}{integer, number of bases of UMI}
#'   \item{\code{Lib_adapter_size}}{integer, the total size of the two adapters}
#'   \item{\code{Lib_conc}}{double, DNA concentration for the library (ng/uL)}
#'   \item{\code{Lib_frag_size}}{integer, average library fragment size}
#'   \item{\code{Lib_molarity}}{double, library molarity (nM)}
#'   \item{\code{Seq_platform}}{character, sequencing platform}
#'   \item{\code{Seq_date}}{integer, sequencing date}
#'   \item{\code{Seq_machine_ID}}{character, serial number of the sequencing machine}
#'   \item{\code{Seq_flowcell_ID}}{character, flow cell ID}
#'   \item{\code{Seq_flowcell_run}}{integer, flow cell run}
#'   \item{\code{Seq_flowcell_lane}}{character, flow cell lane}
#'   \item{\code{Seq_flowcell_type}}{character, flow cell type, e.g., "S4"}
#'   \item{\code{Seq_length}}{integer, read length}
#'   \item{\code{Seq_end_type}}{integer, 1=single-end, 2=paired-end}
#'   \item{\code{reads_raw}}{integer, number of raw read pairs}
#'   \item{\code{pct_adapter_detected}}{double, percent of reads with adapter detected}
#'   \item{\code{pct_trimmed}}{double, percent of trimmed reads from the initial trimming}
#'   \item{\code{pct_no_MSPI}}{double, percent of reads with no MSPI present among the trimmed reads}
#'   \item{\code{pct_trimmed_bases}}{double, percent of bases that were trimmed}
#'   \item{\code{pct_removed}}{double, percent of reads that were removed due to adapter trimming or no presence of MSPI}
#'   \item{\code{reads}}{integer, number of read pairs in the trimmed FASTQ files}
#'   \item{\code{pct_GC}}{double, percent GC content in the trimmed FASTQ files}
#'   \item{\code{pct_dup_sequence}}{double, percent of duplicated sequences in trimmed FASTQ files}
#'   \item{\code{pct_phix}}{double, percent of phix reads in trimmed FASTQ files}
#'   \item{\code{pct_chrX}}{double, percent of reads mapped to chromosome X}
#'   \item{\code{pct_chrY}}{double, percent of reads mapped to chromosome Y}
#'   \item{\code{pct_chrM}}{double, percent of reads mapped to the mitochondrial genome}
#'   \item{\code{pct_chrAuto}}{double, percent of reads mapped to autosomal chromosomes}
#'   \item{\code{pct_contig}}{double, percent of reads mapped to contigs}
#'   \item{\code{pct_Uniq}}{double, percent of uniquely mapped reads}
#'   \item{\code{pct_Unaligned}}{double, percent of unaligned reads}
#'   \item{\code{pct_Ambi}}{double, percent of ambiguously mapped reads}
#'   \item{\code{pct_OT}}{double, percent of mapped reads aligned to the original top stand}
#'   \item{\code{pct_OB}}{double, percent of mapped reads aligned to the original bottom stand}
#'   \item{\code{pct_CTOT}}{double, percent of mapped reads aligned to the complementary to original top strand}
#'   \item{\code{pct_CTOB}}{double, percent of mapped reads aligned to the complementary to original bottom strand}
#'   \item{\code{pct_umi_dup}}{double, PCR duplication rate assessed using UMIs (Unique Molecular Identifiers)}
#'   \item{\code{pct_CpG}}{double, global CpG methylation level based on the deduplicated data}
#'   \item{\code{pct_CHG}}{double, global CHG methylation level  based on the deduplicated data}
#'   \item{\code{pct_CHH}}{double, global CHH methylation level  based on the deduplicated data}
#'   \item{\code{lambda_pct_Uniq}}{double, percent of uniquely mapped reads to lambda}
#'   \item{\code{lambda_pct_Ambi}}{double, percent of ambiguously mapped reads to lambda}
#'   \item{\code{lambda_pct_umi_dup}}{double, PCR duplication rate assessed using UMIs (Unique Molecular Identifiers)}
#'   \item{\code{lambda_pct_CpG}}{double, global CpG methylation level based on the deduplicated data}
#'   \item{\code{lambda_pct_CHG}}{double, global CHG methylation level  based on the deduplicated data}
#'   \item{\code{lambda_pct_CHH}}{double, global CHH methylation level  based on the deduplicated data}
#'   \item{\code{Seq_batch}}{character, unique identifier for sequencing batch}
#' }
#' @source <gs://motrpac-data-freeze-pass/pass1b-06/v1.1/results/epigenomics/qa-qc/motrpac_pass1b-06_epigen-rrbs_qa-qc-metrics.csv>
"METHYL_META"


#' @title Multiplexed immunoassay metadata and QC
#' @description Multiplexed immunoassay (IMMUNO) experimental and quantification QC metrics 
#' @format A data frame with 1511 rows and 23 variables:
#' \describe{
#'   \item{\code{viallabel}}{`r viallabel()`}
#'   \item{\code{tissue}}{`r tissue()`}
#'   \item{\code{tissue_code}}{`r tissue_code()`}
#'   \item{\code{bid}}{integer, biospecimen ID}
#'   \item{\code{luminex_sample_name}}{character, sample name used by HIMC in raw data files}
#'   \item{\code{panel_name}}{character, LUMINEX panel name}
#'   \item{\code{plate_id}}{character, plate ID in format "`[date]-[tissues]_[panel_name]`"}
#'   \item{\code{ppt_type}}{character, species}
#'   \item{\code{weight_mg}}{integer, sample weight in mg (solid tissue only)}
#'   \item{\code{volume_ul}}{integer, sample volume in uL (plasma only)}
#'   \item{\code{date_extracted}}{character, date of protein extraction}
#'   \item{\code{operators}}{character, initials of experimentalists}
#'   \item{\code{protein_conc_mg_ml}}{double, protein concentration in mg/mL}
#'   \item{\code{vol_plate1_ul}}{integer, volume (uL) of sample submitted in first master plate (plate used for LUMINEX assays)}
#'   \item{\code{min_vol_plate2_ul}}{integer, volume (uL) of remaining sample submitted in second master plate (backup plate)}
#'   \item{\code{comments}}{character, comments about protein extraction}
#'   \item{\code{sample_well_position}}{character, well position in plates submitted to the HIMC}
#'   \item{\code{luminex_well}}{character, well position for LUMINEX assay}
#'   \item{\code{sex}}{`r sex()`}
#'   \item{\code{group}}{character, intervention group, one of "1w", "2w", "4w", "8w", "control"}
#'   \item{\code{log2_CHEX1}}{double, custom Assay CHEX control beads to monitor instrument performance, log-transformed}
#'   \item{\code{log2_CHEX2}}{double, custom Assay CHEX control beads to monitor application of the detection antibody, log-transformed}
#'   \item{\code{log2_CHEX3}}{double, custom Assay CHEX control beads to monitor application of the fluorescent reporter, log-transformed}
#'   \item{\code{log2_CHEX4}}{double, custom Assay CHEX control beads to monitor nonspecific binding, log-transformed} 
#' }
#' @details Protein extractions from the same sample/vial label were used for multiple
#'   panels, so metadata must be matched to unique samples using both \code{viallabel} and \code{panel_name}. 
"IMMUNO_META"


#' @title Proteomics sample-level metadata
#' @description TMT channel and plex information for each proteomics assay
#' @format A data frame with 5 variables:
#' \describe{
#'   \item{\code{tmt11_channel}}{character, TMT channel}
#'   \item{\code{tmt_plex}}{character, TMT plex}
#'   \item{\code{tissue}}{`r tissue()`}
#'   \item{\code{viallabel}}{`r viallabel()`}
#'   \item{\code{assay}}{`r assay()`} 
#' }
#' @details Sample information specific to each proteomics assay. This includes the TMT channel assigned to the 
#'   sample (\code{tmt11_channel}) and the plex group in which the sample was assigned (\code{tmt_plex}).
#' @name PROTEOMICS_META
"PHOSPHO_META"

#' @rdname PROTEOMICS_META
"PROT_META"

#' @rdname PROTEOMICS_META
"ACETYL_META"

#' @rdname PROTEOMICS_META
"UBIQ_META"


### Outliers ####

#' @title Sample outliers
#' @description Outliers excluded during differential analysis
#' @format A data frame with 27 rows and 9 variables:
#' \describe{
#'   \item{\code{viallabel}}{`r viallabel()`}
#'   \item{\code{tissue}}{`r tissue()`}
#'   \item{\code{tissue_code}}{`r tissue_code()`}
#'   \item{\code{assay}}{`r assay()`}
#'   \item{\code{assay_code}}{`r assay_code()`}
#'   \item{\code{platform}}{character, LUMINEX panel if \code{assay} is IMMUNO}
#'   \item{\code{pid}}{integer, participant ID, one per animal}
#'   \item{\code{group}}{character, combination of sex and training time point that the sample belongs to, e.g., "female_1w"}
#'   \item{\code{reason}}{character, reason(s) the sample was called an outlier} 
#' }
#' @details If the sample was an outlier in principal component (PC) space, \code{reason}
#'   lists the PC(s) in which it was an outlier. See ome-specific details of outlier calling 
#'   in the supplementary methods of the manuscript. 
#' @source <https://docs.google.com/spreadsheets/d/1DetAMovcmMJqulEA41yBLd4I9Q1XSHPH1U8U2aVMShg/edit#gid=2057863805>
"OUTLIERS"


### Differential analysis results ####

#' @title Differential analysis of RNA-seq datasets 
#' @description Timewise summary statistics and training FDR from 
#'     differential analysis (DA) that tests the effect of training on the 
#'     expression of each gene within each sex. One data frame per tissue.   
#' @format A data frame with 22 variables:
#' \describe{
#'   \item{\code{feature}}{`r feature()`}
#'   \item{\code{assay}}{`r assay()`}
#'   \item{\code{assay_code}}{`r assay_code()`}
#'   \item{\code{tissue}}{`r tissue()`}
#'   \item{\code{tissue_code}}{`r tissue_code()`}
#'   \item{\code{feature_ID}}{`r feature_ID()`}
#'   \item{\code{sex}}{`r sex()`}
#'   \item{\code{comparison_group}}{`r comparison_group()`}
#'   \item{\code{p_value}}{`r p_value_da()`}
#'   \item{\code{adj_p_value}}{`r adj_p_value_da()`}
#'   \item{\code{logFC}}{`r logFC()`}
#'   \item{\code{logFC_se}}{`r logFC_se()`}
#'   \item{\code{shrunk_logFC}}{double, log fold-change with shrinkage applied}
#'   \item{\code{shrunk_logFC_se}}{double, standard error of the shrunken log fold-change}
#'   \item{\code{zscore}}{`r zscore()`}
#'   \item{\code{covariates}}{`r covariates()`}
#'   \item{\code{removed_samples}}{`r removed_samples()`}
#'   \item{\code{comparison_average_intensity}}{`r comparison_average_intensity()`}
#'   \item{\code{comparison_average_intensity_se}}{`r comparison_average_intensity_se()`}
#'   \item{\code{reference_average_intensity}}{`r reference_average_intensity()`}
#'   \item{\code{reference_average_intensity_se}}{`r reference_average_intensity_se()`}
#'   \item{\code{selection_fdr}}{`r selection_fdr()`} 
#' }
#' @details 
#'   Reproduce our analysis with [MotrpacRatTraining6mo::transcript_timewise_da()]
#'   and [MotrpacRatTraining6mo::transcript_training_da()]. 
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
#' @format A data frame 18 variables:
#' \describe{
#'   \item{\code{feature}}{`r feature()`}
#'   \item{\code{assay}}{`r assay()`}
#'   \item{\code{assay_code}}{`r assay_code()`}
#'   \item{\code{tissue}}{`r tissue()`}
#'   \item{\code{tissue_code}}{`r tissue_code()`}
#'   \item{\code{feature_ID}}{`r feature_ID()`}
#'   \item{\code{sex}}{`r sex()`}
#'   \item{\code{comparison_group}}{`r comparison_group()`}
#'   \item{\code{p_value}}{`r p_value_da()`}
#'   \item{\code{adj_p_value}}{`r adj_p_value_da()`}
#'   \item{\code{logFC}}{`r logFC()`}
#'   \item{\code{logFC_se}}{`r logFC_se()`}
#'   \item{\code{tscore}}{`r tscore()`}
#'   \item{\code{covariates}}{`r covariates()`}
#'   \item{\code{numNAs}}{`r numNAs()`}
#'   \item{\code{comparison_average_intensity}}{`r comparison_average_intensity()`}
#'   \item{\code{reference_average_intensity}}{`r reference_average_intensity()`}
#'   \item{\code{selection_fdr}}{`r selection_fdr()`} 
#' }
#' @details 
#'   Reproduce our analysis with [MotrpacRatTraining6mo::proteomics_timewise_da()]
#'   and [MotrpacRatTraining6mo::proteomics_training_da()]. 
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
#' @format A data frame with 27 variables:
#' \describe{
#'   \item{\code{feature}}{`r feature()`}
#'   \item{\code{assay}}{`r assay()`}
#'   \item{\code{assay_code}}{`r assay_code()`}
#'   \item{\code{tissue}}{`r tissue()`}
#'   \item{\code{tissue_code}}{`r tissue_code()`}
#'   \item{\code{feature_ID}}{`r feature_ID()`}
#'   \item{\code{dataset}}{`r dataset_metab()`}
#'   \item{\code{site}}{character, Chemical Analysis Site (CAS) name}
#'   \item{\code{is_targeted}}{logical, is this a targeted platform?}
#'   \item{\code{sex}}{`r sex()`}
#'   \item{\code{comparison_group}}{`r comparison_group()`}
#'   \item{\code{p_value}}{`r p_value_da()`}
#'   \item{\code{adj_p_value}}{`r adj_p_value_da()`}
#'   \item{\code{logFC}}{`r logFC()`}
#'   \item{\code{logFC_se}}{`r logFC_se()`}
#'   \item{\code{tscore}}{`r tscore()`}
#'   \item{\code{covariates}}{`r covariates()`}
#'   \item{\code{comparison_average_intensity}}{`r comparison_average_intensity()`}
#'   \item{\code{reference_average_intensity}}{`r reference_average_intensity()`}
#'   \item{\code{metabolite_refmet}}{character, RefMet name of metabolite}
#'   \item{\code{cv}}{double, feature coefficient of variation in the dataset}
#'   \item{\code{metabolite}}{character, name of metabolite as appears in the CAS's data}
#'   \item{\code{control_cv}}{double, feature coefficient of variation in the dataset}
#'   \item{\code{mz}}{double, mass over charge}
#'   \item{\code{rt}}{double, retention time}
#'   \item{\code{neutral_mass}}{double, neutral mass}
#'   \item{\code{selection_fdr}}{`r selection_fdr()`} 
#' }
#' @details 
#'   Reproduce our analysis with [MotrpacRatTraining6mo::metab_training_da()]
#'   and [MotrpacRatTraining6mo::metab_timewise_da()]. 
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


#' @title Meta-regression of metabolomics differential analysis results
#' @description Timewise summary statistics and training FDR from 
#'   differential analysis (DA) that tests the effect of training on each 
#'   metabolite within each sex. One data frame per tissue.
#'   
#'   For metabolites measured on more than one 
#'   platform, the consensus result may be provided. For metabolites measured 
#'   on a single platform, these results are identical to [METAB_DA]. 
#'      
#' @format A data frame with 30 variables:
#' \describe{
#'   \item{\code{feature}}{`r feature()`}
#'   \item{\code{assay}}{`r assay()`}
#'   \item{\code{assay_code}}{`r assay_code()`}
#'   \item{\code{tissue}}{`r tissue()`}
#'   \item{\code{tissue_code}}{`r tissue_code()`}
#'   \item{\code{feature_ID}}{`r feature_ID()`}
#'   \item{\code{dataset}}{`r dataset_metab()`}
#'   \item{\code{site}}{character, Chemical Analysis Site (CAS) name}
#'   \item{\code{is_targeted}}{logical, is this a targeted platform?}
#'   \item{\code{sex}}{`r sex()`}
#'   \item{\code{comparison_group}}{`r comparison_group()`}
#'   \item{\code{p_value}}{`r p_value_da()`}
#'   \item{\code{adj_p_value}}{`r adj_p_value_da()`}
#'   \item{\code{logFC}}{`r logFC()`}
#'   \item{\code{logFC_se}}{`r logFC_se()`}
#'   \item{\code{tscore}}{`r tscore()`}
#'   \item{\code{zscore}}{`r zscore()`}
#'   \item{\code{covariates}}{`r covariates()`}
#'   \item{\code{comparison_average_intensity}}{`r comparison_average_intensity()`}
#'   \item{\code{reference_average_intensity}}{`r reference_average_intensity()`}
#'   \item{\code{metabolite_refmet}}{character, RefMet name of metabolite}
#'   \item{\code{cv}}{double, feature coefficient of variation in the dataset}
#'   \item{\code{metabolite}}{character, name of metabolite as appears in the CAS's data}
#'   \item{\code{control_cv}}{double, feature coefficient of variation in the dataset}
#'   \item{\code{mz}}{double, mass over charge}
#'   \item{\code{rt}}{double, retention time}
#'   \item{\code{neutral_mass}}{double, neutral mass}
#'   \item{\code{meta_reg_het_p}}{`r meta_reg_het_p()`}
#'   \item{\code{meta_reg_pvalue}}{`r meta_reg_pvalue()`}
#'   \item{\code{selection_fdr}}{`r selection_fdr()`} 
#' }
#' @details TODO
#' 
#'   Reproduce our analysis with TODO. 
#' 
#' @name METAB_DA_METAREG
"METAB_ADRNL_DA_METAREG"

#' @rdname METAB_DA_METAREG
"METAB_BAT_DA_METAREG"

#' @rdname METAB_DA_METAREG
"METAB_COLON_DA_METAREG"

#' @rdname METAB_DA_METAREG
"METAB_CORTEX_DA_METAREG"

#' @rdname METAB_DA_METAREG
"METAB_HEART_DA_METAREG"

#' @rdname METAB_DA_METAREG
"METAB_HIPPOC_DA_METAREG"

#' @rdname METAB_DA_METAREG
"METAB_HYPOTH_DA_METAREG"

#' @rdname METAB_DA_METAREG
"METAB_KIDNEY_DA_METAREG"

#' @rdname METAB_DA_METAREG
"METAB_LIVER_DA_METAREG"

#' @rdname METAB_DA_METAREG
"METAB_LUNG_DA_METAREG"

#' @rdname METAB_DA_METAREG
"METAB_OVARY_DA_METAREG"

#' @rdname METAB_DA_METAREG
"METAB_PLASMA_DA_METAREG"

#' @rdname METAB_DA_METAREG
"METAB_SKMGN_DA_METAREG"

#' @rdname METAB_DA_METAREG
"METAB_SKMVL_DA_METAREG"

#' @rdname METAB_DA_METAREG
"METAB_SMLINT_DA_METAREG"

#' @rdname METAB_DA_METAREG
"METAB_SPLEEN_DA_METAREG"

#' @rdname METAB_DA_METAREG
"METAB_TESTES_DA_METAREG"

#' @rdname METAB_DA_METAREG
"METAB_VENACV_DA_METAREG"

#' @rdname METAB_DA_METAREG
"METAB_WATSC_DA_METAREG"


#' @title Differential analysis of multiplexed immunoassays 
#' @description Timewise summary statistics and training FDR from 
#'     differential analysis (DA) that tests the effect of training on each 
#'     immunoassay analyte within each sex. One data frame per tissue. 
#' @format A data frame with 17 variables:
#' \describe{
#'   \item{\code{feature}}{`r feature()`}
#'   \item{\code{assay}}{`r assay()`}
#'   \item{\code{assay_code}}{`r assay_code()`}
#'   \item{\code{tissue}}{`r tissue()`}
#'   \item{\code{tissue_code}}{`r tissue_code()`}
#'   \item{\code{feature_ID}}{`r feature_ID()`}
#'   \item{\code{dataset}}{character, name of LUMINEX panel, e.g., 'rat-mag27plex'}
#'   \item{\code{sex}}{`r sex()`}
#'   \item{\code{comparison_group}}{`r comparison_group()`}
#'   \item{\code{p_value}}{`r p_value_da()`}
#'   \item{\code{adj_p_value}}{`r adj_p_value_da()`}
#'   \item{\code{logFC}}{`r logFC()`}
#'   \item{\code{logFC_se}}{`r logFC_se()`}
#'   \item{\code{covariates}}{`r covariates()`}
#'   \item{\code{comparison_average_intensity}}{`r comparison_average_intensity()`}
#'   \item{\code{reference_average_intensity}}{`r reference_average_intensity()`}
#'   \item{\code{selection_fdr}}{`r selection_fdr()`} 
#' }
#' @details 
#'   Reproduce our analysis with [MotrpacRatTraining6mo::immuno_timewise_da()]
#'   and [MotrpacRatTraining6mo::immuno_training_da()]. 
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


#' @title Differential analysis of ATAC-seq data
#' @description Timewise summary statistics and training FDR from 
#'     differential analysis (DA) that tests the effect of training on each 
#'     chromatin accessibility peak within each sex. One data frame per tissue. 
#' @format A data frame with 21 variables:
#' \describe{
#'   \item{\code{feature}}{`r feature()`}
#'   \item{\code{assay}}{`r assay()`}
#'   \item{\code{assay_code}}{`r assay_code()`}
#'   \item{\code{tissue}}{`r tissue()`}
#'   \item{\code{tissue_code}}{`r tissue_code()`}
#'   \item{\code{feature_ID}}{`r feature_ID()`}
#'   \item{\code{dataset}}{character, tissue code and Chemical Analysis Site that
#'     generated the data. Note only ATAC data generated by Stanford are included
#'     in this repository.}
#'   \item{\code{sex}}{`r sex()`}
#'   \item{\code{comparison_group}}{`r comparison_group()`}
#'   \item{\code{p_value}}{`r p_value_da()`}
#'   \item{\code{adj_p_value}}{`r adj_p_value_da()`}
#'   \item{\code{logFC}}{`r logFC()`}
#'   \item{\code{logFC_se}}{`r logFC_se()`}
#'   \item{\code{tscore}}{`r tscore()`}
#'   \item{\code{covariates}}{`r covariates()`}
#'   \item{\code{removed_samples}}{`r removed_samples()`}
#'   \item{\code{comparison_average_intensity}}{`r comparison_average_intensity()`}
#'   \item{\code{comparison_average_intensity_se}}{`r comparison_average_intensity_se()`}
#'   \item{\code{reference_average_intensity}}{`r reference_average_intensity()`}
#'   \item{\code{reference_average_intensity_se}}{`r reference_average_intensity_se()`}
#'   \item{\code{selection_fdr}}{`r selection_fdr()`} 
#' }
#' @details 
#'   Unfiltered ATAC differential analysis results are only available via download from Google Cloud Storage. 
#'   For example, <https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/ATAC_BAT_DA.rda> 
#'   is the file for brown adipose tissue (BAT) results. You can change the name of the file to specify other tissues including:
#'   HEART, HIPPOC, KIDNEY, LIVER, LUNG, SKMGN (gastrocnemius skeletal muscle), and WATSC (subcutaneous white adipose tissue).
#'   You can also use [MotrpacRatTraining6mo::get_rdata_from_url()] 
#'   to download differential analysis results for ATAC and METHYL.
#'   For more details about these files see the readme of this repository at 
#'   <https://github.com/MoTrPAC/MotrpacRatTraining6mo/blob/main/README.md>. 
#'   
#'   For training-regulated ATAC features at 5% FDR, see [TRAINING_REGULATED_FEATURES].
#' 
#'   Reproduce our analysis with [MotrpacRatTraining6mo::atac_timewise_da()]
#'   and [MotrpacRatTraining6mo::atac_training_da()]. 
#' 
#' @name ATAC_DA
NULL


#' @title Differential analysis results of RRBS data
#' @description Timewise summary statistics and training FDR from 
#'     differential analysis (DA) that tests the effect of training on each 
#'     methylation feature within each sex. One data frame per tissue. 
#' @format A data frame with 20 variables:
#' \describe{
#'   \item{\code{feature}}{`r feature()`}
#'   \item{\code{assay}}{`r assay()`}
#'   \item{\code{assay_code}}{`r assay_code()`}
#'   \item{\code{tissue}}{`r tissue()`}
#'   \item{\code{tissue_code}}{`r tissue_code()`}
#'   \item{\code{feature_ID}}{`r feature_ID()`}
#'   \item{\code{sex}}{`r sex()`}
#'   \item{\code{comparison_group}}{`r comparison_group()`}
#'   \item{\code{p_value}}{`r p_value_da()`}
#'   \item{\code{adj_p_value}}{`r adj_p_value_da()`}
#'   \item{\code{logFC}}{`r logFC()`}
#'   \item{\code{zscore}}{`r zscore()`}
#'   \item{\code{covariates}}{`r covariates()`}
#'   \item{\code{removed_samples}}{`r removed_samples()`}
#'   \item{\code{Chr}}{integer, chromosome}
#'   \item{\code{Locus}}{character, base pair range of feature}
#'   \item{\code{EntrezID}}{character, Entrez ID of closest gene}
#'   \item{\code{Symbol}}{character, gene symbol of closest gene}
#'   \item{\code{fscore}}{the LRT fscore}
#'   \item{\code{selection_fdr}}{`r selection_fdr()`}
#' }
#' @details 
#'   Unfiltered METHYL differential analysis results are only available via download from Google Cloud Storage. 
#'   For example, <https://storage.googleapis.com/motrpac-rat-training-6mo-extdata/METHYL_BAT_DA.rda> 
#'   is the file for brown adipose tissue (BAT) results. You can change the name of the file to specify other tissues including:
#'   HEART, HIPPOC, KIDNEY, LIVER, LUNG, SKMGN (gastrocnemius skeletal muscle), and WATSC (subcutaneous white adipose tissue).
#'   You can also use [MotrpacRatTraining6mo::get_rdata_from_url()] 
#'   to download differential analysis results for ATAC and METHYL.
#'   For more details about these files see the readme of this repository at 
#'   <https://github.com/MoTrPAC/MotrpacRatTraining6mo/blob/main/README.md>. 
#'   
#'   The differential analysis results in this data frame were computed using [edgeR::glmQLFTest()] after
#'   removing a few outlier low quality samples from the kidney dataset.
#'   
#'   For training-regulated METHYL features at 5% FDR, see [TRAINING_REGULATED_FEATURES].
#'
#'   Reproduce our analysis with [MotrpacRatTraining6mo::rrbs_differential_analysis()].
#'
#' @name METHYL_DA
NULL


#' @title Training-regulated features
#' @description Differential analysis results for training-regulated features at 5% FDR
#' @format A data frame with 279002 rows and 18 variables:
#' \describe{
#'   \item{\code{feature}}{`r feature()`}
#'   \item{\code{assay}}{`r assay()`}
#'   \item{\code{assay_code}}{`r assay_code()`}
#'   \item{\code{tissue}}{`r tissue()`}
#'   \item{\code{tissue_code}}{`r tissue_code()`}
#'   \item{\code{feature_ID}}{`r feature_ID()`}
#'   \item{\code{non_redundant_feature_ID}}{character, non-redundant feature identifier defined for \code{feature_ID}s 
#'       with measurements from multiple platforms. \code{feature} uses \code{non_redundant_feature_ID} 
#'       when available and \code{feature_ID} otherwise. See [REPEATED_FEATURES] for more details.}
#'   \item{\code{platform}}{`r dataset()`}
#'   \item{\code{sex}}{`r sex()`}
#'   \item{\code{training_group}}{character, training time point in weeks, one of "1w", "2w", "4w", "8w". 
#'       Corresponds to \code{comparison_group} in the \code{*_DA} tables.}
#'   \item{\code{timewise_logFC}}{double, log fold-change of the training group specified by 
#'       \code{sex} and \code{training_group} (e.g., 1-week females) relative to 
#'       the sex-matched sedentary controls. Corresponds to \code{logFC} in the \code{*_DA} tables.}
#'   \item{\code{timewise_logFC_se}}{double, standard error of \code{timewise_logFC}.
#'       Corresponds to \code{logFC_se} in the \code{*_DA} tables.}
#'   \item{\code{timewise_p_value}}{double, timewise p-value corresponding to the 
#'       contrast between the training group (e.g., 1-week females) and the sex-matched 
#'       sedentary controls. Corresponds to \code{p_value} in the \code{*_DA} tables.}
#'   \item{\code{timewise_zscore}}{double, smoothed z-score corresponding to 
#'       \code{timewise_p_value} used for graphical clustering.
#'       Corresponds to \code{zscore} or \code{tscore} in the \code{*_DA} tables.}
#'   \item{\code{meta_reg_het_p}}{`r meta_reg_het_p()`}
#'   \item{\code{meta_reg_pvalue}}{`r meta_reg_pvalue()`}
#'   \item{\code{training_p_value}}{double, combined training p-value determined from 
#'       F-tests or LRTs comparing a full model with ome-specific technical covariates 
#'       and training group as a factor variable against a reduced model with only 
#'       technical covariates. Sex-specific training p-values were merged using the 
#'       sum of logs. One value per features (not per training group)}
#'   \item{\code{training_q}}{double, IHW FDR-adjusted \code{training_p_value}; 
#'       training-regulated features were defined by \code{training_q < 0.05}.
#'       Corresponds to \code{selection_fdr} in the \code{*_DA} tables.}
#' }
"TRAINING_REGULATED_FEATURES"


## Graphical analysis ####

#' @title \code{repfdr} state assignments 
#' @description \code{repdfr} state assignments (up, down, or null) for each 
#'   training-regulated feature (5% FDR) for each sex at each time point, 
#'   which specify node assignments for each differential feature. 
#'   Missing values indicate that the \code{repfdr} posterior probabilities 
#'   did not meet the cutoff for that feature.
#' @format A data frame with 34244 rows and 10 variables:
#' \describe{
#'   \item{\code{feature}}{`r feature()`}
#'   \item{\code{ome}}{`r assay()`}
#'   \item{\code{tissue}}{`r tissue()`. Note that VENACV, OVARY, TESTES, 
#'     were not included in the graphical representation of differential 
#'     features due to missing groups (e.g., females trained for 1 week).}
#'   \item{\code{feature_ID}}{`r feature_ID()`}
#'   \item{\code{state_1w}}{character, state (1, up-regulated; 0, null; -1, down-regulated) 
#'     of the feature in each sex (F, females; M, males) at the 1-week training time point, 
#'     relative to sex-matched untrained animals}
#'   \item{\code{state_2w}}{character, state of the feature in each sex at the 2-week training time point}
#'   \item{\code{state_4w}}{character, state of the feature in each sex at the 4-week training time point}
#'   \item{\code{state_8w}}{character, state of the feature in each sex at the 8-week training time point}
#'   \item{\code{path}}{character, assigned states from weeks 1-8, separated by "->". 
#'     This represents a feature's full path through the graph. 
#'     NA if the state at any of the four time points is NA.}
#'   \item{\code{tissue_path}}{character, assigned states from weeks 1-8, separated by "->". 
#'     This represents a feature's full path through the graph. 
#'     NA if the state at any of the four time points is NA.} 
#' }
#' @details 
#'   Given the posteriors Pr(h|z_i) computed using [repfdr::repfdr()] where h is a configuration
#'   vector in {-1,0,1}^8 (specifying the 8 analyzed groups, 4 time points in males and females),
#'   and z_i is the vector of z-scores of analyte i, we assign analytes to "states". 
#'   A state is a tuple (s_{m,j}, s_{f,j}), where s_{m,j} is the differential abundance 
#'   state null, up, or down (0,1, and -1 in the notation above, respectively) in males 
#'   at time point j, and s_{f,j} is defined similarly for females (at time point j). 
#'   Thus, we have nine possible states in each time point. 
#'   For example, assume we inspect analyte i in time point j, asking if the abundance is 
#'   up-regulated in males while null in females. Then, we sum over all posteriors Pr(h|z_i) 
#'   such that males are up-regulated and females have 0. 
#'   If the resulting value is greater than 0.5, then we say that analyte i belongs to the 
#'   node set S(s_{m,j}, s_{f,j}). Thus, we use S(s_{m,j}, s_{f,j}) to denote all analytes 
#'   that belong to a state (s_{m,j}, s_{f,j}). 
#'   Then, for every pair of states from adjacent time points j and j+1 we define their 
#'   edge set as the intersection of S(s_{m,j}, s_{f,j}) and S(s_{m,j+1}, s_{f,j+1}). 
#'   Thus, thenode sets edge sets together define a tree structure that represent different 
#'   differential patterns over sex and time.
"GRAPH_STATES"


#' @title KEGG and Reactome parent pathways 
#' @description List where names are KEGG or Reactome terms IDs and values are 
#'   the parent pathway(s)
#' @format List of length 21697
#' @details 
#'   Pathway hierarchies for KEGG were retrieved with [KEGGREST::keggGet()], e.g. \code{keggGet(kegg_query)[[1]]$CLASS}. 
#'   Reactome pathway parents were extracted from <https://reactome.org/download/current/ReactomePathwaysRelation.txt>. 
"PATHWAY_PARENTS"


#' @title Graph pathway enrichment results 
#' @description Pathway enrichment results for graphical clusters (nodes, edges, and paths) of interest
#' @format A data frame with 156906 rows and 22 variables:
#' \describe{
#'   \item{\code{query}}{character, not used, carried over from [gprofiler2::gost()] output}
#'   \item{\code{significant}}{logical, not used, carried over from [gprofiler2::gost()] output}
#'   \item{\code{term_size}}{double, effective pathway size from [gprofiler2::gost()] output}
#'   \item{\code{query_size}}{integer, size of input, i.e. list of Ensembl genes associated with differential features}
#'   \item{\code{intersection_size}}{double, size of the intersection between the input and the pathway members}
#'   \item{\code{precision}}{double, the proportion of genes in the input list that are annotated to the function (defined as \code{intersection_size/query_size})}
#'   \item{\code{recall}}{double, the proportion of functionally annotated genes that the query recovers (defined as \code{intersection_size/term_size})}
#'   \item{\code{term_id}}{character, pathway term ID}
#'   \item{\code{source}}{character, database corresponding to the pathway, one of: "KEGG", "REAC"}
#'   \item{\code{term_name}}{character, pathway name}
#'   \item{\code{effective_domain_size}}{integer, size of the custom background Ensembl gene set}
#'   \item{\code{source_order}}{integer, not used, carried over from [gprofiler2::gost()] output}
#'   \item{\code{parents}}{list, pathway parent(s)}
#'   \item{\code{evidence_codes}}{character, not used, carried over from [gprofiler2::gost()] output}
#'   \item{\code{intersection}}{character, intersection between input and pathway (Ensembl IDs). NA for metabolomics enrichments}
#'   \item{\code{gost_adj_p_value}}{double, BH-adjusted p-value returned by [gprofiler2::gost()], 
#'     ignored because p-values are only adjusted within each tissue/ome/cluster combination. Use the \code{adj_p_value} column instead.}
#'   \item{\code{computed_p_value}}{double, nominal hypergeometric p-value, computed from the [gprofiler2::gost()] output}
#'   \item{\code{cluster}}{character, graphical cluster (node, edge, or path) name}
#'   \item{\code{tissue}}{`r tissue()`. Note that VENACV, OVARY, TESTES, 
#'     were not included in the graphical representation of differential 
#'     features due to missing groups (e.g., females trained for 1 week).}
#'   \item{\code{ome}}{`r assay()`}
#'   \item{\code{kegg_id}}{character, pathway ID returned from [FELLA::enrich()]}
#'   \item{\code{adj_p_value}}{double, IHW FDR, calculated using [IHW::ihw()] with \code{tissue} as a covariate} 
#' }
#' @details All non-metabolite training-regulated features (5% FDR) were mapped 
#'   to Ensembl gene symbols using [FEATURE_TO_GENE]. Training-regulated metabolites 
#'   were mapped to KEGG IDs. For each graphical cluster of interest (i.e., 
#'   the ten largest paths, two largest nodes, and two largest single edges with 
#'   at least 20 features in each tissue, as well as all 8-week nodes), 
#'   we performed pathway enrichment analysis separately for the Ensembl genes 
#'   (or KEGG IDs for metabolites) associated with differential features in each ome. 
#'   
#'   For gene-centric omes (i.e., all but metabolomics) 
#'   we performed enrichment analysis of KEGG and REACTOME rat pathways (organism "rnorvegicus") 
#'   using [gprofiler2::gost()] with custom backgrounds defined by [GENE_UNIVERSES].  
#'   Only pathways with at least 10 and up to 200 members were tested. Because [gprofiler2::gost()]
#'   only returns adjusted p-values, we recalculated nominal p-values using a one-tailed hypergeometric test, 
#'   which is consistent with how [gprofiler2::gost()] calculates enrichments. 
#'   See [MotrpacRatTraining6mo::cluster_pathway_enrichment()] for implementation. 
#'   
#'   For metabolites, 
#'   we performed enrichment of KEGG pathways using the hypergeometric method in [FELLA::enrich()] 
#'   with custom backgrounds defined by [GENE_UNIVERSES]. See [MotrpacRatTraining6mo::run_fella()] for implementation. 
#'   
#'   Pathway enrichment analysis p-values were adjusted across all results using Independent 
#'   Hypothesis Weighting (IHW) with tissue as a covariate.
#'    
"GRAPH_PW_ENRICH"


#' @title Gene-centric universes
#' @description Nested vectors of the set of genes tested for each dataset used to specify custom backgrounds for enrichment analyses 
#' @format A list of lists of lists of vectors. Access a vector of measured/tested genes with \code{GENE_UNIVERSES[[gene_identifier]][[ome]][[tissue]]}, where:
#' \describe{
#'   \item{\code{gene_identifier}}{one of "entrez_gene", "gene_symbol", "ensembl_gene", "rgd_gene"}
#'   \item{\code{ome}}{`r assay()`}
#'   \item{\code{tissue}}{`r tissue()`}
#' }
#' @details If you are performing any kind of enrichment test that does not accept the 
#'   full set of measured features as the input, you should specify a custom background. 
#'   For example, if you are performing enrichment for an unordered list of differential 
#'   features, then the background should specify the full set of features tested during 
#'   differential analysis. This object provides the gene-centric lists of "expressed" 
#'   features, i.e. those tested during differential analysis, for all assay and 
#'   tissue combinations. 
#' 
#'   For convenience, universes are specified with four different types of gene identifers: 
#'   Entrez/NCBI IDs ("entrez_gene"), gene symbols ("gene_symbol"), Ensembl IDs ("ensembl_gene"), 
#'   and RGD IDs ("rgd_gene"). This object is a list of lists, with one list per type of gene identifier. 
#'   List structure is \emph{gene identifier > assay abbreviation > tissue abbreviation}, 
#'   e.g. \code{GENE_UNIVERSES$entrez_gene$PHOSPHO$HEART} is the unique list of Entrez genes associated 
#'   with all measured phosphosites in the heart.
#' 
#'   Universes for ATAC and METHYL correspond to the universe of expressed genes, 
#'   i.e. TRNSCRPT, not actually the full set of genes associated with any features 
#'   measured in these epigenetic assays.
#' 
#'   Regardless of the gene identifier, features in METAB universe lists are always KEGG IDs, 
#'   not actually gene identifiers.
#'   
"GENE_UNIVERSES"


#' @title \code{repfdr} inputs
#' @description Feature by state data frames that define the input for [repfdr::repfdr()] 
#' @format List of data frames:
#' \describe{
#'   \item{\code{zs_info}}{feature-level metadata for training-regulated features included in graphical analysis}
#'   \item{\code{zs_smoothed}}{feature by state data frame where values are smoothed timewise z-scores} 
#'   \item{\code{nominal_ps}}{feature by state data frame where values are timewise nominal p-values}
#'   \item{\code{effects}}{feature by state data frame where values are timewise effects or log fold-changes}
#' } 
#'
"REPFDR_INPUTS"


#' @title Graph components 
#' @description The feature sets selected using the \code{repfdr} results. Each 
#'   value is a named list, where the name is the label for the edge or node, 
#'   and the members are all of the features that belong to that label, 
#'   in the format \code{\link{ASSAY_ABBREV};\link{TISSUE_ABBREV};[feature_ID]}. 
#' @format List of lists: 
#' \describe{
#'   \item{\code{edge_sets}}{named list, where the name is the label for the edge in the 
#'     format \code{[fromNode]---[toNode]} and the value is a vector of features that belong to that edge}
#'   \item{\code{node_sets}}{named list, where the name is the label for the node in the 
#'     format \code{[week]w_F[state]_M[state]} and the value is a vector of features that belong to that edge} 
#' }
#' @details Nodes are named in the format \code{[week]w_F[state]_M[state]}, where 
#'   \code{[week]} is one of 1, 2, 4, or 8 to specify the training point, \code{[state]} 
#'   takes a value of 0 (null), 1 (up-regulated), or -1 (down-regulated). Altogether, 
#'   then name of a node specifies the state of a feature in each sex at a specific 
#'   time point. For example, the node \code{1w_F-1_M-1} specifies all features that 
#'   are down-regulated in females (\code{F-1}) and males (\code{M-1}) at the 1-week
#'   training time point (\code{1w}). 
#'   
"GRAPH_COMPONENTS"


#' @title \code{repfdr} results  
#' @description Raw [repfdr::repfdr()] results from which the graphical state assignments were determined
#' @format List:
#' \describe{
#'   \item{\code{repfdr_em_res}}{a list with \code{repfdr}'s EM results}
#'   \item{\code{repfdr_clusters}}{\code{repfdr}'s configurations} 
#'   \item{\code{repfdr_clusters_str}}{\code{repfdr}'s configurations, string representation}
#'   \item{\code{repfdr_clusters_pi}}{configuration's inferred priors} 
#' }
#' @details [repfdr::repfdr()] is an algorithm suggested by Yekutieli and Heller in 2014 
#'   (Bioinformatics) for analysis of p-values or z-scores from different resources. 
#'   It is based on the assumption that z-scores from each resource (a time point 
#'   from a specific sex in our case) is either positive, null, or negative. 
#'   Then the algorithm internally learns the mixture distribution in each resource 
#'   and the dependencies among them. In the process, a simplifying assumption is 
#'   made about conditional independence between the z-scores given their state 
#'   (so for example, two high scores from weeks 4 and 8 in males are independent 
#'   given that we know that they are positive, non-null cases). 
#'   
#'   Our analysis pipeline is as follows: we run [repfdr::repfdr()] on the z-score data matrix 
#'   of all training-regulated features at 5% IWH FDR (see [REPFDR_INPUTS]). Then we use the output to extract 
#'   the analytes of each state in each time point. Here a state means one of 
#'   (male up, male null, male down) x (female up, female null, female down) 
#'   for each time point. Once these are selected we call them the \code{node_sets}. 
#'   Then, for each pair of node \code{(x,y)} such that \code{y} is from a time point that is adjacent 
#'   and after \code{x} (e.g., \code{x} is a node from week 4 and \code{y} is a node from week 8), 
#'   we define their edge set as the intersection of their analytes.
#'   
#'   Reproduce our analysis with [MotrpacRatTraining6mo::bayesian_graphical_clustering()]
#'   and [MotrpacRatTraining6mo::repfdr_wrapper()]. 
#'   
"REPFDR_RES"
