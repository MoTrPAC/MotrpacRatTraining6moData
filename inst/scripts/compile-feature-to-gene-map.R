#!/bin/R

# FOR REFERENCE ONLY 
# This script was **not** re-run to generate `FEATURE_TO_GENE`.
# Instead, `master_feature_to_gene` from gs://motrpac-data-freeze-pass/pass1b-06/v1.1/analysis/resources/master_feature_to_gene_20211116.RData
# was used to make `FEATURE_TO_GENE`, which was output by the original version of this script: 
# <https://github.com/MoTrPAC/motrpac-mawg/blob/master/pass1b-06/tools/feature-to-gene/compile-feature-to-gene-map.R> 

library(MotrpacBicQC)
# library(gprofiler2)
# library(ChIPseeker)
# library(GenomicFeatures)
library(data.table)
library(devtools)
load_all()

feature_to_gene_map_list = list() # all features included in clustering input that map to a gene or KEGG ID 

## functions ########################################################################################################

#' Read a single file from Google Cloud into a data.table
#'
#' @param path GCP path, i.e. starts with "gs://"
#' @param sep column separator to use with "fread"
#' @param tmpdir scratch path to download files from GCP
#' @param GSUTIL_PATH path to "gsutil" on your computer
#' @param check_first check if file exists before downloading it. read in existing file if it exists. should be set to TRUE if you are running this function in parallel
#'
#' @return A data.table
dl_read_gcp = function(path,sep='\t',tmpdir='/tmp',GSUTIL_PATH='gsutil',check_first=F){
  system(sprintf('mkdir -p %s',tmpdir))
  # download
  new_path = sprintf('%s/%s',tmpdir,basename(path))
  # only download if it doesn't exist to avoid conflicts when running this script in parallel
  if(check_first){
    if(!file.exists(new_path)){
      cmd = sprintf('%s cp %s %s', GSUTIL_PATH, path, tmpdir)
      system(cmd,ignore.stdout = T,ignore.stderr = T)
    }
  }else{
    cmd = sprintf('%s cp %s %s', GSUTIL_PATH, path, tmpdir)
    system(cmd,ignore.stdout = T,ignore.stderr = T)
  }
  # read in the data as a data.table
  if(file.exists(new_path)){
    dt = fread(new_path,sep=sep,header=T)
    return(dt)
  }
  warning(sprintf("gsutil file %s does not exist.\n",path))
  return()
}

# usage: load_file_pattern('gs://mawg-data/pass1b-06/transcript-rna-seq/dea/', 'timewise-dea')
load_file_pattern = function(bucket, pattern, verbose=T, .rbind=T, add_tissue=T, tmpdir="/tmp"){
  data_list = list()
  # first, download all files
  orig_file_list = system(sprintf('gsutil ls %s | grep %s', bucket, pattern), intern=T)
  if(!all(file.exists(paste0(tmpdir, '/', basename(orig_file_list))))){
    system(sprintf("gsutil -m cp %s %s", paste0(orig_file_list, collapse=' '), tmpdir))
  }
  
  file_list = basename(orig_file_list)
  for (f in file_list){
    dt = fread(sprintf("%s/%s", tmpdir, f))
    if(add_tissue){
      if(!"tissue" %in% colnames(dt)){
        splits = unname(unlist(strsplit(f, "_"))) 
        curr_tissue = splits[grepl("^t[0-9]", splits)]
        dt[,tissue := curr_tissue]
      }
      data_list[[f]] = dt
    }
    if(verbose){print(f)}
  }
  
  if(!.rbind){
    return(data_list)
  }
  
  data = rbindlist(data_list, fill=T)
  return(data)
}

############################################################################################################################## ~
# feature_to_gene_map ####
############################################################################################################################## ~

#### transcript-rna-seq ########################################################################################################

transcript_rna_seq = dl_read_gcp("gs://motrpac-data-freeze-pass/pass1b-06/v1.0/analysis/transcriptomics/transcript-rna-seq/normalized-data/motrpac_pass1b-06_transcript-rna-seq_feature-mapping_20210721.txt", sep='\t')
# only keep ensembl
feature_to_gene_map_list[['transcript-rna-seq']] = unique(transcript_rna_seq[,.(assay, feature_ID, ensembl_gene)])

#### epigen-rrbs ########################################################################################################

# epigen_rrbs = load_file_pattern("gs://motrpac-data-freeze-pass/pass1b-06/v1.0/analysis/epigenomics/epigen-rrbs/normalized-data/", "normalized-data-feature-annot.txt")
# epigen_rrbs = unique(epigen_rrbs[,.(featureID, Chr, LocStart, LocEnd)])
# epigen_rrbs[,Chr := as.character(Chr)]
# epigen_rrbs[grepl("^chrX", featureID), Chr := 'X']
# epigen_rrbs[grepl("^chrY", featureID), Chr := 'Y']
# 
# setnames(epigen_rrbs, c("Chr", "LocStart", "LocEnd", "featureID"), c("chrom","start","end","feature_ID"))
# pa = get_peak_annotations(epigen_rrbs)
# pa[,assay := "epigen-rrbs"]
# 
# # check annotations
# for(anno in unique(pa[,custom_annotation])){
#   print(anno)
#   hist(pa[custom_annotation==anno & geneStrand ==1, relationship_to_gene], xlim=c(-8000,8000), breaks=10000, main=sprintf("%s, positive strand", anno))
#   hist(pa[custom_annotation==anno & geneStrand ==2, relationship_to_gene], xlim=c(-8000,8000), breaks=10000, main=sprintf("%s, negative strand", anno))
# }
# 
# pa = pa[,.(assay, feature_ID, chrom, start, end, width,
#            chipseeker_annotation, custom_annotation, distanceToTSS, relationship_to_gene, 
#            ensembl_gene, geneStart, geneEnd, geneLength, geneStrand)]
# 
# date='20211110'
# outfile = sprintf('pass1b-06_epigen-rrbs_feature-mapping_%s.txt',date)
# write.table(pa, outfile,
#             col.names=T, row.names=F, sep='\t', quote=F)
# system(sprintf("gsutil cp %s gs://mawg-data/pass1b-06/epigen-rrbs/mapping/", outfile))

epigen_rrbs = dl_read_gcp("gs://mawg-data/pass1b-06/epigen-rrbs/mapping/pass1b-06_epigen-rrbs_feature-mapping_20211110.txt", sep='\t')
feature_to_gene_map_list[['epigen-rrbs']] = unique(epigen_rrbs[,.(assay, feature_ID, ensembl_gene, custom_annotation, relationship_to_gene)])

#### epigen-atac-seq ########################################################################################################

epigen_atac_seq = dl_read_gcp("gs://mawg-data/pass1b-06/epigen-atac-seq/mapping/pass1b-06_epigen-atac-seq_feature-mapping_20211110.txt", sep='\t')
feature_to_gene_map_list[['epigen-atac-seq']] = epigen_atac_seq[,.(assay, feature_ID, ensembl_gene, custom_annotation, relationship_to_gene)]

#### prot-pr ########################################################################################################

prot_pr = load_file_pattern("gs://motrpac-data-freeze-pass/pass1b-06/v1.0/analysis/proteomics-untargeted/prot-pr/normalized-data/", "feature-annot.txt")
prot_pr[,assay := 'prot-pr']
prot_pr = prot_pr[is_contaminant == FALSE]
setnames(prot_pr, "entrez_id", "entrez_gene")
prot_pr = unique(prot_pr[,.(assay, feature_ID, gene_symbol, entrez_gene, protein_id)])

feature_to_gene_map_list[['prot-pr']] = prot_pr

#### prot-ph ########################################################################################################

prot_ph = load_file_pattern("gs://motrpac-data-freeze-pass/pass1b-06/v1.0/analysis/proteomics-untargeted/prot-ph/normalized-data/", "feature-annot.txt")
prot_ph[,assay := 'prot-ph']
prot_ph = prot_ph[is_contaminant == FALSE]
setnames(prot_ph, "entrez_id", "entrez_gene")
prot_ph = unique(prot_ph[,.(assay, feature_ID, gene_symbol, entrez_gene, protein_id)])

feature_to_gene_map_list[['prot-ph']] = prot_ph

#### prot-ac ########################################################################################################

prot_ac = load_file_pattern("gs://motrpac-data-freeze-pass/pass1b-06/v1.0/analysis/proteomics-untargeted/prot-ac/normalized-data/", "feature-annot.txt")
prot_ac[,assay := 'prot-ac']
prot_ac = prot_ac[is_contaminant == FALSE]
setnames(prot_ac, "entrez_id", "entrez_gene")
prot_ac = unique(prot_ac[,.(assay, feature_ID, gene_symbol, entrez_gene, protein_id)])

feature_to_gene_map_list[['prot-ac']] = prot_ac

#### prot-ub ########################################################################################################

prot_ub = load_file_pattern("gs://motrpac-data-freeze-pass/pass1b-06/v1.0/analysis/proteomics-untargeted/prot-ub/normalized-data/", "feature-annot.txt")
prot_ub[,assay := 'prot-ub']
prot_ub = prot_ub[is_contaminant == FALSE]
setnames(prot_ub, "entrez_id", "entrez_gene")
prot_ub = unique(prot_ub[,.(assay, feature_ID, gene_symbol, entrez_gene, protein_id)])

feature_to_gene_map_list[['prot-ub']] = prot_ub

#### metab ########################################################################################################

metab = data.table(MotrpacBicQC::get_and_validate_mdd())
# define it both ways 
metab1 = unique(metab[,.(metabolite_name, kegg_id)])
setnames(metab1, "metabolite_name", "feature_ID")
metab2 = unique(metab[,.(refmet_name, kegg_id)])
setnames(metab2, "refmet_name", "feature_ID")
metab = unique(rbindlist(list(metab1, metab2)))

metab[,assay := 'metab']
metab = metab[!is.na(kegg_id)]
metab[,gene_symbol := kegg_id] # a bit of a misnomer, but this will work 

feature_to_gene_map_list[['metab']] = metab

#### immunoassay ########################################################################################################

immunoassay = dl_read_gcp("gs://mawg-data/pass1b-06/immunoassay/mapping/pass1b-06_immunoassay_feature-mapping_20211109.txt", sep='\t')
feature_to_gene_map_list[['immunoassay']] = immunoassay

### save feature_to_gene_map ########################################################################################################

#save(feature_to_gene_map_list, file="feature_to_gene_map_list_20211116.RData")
#load("feature_to_gene_map_list_20211116.RData")

############################################################################################################################## ~
# add entrez, ensembl, symbol ####
############################################################################################################################## ~

feature_to_gene_map = rbindlist(feature_to_gene_map_list, fill=T)
# make sure all rows have assay, feature_ID
stopifnot(nrow(feature_to_gene_map[is.na(feature_ID)])==0)
stopifnot(nrow(feature_to_gene_map[is.na(assay)])==0)

# table(is.na(feature_to_gene_map[,gene_symbol]))
# table(is.na(feature_to_gene_map[,entrez_gene]))
# table(is.na(feature_to_gene_map[,ensembl_gene]))
# table(is.na(feature_to_gene_map[,ensembl_gene]) & is.na(feature_to_gene_map[,entrez_gene]))

# use RGD for all conversions 
# 11/15/21 RGD mapping
rgd = dl_read_gcp("nicolerg/MOTRPAC/RGD/GENES_RAT.txt", sep='\t')
rgd = unique(rgd[,.(GENE_RGD_ID, SYMBOL, NCBI_GENE_ID, OLD_SYMBOL, ENSEMBL_ID)])

# rows often contain semicolons, which complicates things 
rgd0 = rgd[!grepl(';', SYMBOL) & !grepl(';', NCBI_GENE_ID) & !grepl(';', ENSEMBL_ID)]
# for remaining rows, make it into long format
rgd1 = rgd[grepl(';', SYMBOL) | grepl(';', NCBI_GENE_ID) | grepl(';', ENSEMBL_ID)]
dt_list = list()
for(i in 1:nrow(rgd1)){
  entrez = unname(unlist(strsplit(as.character(rgd1[i, NCBI_GENE_ID]), ';')))
  ensembl = unname(unlist(strsplit(as.character(rgd1[i, ENSEMBL_ID]), ';')))
  symbol = unname(unlist(strsplit(as.character(rgd1[i, SYMBOL]), ';')))
  # expand all 
  symbol_rep = rep(rep(symbol, each=length(ensembl)),
                   each = length(entrez))
  entrez_rep = rep(entrez, each = length(symbol)*length(ensembl))
  ensembl_rep = rep(ensembl, length(symbol)*length(entrez))
  
  dt = data.table(GENE_RGD_ID = rgd1[i, GENE_RGD_ID],
                  SYMBOL = symbol_rep,
                  NCBI_GENE_ID = entrez_rep,
                  OLD_SYMBOL = rgd1[i, OLD_SYMBOL],
                  ENSEMBL_ID = ensembl_rep)
  if(!nrow(dt) == nrow(unique(dt))){
    print(i)
    print(dt)
    stop("Rows are not unique")
  }
  dt_list[[i]] = dt
}
rgd1_long = rbindlist(dt_list)
rgd_long = rbindlist(list(rgd0, rgd1_long))
rgd_long[,NCBI_GENE_ID := as.numeric(NCBI_GENE_ID)]
# rgd_long now contains a unique mapping between NCBI_GENE_ID, ENSEMBL_ID, SYMBOL per row 
# no row contains a semicolon in the NCBI_GENE_ID, ENSEMBL_ID, SYMBOL columns 

feature_to_gene_map_list_v2 = list()

# for proteomics, save entrez_gene 
prot = unique(feature_to_gene_map[assay %in% c("prot-pr","prot-ph","prot-ac","prot-ub"), .(feature_ID, entrez_gene)])
prot = merge(prot, rgd_long[!is.na(NCBI_GENE_ID)], by.x='entrez_gene', by.y='NCBI_GENE_ID', all.x=T)
setnames(prot, 
         c('GENE_RGD_ID','SYMBOL','OLD_SYMBOL','ENSEMBL_ID'),
         c('rgd_gene','gene_symbol','old_gene_symbol','ensembl_gene'))
feature_to_gene_map_list_v2[['prot']] = prot

# for RNA, ATAC, RRBS, save ensembl 
genom = unique(feature_to_gene_map[assay %in% c("transcript-rna-seq","epigen-atac-seq","epigen-rrbs"), 
                                   .(feature_ID, ensembl_gene, relationship_to_gene, custom_annotation)])
genom = merge(genom, rgd_long[!is.na(ENSEMBL_ID)], by.x='ensembl_gene', by.y='ENSEMBL_ID', all.x=T)
setnames(genom, 
         c('GENE_RGD_ID','SYMBOL','OLD_SYMBOL','NCBI_GENE_ID'),
         c('rgd_gene','gene_symbol','old_gene_symbol','entrez_gene'))
feature_to_gene_map_list_v2[['genomics']] = genom

# for immunoassay, save gene_symbol 
immuno = unique(feature_to_gene_map[assay %in% c("immunoassay"), .(feature_ID, gene_symbol)])
# fill in missing values 
for(s in immuno[!gene_symbol %in% rgd_long[,SYMBOL],gene_symbol]){
  new_symbol = unique(rgd_long[grepl(sprintf(";%s;|;%s$|^%s;|^%s$", s,s,s,s), OLD_SYMBOL), SYMBOL])
  if(length(new_symbol) == 1){
    immuno[gene_symbol==s, gene_symbol := new_symbol]
    print(c(s, new_symbol))
  }else{
    stop(sprintf("More than one gene symbol for %s. Not sure what to do: %s", s, paste(new_symbol, collapse=', ')))
  }
}

immuno = merge(immuno, rgd_long[!is.na(SYMBOL)], by.x='gene_symbol', by.y='SYMBOL', all.x=T)
immuno[is.na(GENE_RGD_ID)]

setnames(immuno, 
         c('GENE_RGD_ID','NCBI_GENE_ID','OLD_SYMBOL','ENSEMBL_ID'),
         c('rgd_gene','entrez_gene','old_gene_symbol','ensembl_gene'))
feature_to_gene_map_list_v2[['immunoassay']] = immuno

# ignore metab for now

# merge
feature_to_gene_map_v2 = rbindlist(feature_to_gene_map_list_v2, fill=T)

# how many features have no ID?
writeLines(feature_to_gene_map_v2[is.na(entrez_gene) & is.na(gene_symbol) & is.na(ensembl_gene), feature_ID])

# add some manual annotation
manual_cur = c('AP_004892.1'='Mt-nd1',
               'AP_004893.1'='Mt-nd2',
               'AP_004894.1'='Mt-co1',
               'AP_004895.1'='Mt-co2',
               'AP_004896.1'='Mt-atp8',
               'AP_004897.1'='Mt-atp6',
               'AP_004898.1'='Mt-co3',
               'AP_004899.1'='Mt-nd3',
               'AP_004900.1'='Mt-nd4l',
               'AP_004901.1'='Mt-nd4',
               'AP_004904.1'='Mt-cyb',
               'AP_004902.1'='Mt-nd5',
               'AP_004898.1'='Mt-co3',
               'AP_004900.1'='Mt-nd4l', 
               'AP_004895.1_Y218y'='Mt-co2',
               'AP_004896.1_S63s'='Mt-atp8',
               'AP_004902.1_S565s'='Mt-nd5',
               'AP_004896.1_T47t'='Mt-atp8',
               'AP_004896.1_S53s'='Mt-atp8',
               'AP_004902.1_T571t'='Mt-nd5',
               'AP_004895.1_S225s'='Mt-co2',
               'AP_004902.1_S567s'='Mt-nd5',
               'AP_004896.1_S30sT32t'='Mt-atp8',
               'AP_004902.1_S545s'='Mt-nd5',
               'AP_004896.1_K46k'='Mt-atp8',
               'AP_004896.1_K54k'='Mt-atp8',
               'AP_004896.1_K57k'='Mt-atp8',
               'AP_004902.1_K547k'='Mt-nd5')
# for(g in unique(unname(manual_cur))){
#   print(g)
#   r1 = rgd_long[SYMBOL==g]
#   if(nrow(r1) == 0){
#     print(rgd_long[grepl(sprintf(";%s;|;%s$|^%s;|^%s$", g,g,g,g), OLD_SYMBOL)])
#   }else{
#     print(r1[,SYMBOL])
#   }
# }

# add RGD for these
feature_to_gene_map_v2[feature_ID %in% names(manual_cur), gene_symbol := manual_cur[feature_ID]]
tmp = feature_to_gene_map_v2[is.na(rgd_gene) & !is.na(gene_symbol)]
tmp = merge(tmp[,.(feature_ID, gene_symbol)], rgd_long, by.x='gene_symbol', by.y='SYMBOL')
setnames(tmp, 
         c('GENE_RGD_ID','NCBI_GENE_ID','OLD_SYMBOL','ENSEMBL_ID'),
         c('rgd_gene','entrez_gene','old_gene_symbol','ensembl_gene'))
feature_to_gene_map_v2 = feature_to_gene_map_v2[!(is.na(rgd_gene) & !is.na(gene_symbol))]
feature_to_gene_map_v2 = rbindlist(list(feature_to_gene_map_v2, tmp), fill=T)

# replace '' with NA 
feature_to_gene_map_v2 = feature_to_gene_map_v2[,lapply(.SD, function(x){
  x[x==''] = NA 
  x
})]

# how many features have no ID?
writeLines(feature_to_gene_map_v2[is.na(entrez_gene) & is.na(gene_symbol) & is.na(ensembl_gene), feature_ID])
length(unique(feature_to_gene_map_v2[,feature_ID]))
feature_to_gene_map_v2 = feature_to_gene_map_v2[!is.na(entrez_gene) | !is.na(gene_symbol) | !is.na(ensembl_gene)]
length(unique(feature_to_gene_map_v2[,feature_ID]))
# # we lose a few features at this step: 
# XP_008758695.1
# XP_006223321.1
# XP_008758695.1_S499s
# XP_008758695.1_S213s
# XP_008758695.1_T215t
# XP_008758695.1_S420s

# how many features don't have good gene mapping?
nrow(feature_to_gene_map_v2[is.na(rgd_gene)])
length(unique(feature_to_gene_map_v2[, feature_ID]))
length(unique(feature_to_gene_map_v2[is.na(rgd_gene), feature_ID]))

unique(feature_to_gene_map_v2[is.na(rgd_gene) & !is.na(entrez_gene), entrez_gene]) # 225
unique(feature_to_gene_map_v2[is.na(rgd_gene) & !is.na(ensembl_gene), ensembl_gene]) # 157
length(unique(feature_to_gene_map_v2[is.na(rgd_gene) & !is.na(entrez_gene), feature_ID])) # 1480
length(unique(feature_to_gene_map_v2[is.na(rgd_gene) & !is.na(ensembl_gene), feature_ID])) # 413

# add back metabolites 
feature_to_gene_map_v3 = rbindlist(list(feature_to_gene_map_v2, 
                                        unique(feature_to_gene_map[assay=='metab',.(feature_ID, kegg_id)])), 
                                   fill=T)

############################################################################################################################## ~
# rename duplicate features
############################################################################################################################## ~

# from load-all-dea-rdata.R
all_timewise_dea = dl_load_gcp("gs://mawg-data/pass1b-06/merged/pass1b-06_all-sig-timewise-dea_20211018.RData", "all_timewise_dea")
my_universe = unique(all_timewise_dea[,.(feature_ID, assay_abbr, tissue_abbreviation, dataset, panel, selection_fdr, metabolite_refmet)])

# FIRST - ID the 91 repeated features and fix the names

# make new feature column 
my_universe[, feature := sprintf("%s;%s;%s", assay_abbr, tissue_abbreviation, feature_ID)]

# check for duplicate features
# this will be true because of how metabolites with high heterogeneity in the meta-regression were handled
# there are also some duplicate analytes in the cytokine panels 
check_duplicate_features = function(all_timewise){
  if(!"feature" %in% colnames(all_timewise)){
    all_timewise[, feature := sprintf("%s;%s;%s", assay_abbr, tissue_abbreviation, feature_ID)]
  }
  unique_features = unique(all_timewise[,.(feature, feature_ID, tissue_abbreviation, assay_abbr, panel, dataset, selection_fdr, metabolite_refmet)])
  n_rows = nrow(unique_features)
  n_features = length(unique(unique_features[,feature])) 
  if(n_rows > n_features){
    mm_features = unique(unique_features[duplicated(feature), feature_ID])
    writeLines(sprintf("%s duplicate feature identifiers.", length(mm_features)))
    return(mm_features)
  }else if(n_features == n_rows){
    writeLines("No duplicate feature identifiers.")
    return()
  }
}

multimeasured_features = check_duplicate_features(my_universe)
writeLines(multimeasured_features)

# make unique feature_IDs, but also keep the old ones
repeated_feature_map = unique(my_universe[feature_ID %in% multimeasured_features])
repeated_feature_map[feature_ID %in% multimeasured_features & !is.na(dataset), new_feature_ID := paste0(dataset, ":", feature_ID)]
repeated_feature_map[feature_ID %in% multimeasured_features & !is.na(panel), new_feature_ID := paste0(panel, ":", feature_ID)]
repeated_feature_map[, new_feature := sprintf("%s;%s;%s", assay_abbr, tissue_abbreviation, new_feature_ID)]
repeated_feature_map = unique(repeated_feature_map)
tmp_dup = copy(repeated_feature_map)
tmp_dup[,c('feature','feature_ID') := NULL]
setnames(tmp_dup, c('new_feature_ID','new_feature'), c('feature_ID','feature'))

# add duplicated feature name to my_universe
my_universe = rbindlist(list(my_universe, tmp_dup), use.names=T)

# add to master_feature_to_gene for duplicated features 
tmp_f2g = feature_to_gene_map_v3[feature_ID %in% multimeasured_features]
tmp_f2g = merge(tmp_f2g, repeated_feature_map[,.(feature_ID, new_feature_ID)], by='feature_ID')
tmp_f2g[,feature_ID := NULL]
setnames(tmp_f2g, 'new_feature_ID', 'feature_ID')

master_feature_to_gene = unique(rbindlist(list(feature_to_gene_map_v3, tmp_f2g), fill=T))
master_feature_to_gene[grepl("::", feature_ID)]

# check entrez, ensembl columns
master_feature_to_gene[!is.na(ensembl_gene) & !grepl("^ENSRNOG", ensembl_gene)]
master_feature_to_gene[!is.na(entrez_gene) & !grepl("^[0-9]+$", entrez_gene)]
master_feature_to_gene[,entrez_gene := as.numeric(entrez_gene)]

writeLines(sprintf("N total features: %s", length(unique(master_feature_to_gene[,feature_ID]))))
writeLines(sprintf("N features with gene symbols: %s", length(unique(master_feature_to_gene[!is.na(gene_symbol), feature_ID]))))
writeLines(sprintf("N features with Entrez IDs: %s", length(unique(master_feature_to_gene[!is.na(entrez_gene), feature_ID]))))
writeLines(sprintf("N features with Ensembl IDs: %s", length(unique(master_feature_to_gene[!is.na(ensembl_gene), feature_ID]))))
writeLines(sprintf("N features with symbols, Entrez IDs, and Ensembl IDs: %s", length(unique(master_feature_to_gene[!is.na(ensembl_gene) & !is.na(gene_symbol) & !is.na(entrez_gene), feature_ID])) ))
writeLines(sprintf("N features with symbols and Entrez IDs only: %s", length(unique(master_feature_to_gene[is.na(ensembl_gene) & !is.na(gene_symbol) & !is.na(entrez_gene), feature_ID]))))
writeLines(sprintf("N features with symbols and Ensembl IDs only: %s", length(unique(master_feature_to_gene[!is.na(ensembl_gene) & !is.na(gene_symbol) & is.na(entrez_gene), feature_ID]))))

### save master annotation, which includes all ATAC and RRBS features 
save(master_feature_to_gene, repeated_feature_map, file="master_feature_to_gene_20211116.RData")

############################################################################################################################## ~
# make universe ####
############################################################################################################################## ~

# now, make universe for each tissue and ome 

my_universe2 = merge(my_universe[,.(feature_ID, assay_abbr, tissue_abbreviation)], master_feature_to_gene, by='feature_ID')
my_universe2[, feature := sprintf("%s;%s;%s", assay_abbr, tissue_abbreviation, feature_ID)]
# my_universe2 has more lines than my_universe because master_feature_to_gene can have multiple lines per feature_ID

# remove sites that do not have a strong relationship with a gene
feature_to_gene = master_feature_to_gene[custom_annotation != 'Distal Intergenic' | is.na(custom_annotation)]
feature_to_gene = unique(feature_to_gene[,.(feature_ID, entrez_gene, gene_symbol, ensembl_gene, kegg_id)])

universes_list = list()
for(identifier in c('entrez_gene', 'gene_symbol', 'ensembl_gene', 'rgd_gene')){
  universe_list = list()
  for(assay in unique(my_universe2[,assay_abbr])){
    # for ATAC and RRBS, use expressed genes as the universe
    if(assay %in% c('ATAC','METHYL')){
      curr_universe = my_universe2[assay_abbr == 'TRNSCRPT']
    }else{
      curr_universe = my_universe2[assay_abbr == assay]
    }
    universe_list[[assay]] = list()
    for(tissue in unique(curr_universe[,tissue_abbreviation])){
      if(assay == 'METAB'){
        # always KEGG IDs for METAB
        ids = unique(curr_universe[tissue_abbreviation == tissue, kegg_id])
        ids = ids[!is.na(ids)]
        universe_list[[assay]][[tissue]] = ids 
      }else{
        ids = unique(curr_universe[tissue_abbreviation == tissue, get(identifier)])
        ids = ids[!is.na(ids)]
        universe_list[[assay]][[tissue]] = ids 
      }
    }
  }
  universes_list[[identifier]] = universe_list
}

lapply(universes_list, function(x) lapply(x, function(y) lapply(y, function(z) length(z))))

save(feature_to_gene, universes_list, repeated_feature_map, file="cluster-viz-inputs_20211116.RData")
system("gsutil cp cluster-viz-inputs_20211116.RData gs://mawg-data/pass1b-06/merged/")

# save universes separately
save(universes_list, file="universes_list_20211116.RData")
system("gsutil cp universes_list_20211116.RData gs://mawg-data/pass1b-06/merged/")
