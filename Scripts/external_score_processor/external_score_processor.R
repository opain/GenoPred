#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
library("optparse")

option_list = list(
	make_option("--ref_plink_chr", action="store", default=NULL, type='character',
		help="Path to per chromosome reference PLINK files [required]"),
	make_option("--pop_data", action="store", default=NULL, type='character',
		help="File containing the population code and location of the keep file [required]"),
	make_option("--plink2", action="store", default='plink2', type='character',
		help="Path PLINKv2 software binary [optional]"),
	make_option("--output", action="store", default=NULL, type='character',
		help="Path for output files [required]"),
	make_option("--score", action="store", default=NULL, type='character',
		help="Score file with format SNP A1 and then one or more effect sizes [required]"),
	make_option("--test", action="store", default=NA, type='character',
		help="Specify number of SNPs to include [optional]")
)

opt = parse_args(OptionParser(option_list=option_list))

# Load dependencies
library(GenoUtils)
library(data.table)
source('../functions/misc.R')
source_all('../functions')

# Check required inputs
if(is.null(opt$ref_plink_chr)){
  stop('--ref_plink_chr must be specified.\n')
}
if(is.null(opt$pop_data)){
  stop('--pop_data must be specified.\n')
}
if(is.null(opt$output)){
  stop('--output must be specified.\n')
}
if(is.null(opt$score)){
  stop('--score must be specified.\n')
}

# Create output directory
opt$output_dir <- paste0(dirname(opt$output), '/')
system(paste0('mkdir -p ', opt$output_dir))

# Create temp directory
tmp_dir <- tempdir()

# Initiate log file
log_file <- paste0(opt$output,'.log')
log_header(log_file = log_file, opt = opt, script = 'external_score_processor.R', start.time = start.time)

# If testing, change CHROMS to chr value
if(!is.na(opt$test) && opt$test == 'NA'){
  opt$test<-NA
}

if(!is.na(opt$test)){
  CHROMS <- as.numeric(gsub('chr','',opt$test))
}

####
# Read in score file
####

# Check for PGS Catalogue header format (https://www.pgscatalog.org/downloads/#scoring_header)
score_header <- readLines(opt$score, n = 40)
score_header <- score_header[grepl("^#", score_header)]
pgsc_header <- any(grepl('PGS CATALOG SCORING FILE', score_header))

# Use header to determine build if available
if(pgsc_header){
	if(any(grepl('HARMONIZATION', score_header))){
		target_build <- gsub('.*HmPOS_build=', '', score_header[grepl('HmPOS_build', score_header)])
	} else {
		target_build <- gsub('.*genome_build=', '', score_header[grepl('genome_build', score_header)])
	}

	# Change from hg and ncbi notation to GRCh
	if(tolower(target_build) == 'hg18' | tolower(target_build) == 'ncbi36') target_build <- 'GRCh36'
	if(tolower(target_build) == 'hg19' | tolower(target_build) == 'ncbi37') target_build <- 'GRCh37'
	if(tolower(target_build) == 'hg38' | tolower(target_build) == 'ncbi38') target_build <- 'GRCh38'
} else {
	target_build <- NA
}

# Read in the score file
score <- read_score(score = opt$score, chr = CHROMS, log_file = log_file)
n_snp_orig <- nrow(score)

log_add(log_file = log_file, message = paste0('Score file contains ',nrow(score),' variants.'))

#####
# Insert IUPAC codes into target
#####

# Insert IUPAC codes into target
score$IUPAC<-snp_iupac(score$A1, score$A2)

# Retain only non-ambiguous SNPs
score <- remove_ambig(score)

log_add(log_file = log_file, message = paste0('After removal of variants that are not SNPs or are ambiguous, ',nrow(score),' variants remain.'))

#####
# Harmonise per chromosome with reference
#####
# We identify variants present in score and reference, insert missing data (SNP, CHR, BP), and REF.FREQ
targ <- score
ref_rds <- opt$ref_plink_chr
chrs <- CHROMS

# Check whether CHR and BP information are present
chr_bp_avail<-sum(c('CHR','BP') %in% names(targ)) == 2

# Check whether RSIDs are available for majority of SNPs in GWAS
rsid_avail<-(sum(grepl('rs', targ$SNP)) > 0.9*length(targ$SNP))

targ_matched<-NULL
flip_logical_all<-NULL
if(chr_bp_avail){
	log_add(log_file = log_file, message = 'Merging sumstats with reference using CHR, BP, A1, and A2')

	if(!pgsc_header){
		###
		# Determine build
		###

		# If the score file is sparse, try using a few chromosomes until the build can be distinguished
		for(ref_chr in rev(sort(unique(targ$CHR)))){
			target_build <-
				detect_build(ref = readRDS(file = paste0(ref_rds, ref_chr, '.rds')),
											targ = targ[targ$CHR == ref_chr,])

			if(!is.na(target_build)){
				log_add(log_file = log_file, message = paste0('Genome build ', target_build,' detected.'))
				break
			}
		}
	}

	if(!is.na(target_build)){
		for(i in chrs){
			# Read reference data
			ref_i<-readRDS(file = paste0(ref_rds,i,'.rds'))

			# Retain only non-ambiguous SNPs
			ref_i<-remove_ambig(ref_i)

			# Rename columns prior to merging with target
			names(ref_i)<-paste0('REF.',names(ref_i))
			ref_i<-ref_i[, c('REF.CHR','REF.SNP',paste0('REF.BP_',target_build),'REF.A1','REF.A2','REF.IUPAC'), with=F]

			# Subset chromosome i from target
			targ_i<-targ[targ$CHR == i,]

			# Merge target and reference by BP
			ref_target<-merge(targ_i, ref_i, by.x='BP', by.y=paste0('REF.BP_',target_build))

			# Identify targ-ref strand flips, and flip target
			flip_logical<-detect_strand_flip(targ = ref_target$IUPAC, ref = ref_target$REF.IUPAC)
			flip_logical_all<-c(flip_logical_all, flip_logical)

			flipped<-ref_target[flip_logical,]
			flipped$A1<-snp_allele_comp(flipped$A1)
			flipped$A2<-snp_allele_comp(flipped$A2)
			flipped$IUPAC<-snp_iupac(flipped$A1, flipped$A2)

			# Identify SNPs that have matched IUPAC
			matched<-ref_target[ref_target$IUPAC == ref_target$REF.IUPAC,]
			matched<-rbind(matched, flipped)

			# Retain reference SNP and REF.FREQ data
			matched<-matched[, names(matched) %in% c('CHR','BP','A1','A2','effect_weight'), with=F]
			names(matched)[names(matched) == 'REF.SNP']<-'SNP'

			targ_matched<-rbind(targ_matched, matched)
		}
	}

	if(is.na(target_build) & !(rsid_avail)){
		log_add(log_file = log_file, message = 'Error: Target build could not be determined and SNP IDs unavailable.')
		stop('Target build could not be determined and SNP IDs unavailable.\n')
	}
	if(is.na(target_build) & rsid_avail){
		log_add(log_file = log_file, message = 'Target build could not be determined from CHR and BP data.')
	}
}

if(is.na(target_build) & rsid_avail){
	log_add(log_file = log_file, message = 'Using SNP, A1 and A2 to merge with the reference.')

	for(i in chrs){

		# Read reference data
		ref_i<-readRDS(file = paste0(ref_rds,i,'.rds'))

		# Retain only non-ambiguous SNPs
		ref_i<-remove_ambig(ref_i)

		# Rename columns prior to merging with target
		names(ref_i)<-paste0('REF.',names(ref_i))
		ref_i<-ref_i[, c('REF.CHR','REF.SNP','REF.BP_GRCh37','REF.A1','REF.A2','REF.IUPAC'), with=F]

		# Merge target and reference by SNP ID
		ref_target<-merge(targ, ref_i, by.x='SNP', by.y='REF.SNP')

		# Identify targ-ref strand flips, and flip target
		flip_logical<-detect_strand_flip(targ = ref_target$IUPAC, ref = ref_target$REF.IUPAC)
		flip_logical_all<-c(flip_logical_all, flip_logical)

		flipped<-ref_target[flip_logical,]
		flipped$A1<-snp_allele_comp(flipped$A1)
		flipped$A2<-snp_allele_comp(flipped$A2)
		flipped$IUPAC<-snp_iupac(flipped$A1, flipped$A2)

		# Identify SNPs that have matched IUPAC
		matched<-ref_target[ref_target$IUPAC == ref_target$REF.IUPAC,]
		matched<-rbind(matched, flipped)

		# Retain reference CHR and BP_GRCh37 data
		matched<-matched[, names(matched) %in% c('REF.CHR','REF.BP_GRCh37','SNP','A1','A2','effect_weight'), with=F]
		names(matched)[names(matched) == 'REF.CHR']<-'CHR'
		names(matched)[names(matched) == 'REF.BP_GRCh37']<-'BP'

		targ_matched<-rbind(targ_matched, matched)
	}
}

log_add(log_file = log_file, message = paste0('After matching variants to the reference, ',nrow(targ_matched),' variants remain.'))
log_add(log_file = log_file, message = paste0(sum(flip_logical_all), ' variants were flipped to match reference.'))

if(nrow(targ_matched) < 0.75*n_snp_orig){
	log_add(log_file = log_file, message = paste0("<75% of variants in score file are present in the reference (", nrow(targ_matched)," out of ", n_snp_orig, ")"))
	log_add(log_file = log_file, message = 'Skipping')

} else {
	####
	# Format as score file for GenoPred
	####

	score$CHR <- NULL
	score$BP <- NULL
	score <- score[, c('SNP','A1','A2','effect_weight'), with = F]
	names(score)[names(score) == 'effect_weight']<-paste0('SCORE_external')

	# Flip effects to match reference alleles
	ref <- read_pvar(opt$ref_plink_chr, chr = CHROMS)[, c('SNP','A1','A2'), with=F]
	score_new <- map_score(ref = ref, score = score)

	fwrite(score_new, paste0(opt$output,'.score'), col.names=T, sep=' ', quote=F)

	if(file.exists(paste0(opt$output,'.score.gz'))){
		system(paste0('rm ',opt$output,'.score.gz'))
	}

	system(paste0('gzip ',opt$output,'.score'))

	####
	# Calculate mean and sd of polygenic scores
	####

	log_add(log_file = log_file, message = 'Calculating polygenic scores in reference.')

	# Calculate scores in the full reference
	ref_pgs <- plink_score(pfile = opt$ref_plink_chr, chr = CHROMS, plink2 = opt$plink2, score = paste0(opt$output,'.score.gz'))

	# Calculate scale within each reference population
	pop_data <- fread(opt$pop_data)
	pop_data<-data.table(
	  FID=pop_data$`#IID`,
	  IID=pop_data$`#IID`,
	  POP=pop_data$POP
	)

	for(pop_i in unique(pop_data$POP)){
	ref_pgs_scale_i <- score_mean_sd(scores = ref_pgs, keep = pop_data[pop_data$POP == pop_i, c('FID','IID'), with=F])
	fwrite(ref_pgs_scale_i, paste0(opt$output, '-', pop_i, '.scale'), row.names = F, quote=F, sep=' ', na='NA')
	}
}

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = log_file, append = T)
cat('Analysis finished at', as.character(end.time),'\n')
cat('Analysis duration was', as.character(round(time.taken,2)), attr(time.taken, 'units'), '\n')
sink()
