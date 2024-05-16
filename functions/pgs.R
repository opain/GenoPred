#!/usr/bin/Rscript

# Identify score files (name and method combinations)
list_score_files <- function(config){

  combos<-NULL

  # Read in gwas_list
  gwas_list <- read_param(config = config, param = 'gwas_list')

  if(!is.null(gwas_list)){
    # Identify PGS methods to be included
    pgs_methods_list <- read_param(config = config, param = 'pgs_methods', return_obj = F)

    combos <- rbind(combos,
                    expand.grid(name = gwas_list$name[gwas_list$pop == 'EUR'], method = pgs_methods_list))

    # List PGS methods applied to non-EUR populations
    pgs_methods_noneur <- c('ptclump','lassosum','megaprs','prscs','dbslmm')
    pgs_methods_noneur <- pgs_methods_noneur[pgs_methods_noneur %in% pgs_methods_list]

    combos <- rbind(combos,
                    expand.grid(name = gwas_list$name[gwas_list$pop != 'EUR'], method = pgs_methods_noneur))
  }

  # Read in score_list
  score_list <- read_param(config = config, param = 'score_list')

  outdir <- read_param(config = config, param = 'outdir', return_obj = F)

  if(!is.null(score_list)){
    # Read in score_reporter output
    score_reporter <- fread(paste0(outdir, "/reference/pgs_score_files/external/score_report.txt"))
    score_list <- merge(score_list, score_reporter, by='name')

    # Remove scores that did not pass ref harmonisation
    score_list <- score_list[score_list$pass == T,]

    combos <- rbind(combos,
                    data.frame(
                      name = score_list$name,
                      method = 'external'))
  }

  return(combos)
}

# Flip effects in score file to match A1 reference
map_score<-function(ref, score){
  tmp <- merge(ref, score, by = 'SNP', all.x=T, sort = F)
  flip <- which(tmp$A1.x != tmp$A1.y)
  tmp <- as.matrix(tmp[, -1:-5, drop = FALSE])
  tmp[flip, ] <- -tmp[flip, ,drop=F]
  new_score<-cbind(ref, tmp)
  new_score[is.na(new_score)]<-0

  return(new_score)
}

# Calculate mean and sd of scores in file with plink .sscore format
score_mean_sd<-function(scores, keep=NULL){
    if(!is.null(keep)){
        scores<-scores[paste0(scores$FID, '_', scores$FID) %in% paste0(keep$FID, '_', keep$FID),]
    }

    scale<-data.table(  Param=names(scores)[-1:-2],
                        Mean=sapply(scores[,-1:-2, with=F], function(x) mean(x)),
                        SD=sapply(scores[,-1:-2, with=F], function(x) sd(x)))

    return(scale)
}

# Scale the target scores based on the reference mean and sd
score_scale<-function(score, ref_scale){
    score_scaled<-score
    for(i in ref_scale$Param){
        score_scaled[[i]]<-score[[i]]-ref_scale$Mean[ref_scale$Param == i]
        score_scaled[[i]]<-score_scaled[[i]]/ref_scale$SD[ref_scale$Param == i]
        score_scaled[[i]]<-round(score_scaled[[i]],3)
    }
    return(score_scaled)
}

# Read in SNP data from either plink1 binary, bgen or vcf
read_geno <- function(target, format) {
  if (!all(format %in% c('plink1', 'plink2', 'bgen', 'vcf'))) {
    stop('Specified format must be either plink1, plink2, bgen, or vcf.')
  }

  if (format == 'plink1') {
    target_snp <- fread(paste0(target, '.bim'))
    target_snp$V3 <- NULL
    names(target_snp) <- c('CHR', 'SNP', 'BP', 'A1', 'A2')
  }

  if (format == 'plink2') {
    target_snp <- fread(paste0(target, '.pvar'))
    names(target_snp)[1:5] <- c('CHR', 'BP', 'SNP', 'A1', 'A2')
  }

  if (format == 'bgen') {
    connection = dbConnect(RSQLite::SQLite(), paste0(target, '.bgen.bgi'))
    target_snp = dbGetQuery(connection, "SELECT * FROM Variant")
    target_snp <-
      target_snp[, c('chromosome', 'rsid', 'position', 'allele1', 'allele2')]
    names(target_snp) <- c('CHR', 'SNP', 'BP', 'A1', 'A2')
    dbDisconnect(connection)
    target_snp <- data.table(target_snp)
    target_snp$CHR <- as.numeric(gsub('chr', '', target_snp$CHR))
  }

  if (format == 'vcf') {
    target_snp <- fread(cmd = paste0("zcat ", target, ".vcf.gz | cut -f 1-5"), skip='#CHROM')
    names(target_snp) <- c('CHR', 'BP', 'SNP', 'A1', 'A2')
    target_snp$CHR <- as.numeric(gsub('chr', '', target_snp$CHR))
  }

  target_snp <- target_snp[, c('CHR', 'BP', 'SNP', 'A1', 'A2'), with = F]

  return(target_snp)
}

# Read in PLINK .bim file
read_bim<-function(dat, chr = 1:22){
  bim<-NULL
  for(i in chr){
    bim<-rbind(bim, fread(paste0(dat, i,'.bim')))
  }
  bim<-bim[,c('V1','V2','V4','V5','V6')]
  names(bim)<-c('CHR','SNP','BP','A1','A2')

  return(bim)
}

# Read in reference pop_data
read_pop_data <- function(x){
  pop_data<-fread(x)
  names(pop_data)<-gsub('\\#', '', names(pop_data))
  if (ncol(pop_data) == 2) {
    pop_data <- data.table(FID = pop_data$IID,
                           IID = pop_data$IID,
                           POP = pop_data$POP)
  } else {
    pop_data <- data.table(FID = pop_data$FID,
                           IID = pop_data$IID,
                           POP = pop_data$POP)
  }
  return(pop_data)
}

# Read in PLINK2 .pvar file
read_pvar<-function(dat, chr = 1:22){
  pvar<-NULL
  for(i in chr){
    pvar<-rbind(pvar, fread(paste0(dat, i,'.pvar')))
  }
  names(pvar)<-c('CHR','BP','SNP','A2','A1')

  return(pvar)
}

# Remove variants within genomic regions (REF: PMC2443852)
remove_regions<-function(dat, regions){
  exclude<-NULL
  for(i in 1:nrow(regions)){
    exclude<-c(exclude, dat$SNP[(   dat$CHR == regions$CHR[i]  &
                                    dat$BP >= regions$P0[i] &
                                    dat$BP <= regions$P1[i])])
  }

  return(dat[!(dat$SNP %in% exclude),])
}

# Read in GWAS summary statistics
read_sumstats<-function(sumstats, chr = 1:22, log_file = NULL, extract = NULL, req_cols = c('CHR','BP','A1','A2','BETA','SE','P','FREQ')){
  gwas<-fread(sumstats)

  log_add(log_file = log_file, message = paste0('sumstats contains ',nrow(gwas),' variants.'))

  # Retain requested chromosomes
  if(!all(1:22 %in% chr)){
    if(!('CHR' %in% names(gwas))){
      stop('Cannot filter sumstats by chromosome when CHR column is not present.')
    }
    gwas <- gwas[gwas$CHR %in% chr,]
    log_add(log_file = log_file, message = paste0(nrow(gwas), ' variants remain after selecting chromosomes.'))
  }

  # If FREQ is missing, use REF.FREQ
  if('FREQ' %in% req_cols){
    if(all(names(gwas) != 'FREQ')){
      names(gwas)[names(gwas) == 'REF.FREQ']<-'FREQ'
      log_add(log_file = log_file, message = 'REF.FREQ being used as FREQ.')
    }
  }

  # Check for essential columns
  if(!(all(req_cols %in% names(gwas)))){
    stop(paste0('Required column/s missing in sumstats: ', paste(!(names(gwas) %in% req_cols), collapse = ', ')))
  }

  # Remove non-essential columns
  gwas <- gwas[, req_cols, with = F]

  # Remove rows with missing data
  gwas <- gwas[complete.cases(gwas),]

  log_add(log_file = log_file, message = paste0('sumstats contains ',nrow(gwas),' variants with complete data.'))

  if(!is.null(extract)){
    extract <- obj_or_file(extract, return_file = F)
    gwas<-gwas[(gwas$SNP %in% extract),]
    log_add(log_file = log_file, message = paste0('After applying the extract file, ',nrow(gwas),' variants remain.'))
  }

  return(gwas)
}

# Run LDSC
ldsc <- function(sumstats, ldsc, hm3_snplist, munge_sumstats, ld_scores, pop_prev = NULL, sample_prev = NULL, log_file = NULL){
  tmp_dir<-tempdir()

  # Munge the sumstats
  system(paste0(munge_sumstats, ' --sumstats ', sumstats,' --merge-alleles ', hm3_snplist, ' --out ', tmp_dir,'/munged'))

  # Define the file paths
  sumstats_path <- file.path(tmp_dir,'/munged.sumstats.gz')
  output_path <- file.path(tmp_dir, '/ldsc')

  # Run the appropriate LDSC command based on the availability of prevalence data
  if(!is.null(pop_prev) && !is.null(sample_prev)) {
    system(paste0(ldsc, ' --h2 ', sumstats_path, ' --ref-ld ', ld_scores, ' --w-ld ', ld_scores, ' --out ', output_path, ' --samp-prev ', sample_prev, ' --pop-prev ', pop_prev))
    scale_type <- "Liability"
  } else {
    system(paste0(ldsc, ' --h2 ', sumstats_path, ' --ref-ld ', ld_scores, ' --w-ld ', ld_scores, ' --out ', output_path))
    scale_type <- "Observed"
  }

  # Read the log file and extract heritability
  ldsc_log <- readLines(paste0(tmp_dir, '/ldsc.log'))
  pattern <- paste0('Total ', scale_type, ' scale h2: ')
  ldsc_h2 <- gsub(pattern, '', ldsc_log[grepl(pattern, ldsc_log)])

  # Log the heritability estimate
  log_add(log_file = log_file, message = paste0('SNP-heritability estimate on the ', tolower(scale_type), ' scale = ', ldsc_h2, '.'))

  # Process and return the heritability estimate
  ldsc_h2 <- as.numeric(sub(' .*', '', ldsc_h2))
  if(!is.na(ldsc_h2) && ldsc_h2 > 1){
    ldsc_h2 <- 1
    log_add(log_file = log_file, message = paste0('Setting SNP-heritability estimate to 1.'))
  }

  return(ldsc_h2)
}

# Match A1 and A2 match a reference
allele_match<-function(sumstats, ref_bim, chr = 1:22){
  sumstats_ref<-merge(sumstats, ref_bim[, c('SNP','A1','A2'), with=F], by = 'SNP')
  swap_index <- sumstats_ref$A1.x == sumstats_ref$A2.y & sumstats_ref$A2.x == sumstats_ref$A1.y
  sumstats_ref$BETA[swap_index] <- -sumstats_ref$BETA[swap_index]
  sumstats_ref$FREQ[swap_index] <- 1-sumstats_ref$FREQ[swap_index]

  names(sumstats_ref)[names(sumstats_ref) == 'A1.y'] <- 'A1'
  names(sumstats_ref)[names(sumstats_ref) == 'A2.y'] <- 'A2'
  sumstats_ref<-sumstats_ref[, names(sumstats), with = F]

  return(sumstats_ref)
}

# Run DBSLMM
dbslmm <- function(dbslmm, plink = plink, ld_blocks, chr = 1:22, bfile, sumstats, h2, h2f = 1, nsnp, nindiv, log_file = NULL, ncores=1){
  # Create temp directory
  tmp_dir<-tempdir()

  # Make sure there is permission to run dbslmm
  system(paste0('chmod 777 ', dbslmm))

  # Run dbslmm for each chromosome
  for(chr_i in chr){
    cmd <-paste0('Rscript ', dbslmm,'/DBSLMM.R --plink ', plink,' --type t --block ', ld_blocks,'/fourier_ls-chr', chr_i, '.bed --dbslmm ', dbslmm, '/dbslmm --model DBSLMM --h2 ', h2,' --h2f ', paste(h2f, collapse = ','),' --ref ', bfile, chr_i, ' --summary ', sumstats, chr_i, '.assoc.txt --n ', nindiv,' --nsnp ', nsnp, ' --outPath ', tmp_dir, '/ --thread ', ncores)
    exit_status <- system(cmd, intern = F)
  }

  # Read in DBSLMM output
  for(h2f_i in h2f){
    dbslmm_output_i <- list.files(path=tmp_dir, pattern=paste0('h2f', h2f_i, '.dbslmm.txt'))
    if(length(dbslmm_output_i) != 22){
      log_add(log_file = log_file, message = paste0('At least one chromosome did not complete with h2f of', h2f_i, '.'))
    }

    dbslmm_all_i<-NULL
    for(file_i in dbslmm_output_i){
      dbslmm_all_i<-rbind(dbslmm_all_i, fread(paste0(tmp_dir,'/',file_i)))
    }

    if(h2f_i == h2f[1]){
      dbslmm_all<-dbslmm_all_i[,c(1,2,4), with=T]
      names(dbslmm_all)<-c('SNP', 'A1', paste0('SCORE_DBSLMM_',h2f_i))
    } else {
      dbslmm_all_i<-dbslmm_all_i[,c(1,4), with=T]
      names(dbslmm_all_i)<-c('SNP', paste0('SCORE_DBSLMM_',h2f_i))
      dbslmm_all <- merge(dbslmm_all, dbslmm_all_i, by = 'SNP')
    }
  }

  # Insert A2 information
  bim<-read_bim(bfile, chr = chr)
  dbslmm_all<-merge(dbslmm_all, bim[, c('SNP','A1','A2'), with = F], by = c('SNP','A1'), all.x = T)
  if(any(is.na(dbslmm_all$A2))){
    stop('Insertion of A2 in dbslmm score file failed.')
  }

  dbslmm_all <- dbslmm_all[, c('SNP','A1','A2',names(dbslmm_all)[grepl('SCORE_DBSLMM', names(dbslmm_all))]), with=F]

  return(dbslmm_all)
}

# Calculate predictor correlations using LDAK
ldak_pred_cor<-function(bfile, ldak, n_cores, chr = 1:22, keep = NULL){
  tmp_dir<-tempdir()
  tmp_file<-tempfile()

  ldak_opt <- NULL
  if(!is.null(keep)){
    ldak_opt <- paste0(ldak_opt, '--keep ', keep, ' ')
  }
  if(length(chr) != 1){
    for(chr_i in chr){
      system(paste0(ldak, ' --calc-cors ', tmp_dir, '/cors_full', chr_i, ' --bfile ', bfile, ' ', ldak_opt,'--window-cm 3 --chr ', chr_i, ' --max-threads ', n_cores))
    }
    write.table(paste0(tmp_dir, '/cors_full', chr), paste0(tmp_dir, '/cors_list.txt'), col.names=F, row.names=F, quote=F)
    system(paste0(ldak,' --join-cors ', tmp_file, ' --corslist ', tmp_dir, '/cors_list.txt --max-threads ', n_cores))
  } else {
    system(paste0(ldak,' --calc-cors ', tmp_file, ' --bfile ', bfile, ' ', ldak_opt,'--window-cm 3 --chr ', chr, ' --max-threads ', n_cores))
  }

  return(tmp_file)
}

# Read in external score file
read_score <- function(score, chr = 1:22, log_file = NULL){
	# Read in score file
	score <- fread(score)

	# Check whether harmonised columns present
	if(any(c("hm_rsID", "hm_chr", "hm_pos") %in% names(score))){
		log_add(log_file = log_file, message = 'PGSC harmonisation data present.')

    # If other_allele is not present in score file, try to use inferred other allele data in score file
    if(all(names(score) != 'other_allele')){
      score$other_allele <- score$hm_inferOtherAllele
      log_add(log_file = log_file, message = 'A2 is missing so using PGSC inferred A2.')
    }

		# Relabel header
		score <- score[, names(score) %in% c('hm_rsID','hm_chr','hm_pos','effect_allele','other_allele','effect_weight'), with = F]
		names(score)[names(score) == 'hm_rsID'] <- 'SNP'
		names(score)[names(score) == 'hm_chr'] <- 'CHR'
		names(score)[names(score) == 'hm_pos'] <- 'BP'
		names(score)[names(score) == 'effect_allele'] <- 'A1'
		names(score)[names(score) == 'other_allele'] <- 'A2'

	} else {
		# Relabel header
		score <- score[, names(score) %in% c('rsID','chr_name','chr_position','effect_allele','other_allele','effect_weight'), with = F]
		names(score)[names(score) == 'rsID'] <- 'SNP'
		names(score)[names(score) == 'chr_name'] <- 'CHR'
		names(score)[names(score) == 'chr_position'] <- 'BP'
		names(score)[names(score) == 'effect_allele'] <- 'A1'
		names(score)[names(score) == 'other_allele'] <- 'A2'
	}

  if(!('SNP' %in% names(score)) & !(all(c('CHR','BP') %in% names(score)))){
    log_add(log_file = log_file, message = 'Either SNP, or CHR and BP data must be present in the score file')
    stop('Either SNP, or CHR and BP data must be present in the score file.')
  }

  if(!(all(c('A1','A2') %in% names(score)))){
    log_add(log_file = log_file, message = 'Both A1 and A2 data must be present in the score file')
    stop('Both A1 and A2 data must be present in the score file.')
  }

	if(any(names(score) == 'CHR')){
	  score <- score[score$CHR %in% chr,]
	}

	return(score)
}

