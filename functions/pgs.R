#!/usr/bin/Rscript

# Identify score files (name and method combinations)
list_score_files <- function(config, quiet = F){

  combos<-NULL

  # Read in gwas_list
  gwas_list <- read_param(config = config, param = 'gwas_list', quiet = quiet)

  if(!is.null(gwas_list)){
    # Identify PGS methods to be included
    pgs_methods_list <- read_param(config = config, param = 'pgs_methods', return_obj = F, quiet = quiet)

    # Remove methods that are applied to groups of gwas
    pgs_methods_list <- pgs_methods_list[!(pgs_methods_list %in% pgs_group_methods)]

    combos <- rbind(combos,
                    expand.grid(name = gwas_list$name[gwas_list$pop == 'EUR'], method = pgs_methods_list))

    # List PGS methods applied to non-EUR populations
    pgs_methods_noneur <- pgs_methods_noneur[pgs_methods_noneur %in% pgs_methods_list]

    combos <- rbind(combos,
                    expand.grid(name = gwas_list$name[gwas_list$pop != 'EUR'], method = pgs_methods_noneur))
  }

  # Read in score_list
  score_list <- read_param(config = config, param = 'score_list', quiet = quiet)

  outdir <- read_param(config = config, param = 'outdir', return_obj = F, quiet = quiet)

  if(!is.null(score_list)){
    combos <- rbind(combos,
                    data.frame(
                      name = score_list$name,
                      method = 'external'))
  }

  # Read in gwas_groups
  gwas_groups <- read_param(config = config, param = 'gwas_groups', quiet = quiet)

  # Methods implemented when GWAS groups contains only 2 GWAS
  if(!is.null(gwas_groups)){
    # Identify gwas_groups containing >2 GWAS
    gwas_groups_2 <- gwas_groups[sapply(gwas_groups$gwas, function(x) sum(strsplit(x, ",")[[1]] != "") == 2), ]
    
    # Identify PGS methods to be included
    pgs_methods_list <- read_param(config = config, param = 'pgs_methods', return_obj = F, quiet = quiet)

    # Retain methods that are applied to groups with only 2 gwas
    pgs_methods_list <- pgs_methods_list[(pgs_methods_list %in% c('prscsx', 'xwing'))]
    
    # Provide combos for methods applied to groups of gwas
    combos <- rbind(combos, expand.grid(name = gwas_groups_2$name, method = pgs_methods_list))
    
    # For TL-PRS, list combos for tlprs_methods
    tlprs_methods<-read_param(config = config, param = 'tlprs_methods', return_obj = F, quiet = quiet)
    if(length(tlprs_methods) > 1 || !is.na(tlprs_methods)){
      combos <- rbind(combos, expand.grid(name = gwas_groups_2$name, method = paste0('tlprs_', tlprs_methods)))
    }
    
    # For LEOPARD, list combos for leopard_methods
    leopard_methods<-read_param(config = config, param = 'leopard_methods', return_obj = F, quiet = quiet)
    if(length(leopard_methods) > 1 || !is.na(leopard_methods)){
      combos <- rbind(combos, expand.grid(name = gwas_groups_2$name, method = paste0(leopard_methods,'_multi')))
    }
  }

  # Methods implemented when GWAS groups contain >2 GWAS
  if(!is.null(gwas_groups)){
    # Identify gwas_groups with more than 2 gwas
    gwas_groups_more <- gwas_groups[sapply(gwas_groups$gwas, function(x) sum(strsplit(x, ",")[[1]] != "") > 2), ]
    
    # Identify PGS methods to be included
    pgs_methods_list <- read_param(config = config, param = 'pgs_methods', return_obj = F, quiet = quiet)
    
    # Retain methods that are applied to groups with only 2 gwas
    pgs_methods_list <- pgs_methods_list[(pgs_methods_list %in% c('prscsx'))]
    
    # Provide combos for methods applied to groups of gwas
    combos <- rbind(combos, expand.grid(name = gwas_groups_more$name, method = pgs_methods_list))
    
    # For LEOPARD, list combos for leopard_methods
    leopard_methods<-read_param(config = config, param = 'leopard_methods', return_obj = F)
    if(length(leopard_methods) > 1 || !is.na(leopard_methods)){
      combos <- rbind(combos, expand.grid(name = gwas_groups_more$name, method = paste0(leopard_methods,'_multi')))
    }
  }
  
  combos <- data.table(combos)
  combos <- combos[, lapply(.SD, as.character)]
  
  # Remove score files where .score.gz does not exist
  combos_keep <- NULL
  for(i in 1:nrow(combos)){
    score_i <- paste0(outdir, '/reference/pgs_score_files/', combos$method[i],'/', combos$name[i],'/ref-',combos$name[i], '.score.gz')
    if(file.exists(score_i)){
      combos_keep <- rbind(combos_keep, combos[i,])
    } else {
      if(quiet == F){
        cat0('No score file present for ', combos$method[i],' - ', combos$name[i],'. Check logs for reason.\n')
      }
    }
  }
  
  return(combos_keep)
}

# Flip effects in score file to match A1 reference
map_score<-function(ref, score){
  # Check if required columns exist
  required_cols <- c('SNP', 'A1', 'A2')
  if (!all(required_cols %in% names(ref)) | 
      !all(required_cols %in% names(score))) {
    stop('ref and score must contain SNP, A1, and A2 columns.')
  }
  
  # Valid alleles
  valid_alleles <- c('A', 'T', 'C', 'G')
  
  # Check for NA or invalid alleles in A1 and A2 columns
  for (col in c('A1', 'A2')) {
    if (any(is.na(ref[[col]])) | any(!ref[[col]] %in% valid_alleles)) {
      stop(paste('Invalid allele values detected in ref column:', col))
    }
    if (any(is.na(score[[col]])) | any(!score[[col]] %in% valid_alleles)) {
      stop(paste('Invalid allele values detected in score column:', col))
    }
  }
  
  # Check for NA values in SNP column
  if (any(is.na(ref$SNP)) | any(is.na(score$SNP))) {
    stop('NA values detected in SNP column of ref or score.')
  }
  
  ref <- ref[, c('SNP','A1','A2'), with = F]
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
        scores<-scores[paste0(scores$FID, '_', scores$IID) %in% paste0(keep$FID, '_', keep$IID),]
    }

    scale<-data.table(  Param=names(scores)[-1:-2],
                        Mean=sapply(scores[,-1:-2, with=F], function(x) mean(x)),
                        SD=sapply(scores[,-1:-2, with=F], function(x) sd(x)))

    return(scale)
}

# Scale the target scores based on the reference mean and sd
score_scale<-function(score, ref_scale){
    for(i in ref_scale$Param){
        score[[i]]<-score[[i]]-ref_scale$Mean[ref_scale$Param == i]
        score[[i]]<-score[[i]]/ref_scale$SD[ref_scale$Param == i]
        score[[i]]<-round(score[[i]],3)
    }
    return(score)
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

# Read in the .freq file for target population
read_frq<-function(freq_dir, population, chr){
  freq_data<-NULL
  for(i in chr){
    tmp<-fread(paste0(freq_dir,'/',population,'/ref.',population,'.chr',i,'.afreq'))
    tmp<-data.table(
      SNP = tmp$ID,
      A1 = tmp$ALT,
      A2 = tmp$REF,
      FREQ = tmp$ALT_FREQS
    )
    freq_data<-rbind(freq_data, tmp)
  }
  return(freq_data)
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
      gwas$FREQ <- gwas$REF.FREQ
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
ldsc <- function(sumstats, ldsc, hm3_snplist, munge_sumstats, ld_scores, log_file = NULL){
  tmp_dir<-tempdir()

  # Munge the sumstats
  system(paste0(munge_sumstats, ' --sumstats ', sumstats,' --merge-alleles ', hm3_snplist, ' --out ', tmp_dir,'/munged'))

  # Define the file paths
  sumstats_path <- file.path(tmp_dir,'/munged.sumstats.gz')
  output_path <- file.path(tmp_dir, '/ldsc')

  # Run LDSC
  system(paste0(ldsc, ' --h2 ', sumstats_path, ' --ref-ld ', ld_scores, ' --w-ld ', ld_scores, ' --out ', output_path))

  # Read the log file
  ldsc_log <- readLines(paste0(tmp_dir, '/ldsc.log'))
  
  # Extract Total Observed h² and SE
  h2_line <- grep("^Total Observed scale h2:", ldsc_log, value = TRUE)
  h2_parts <- regmatches(h2_line, regexec("h2:\\s+([0-9.]+)\\s+\\(([0-9.]+)\\)", h2_line))[[1]]
  h2 <- as.numeric(h2_parts[2])
  h2_se <- as.numeric(h2_parts[3])
  
  # Extract lambda GC
  lambda_line <- grep("^Lambda GC:", ldsc_log, value = TRUE)
  lambda_gc <- as.numeric(sub("Lambda GC:\\s+", "", lambda_line))
  
  # Extract intercept and its SE
  intercept_line <- grep("^Intercept:", ldsc_log, value = TRUE)
  intercept_parts <- regmatches(intercept_line, regexec("Intercept:\\s+([0-9.]+)\\s+\\(([0-9.]+)\\)", intercept_line))[[1]]
  intercept <- as.numeric(intercept_parts[2])
  intercept_se <- as.numeric(intercept_parts[3])
  
  # Return as a data frame
  results <- data.frame(
    h2 = h2,
    h2_se = h2_se,
    lambda_gc = lambda_gc,
    intercept = intercept,
    intercept_se = intercept_se
  )
  
  # Log the heritability estimate
  log_add(log_file = log_file, message = paste0('SNP-heritability estimate on the observed scale = ', results$h2, " (", results$h2_se, ")."))
  log_add(log_file = log_file, message = paste0('LDSC intercept = ', results$intercept, " (", results$intercept_se, ")."))
  log_add(log_file = log_file, message = paste0('Lambda GC = ', results$lambda_gc, "."))
  
  return(results)
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
	  # Remove for 'chr' string in CHR column
	  score$CHR <- gsub('chr', '', score$CHR)
	  score <- score[score$CHR %in% chr,]
	}

	return(score)
}


quickprs<-function(sumstats, quickprs_ldref, quickprs_multi_ldref = NULL, genomic_control, prs_model, n_cores = 1, ref_subset = NULL){
  tmp_dir<-tempfile()
  dir.create(tmp_dir)
  
  # Check if quickprs_multi_ldref and ref_subset are both NULL or both non-NULL
  if (xor(is.null(quickprs_multi_ldref), is.null(ref_subset))) {
    stop("Both 'quickprs_multi_ldref' and 'ref_subset' must either be NULL or non-NULL.")
  }
  
  ######
  # Estimate Per-Predictor Heritabilities
  ######

  # Calculate Per-Predictor Heritabilities.
  quickprs_ldref_files<-list.files(quickprs_ldref)

  tagging_file<-quickprs_ldref_files[grepl('quickprs.tagging',quickprs_ldref_files)]
  matrix_file<-quickprs_ldref_files[grepl('quickprs.matrix',quickprs_ldref_files)]

  if(opt$genomic_control == F){
    system(paste0(opt$ldak,' --sum-hers ', tmp_dir, '/bld.ldak --tagfile ', quickprs_ldref, '/', tagging_file, ' --summary ', sumstats, ' --matrix ', quickprs_ldref, '/', matrix_file, ' --max-threads ', n_cores, ' --check-sums NO'))
  } else{
    system(paste0(opt$ldak,' --sum-hers ', tmp_dir, '/bld.ldak --genomic-control YES --tagfile ', quickprs_ldref, '/', tagging_file, ' --summary ', sumstats, ' --matrix ', quickprs_ldref, '/', matrix_file, ' --max-threads ', n_cores, ' --check-sums NO'))
  }

  ldak_res_her<-fread(paste0(tmp_dir,'/bld.ldak.hers'))

  ######
  # Estimate effect sizes for training and full prediction models.
  ######

  if(!is.null(ref_subset)){
    quickprs_multi_ldref_files<-list.files(quickprs_multi_ldref)
    ref_dir <- quickprs_multi_ldref
    cor_file_prefix<-gsub('.cors.bin','',quickprs_multi_ldref_files[grepl(paste0('subset_', ref_subset, '.cors.bin'),quickprs_multi_ldref_files)])
  } else {
    cor_file_prefix<-gsub('.cors.bin','',quickprs_ldref_files[grepl('.cors.bin',quickprs_ldref_files) & !grepl('subset', quickprs_ldref_files)])
    ref_dir <- quickprs_ldref
  }

  system(paste0(opt$ldak,' --mega-prs ',tmp_dir,'/mega_full --model ', prs_model,' --cors ', ref_dir, '/', cor_file_prefix, ' --ind-hers ', tmp_dir, '/bld.ldak.ind.hers --summary ', sumstats, ' --high-LD ', quickprs_ldref, '/highld.snps --cv-proportion 0.1 --window-cm 1 --max-threads ', n_cores,' --extract ', sumstats))

  # Identify the best fitting model
  ldak_res_cors <- fread(paste0(tmp_dir, '/mega_full.cors'), nThread = n_cores)
  best_score <- ldak_res_cors[which.max(ldak_res_cors$Correlation),]

  ######
  # Format final score file
  ######

  # Read in the scores
  score <- fread(paste0(tmp_dir,'/mega_full.effects'), nThread = n_cores)
  score <- score[, c(1, 2, 3, 5), with = F]
  names(score) <- c('SNP','A1','A2','BETA')

  return(score)
}

# Derive trans-ancestry PGS models and estimate PGS residual scale
model_trans_pgs<-function(scores=NULL, pcs=NULL, output=NULL){
  if(any(is.null(c(scores, pcs, output)))){
    stop('Error: All parameters must be specified.')
  }
  
  if(is.character(pcs)){
    # Read in the reference PCs, extract PC columns, and update headers
    pcs_dat<-fread(pcs)
    names(pcs_dat)[1]<-'FID'
    pcs_dat<-pcs_dat[,grepl('FID|IID|^PC', names(pcs_dat)), with=F]
  } else {
    pcs_dat<-pcs
  }
  
  # Merge PGS and PCs
  scores_pcs<-merge(scores, pcs_dat, by=c('FID','IID'))
  
  # Calculate PGS residuals
  pcs_noid<-scores_pcs[,grepl('^PC', names(scores_pcs)), with=F]
  
  mod_list<-NULL
  scores_pcs_resid<-scores_pcs
  for(i in names(scores)[-1:-2]){
    mod_list[[i]]<-list()
    
    tmp<-data.table(y=scores_pcs[[i]], pcs_noid)
    
    # Model differences in mean
    pgs_pc_mean_mod<-lm(y ~ ., data=tmp)
    
    # Model differences in variance of residuals
    # Use gamma distribution to constrain predicted variance to be non-negative
    predicted_pgs <- predict(pgs_pc_mean_mod, newdata = tmp)
    residual_pgs <- tmp$y - predicted_pgs
    squared_residuals <- residual_pgs^2
    squared_residuals <- pmax(squared_residuals, 1e-4)

    pgs_pc_var_mod <- glm2(squared_residuals ~ ., data = tmp[, names(tmp) != 'y', with=F], family = Gamma(link = "log"))
    predicted_pgs_var <- exp(predict(pgs_pc_var_mod, newdata = tmp))
    
    scores_pcs_resid[[i]]<-residual_pgs/sqrt(predicted_pgs_var)
    
    mod_list[[i]]$mean_model <- compact_lm(pgs_pc_mean_mod)
    mod_list[[i]]$var_model <- compact_lm(pgs_pc_var_mod)
  }
  
  scores_pcs_resid<-scores_pcs_resid[,grepl('FID|IID|^SCORE', names(scores_pcs_resid)), with=F]
  
  # Save mean and SD of PGS residuals in 'trans' population
  # This should be approximately mean = 0 and SD = 1, but save as a sanity check
  scores_pcs_resid_scale<-score_mean_sd(scores=scores_pcs_resid)
  fwrite(scores_pcs_resid_scale, paste0(output,'-TRANS.scale'), sep=' ', col.names=T, quote=F)
  
  # Save PGS ~ PC models
  saveRDS(mod_list, file = paste0(output,'-TRANS.model.rds'))
  
  # Save TRANS PGS in reference
  fwrite(scores_pcs_resid, paste0(output,'-TRANS.profiles'), sep=' ', na='NA', quote=F)
}

# Remove unused parts of model object for prediction
compact_lm <- function(cm) {
  # just in case we forgot to set
  # y=FALSE and model=FALSE
  cm$y = c()
  cm$model = c()
  
  cm$residuals = c()
  cm$fitted.values = c()
  cm$effects = c()
  cm$qr$qr = c()
  cm$linear.predictors = c()
  cm$weights = c()
  cm$prior.weights = c()
  cm$data = c()
  cm
}

# Adjust PGS for ancestry using reference PC models with parallel processing
score_adjust <- function(score, pcs, ref_model, chunk_size = 10) {
  original_order <- names(score)
  
  # Ensure 'pcs' keyed for fast lookup
  setkey(pcs, FID, IID)
  
  # List of score columns
  score_cols <- setdiff(names(score), c("FID", "IID"))
  
  # Split into chunks for memory efficiency
  score_chunks <- split(score_cols, ceiling(seq_along(score_cols) / chunk_size))
  
  # Match PC rows by FID/IID directly (minimal memory usage)
  matched_idx <- pcs[score[, .(FID, IID)], which = TRUE, nomatch = 0]
  
  # Process each chunk sequentially to manage RAM
  for (chunk in score_chunks) {
    
    # Run in parallel across columns within chunk
    chunk_results <- mclapply(chunk, function(col_name) {
      cat("Processing:", col_name, "\n")
      
      # Fetch models for current score column
      mean_model <- ref_model[[col_name]]$mean_model
      var_model  <- ref_model[[col_name]]$var_model
      
      # Pre-allocate adjusted score vector
      adjusted_score <- rep(NA_real_, nrow(score))
      
      # Predict mean and variance using matched PCs
      if (length(matched_idx) > 0) {
        adjusted_score[matched_idx] <- round(
          (score[[col_name]][matched_idx] - predict(mean_model, newdata = pcs[matched_idx])) / 
            sqrt(exp(predict(var_model, newdata = pcs[matched_idx]))),
          3
        )
      }
      
      adjusted_score
    }, mc.cores = min(getDoParWorkers(), 10))
    
    # Write results directly back to 'score' object by reference
    for (idx in seq_along(chunk)) {
      set(score, j = chunk[idx], value = chunk_results[[idx]])
    }
    
    # Explicitly clean memory after each chunk
    rm(chunk_results)
    gc()
    
  }
  
  setcolorder(score, original_order)
  
  return(score)
}

# Helper function to calculate relative weights for a single file from LEOPARD
cal_avg_rel_weights <- function(path){
  weights_file <- fread(path)
  weights_non_0 <-  ifelse(weights_file$Weights < 0, 0, weights_file$Weights)
  rel_weights <- weights_non_0/(sum(weights_non_0))
  return(rel_weights)
}

# Function to calculate average weights across replications from LEOPARD
calculate_avg_weights <- function(populations, leopard_dir, log_file = NULL) {
  avg_weights <- list()
  
  # Iterate over populations
  for (targ_pop in populations) {
    # Define the file prefix for the population
    weights_prefix <- paste0(leopard_dir, '/weights_', targ_pop, '/output_LEOPARD_weights_rep')
    
    # Calculate relative weights for each replication
    rel_weights_list <- lapply(1:4, function(i) {
      cal_avg_rel_weights(paste0(weights_prefix, i, ".txt"))
    })
    
    # Calculate the average weights across replications
    avg_weights[[targ_pop]] <- as.numeric(colMeans(do.call(rbind, rel_weights_list), na.rm = TRUE))
  }
  
  log_add(log_file = log_file, message = '------------------------')
  for(i in names(avg_weights)){
    log_add(log_file = log_file, message = paste0("LEOPARD weights - ", i, " target: "))
    for(j in populations){
      log_add(log_file = log_file, message = paste0(j, ' = ', avg_weights[[i]][which(populations == j)]))
    }
    log_add(log_file = log_file, message = '------------------------')
  }
  
  return(avg_weights)
}

# Centre SNP-weights
centre_weights <- function(score, freq, ref){
  # Sort and flip according to reference data
  score <- map_score(ref = ref, score = score)
  
  # Sort and flip freq according to reference just in case they are different
  freq <- map_score(ref = ref, score = freq)
  
  # Calculate mean genotype dosage and denominator
  freq$MeanGenotype <- 2 * freq$FREQ
  denominator <- sum(freq$MeanGenotype^2)
  
  for(i in names(score)[!(names(score) %in% c('SNP','A1','A2'))]){
    # Calculate mean of PGS
    mean_pgs <- sum(score[[i]] * freq$MeanGenotype)
    
    # Adjust the SNP-weights so PGS is centered
    score[[i]] <- score[[i]] - (mean_pgs / denominator) * freq$MeanGenotype
  }
  return(score)
}

# Linearly combine scores using mixing weights for target population
calculate_weighted_scores <- function(score, targ_pop, mix_weights) {
  if(!all((names(score) %in% c('SNP','A1','A2', paste0('SCORE_targ_', names(mix_weights)))))){
    stop(paste0('score should only contain columns SNP, A1, A2, ', paste(paste0('SCORE_targ_', names(mix_weights)), collapse=', ')))
  }
  score_weighted<-score
  for(disc_pop in names(mix_weights)){
    score_tmp <- score[[paste0('SCORE_targ_', disc_pop)]]
    weight_tmp <- mix_weights[[targ_pop]][which(names(mix_weights) == disc_pop)]
    score_weighted[[paste0('SCORE_targ_', disc_pop)]] <- score_tmp * weight_tmp
  }
  score_combined <- rowSums(score_weighted[, grepl('SCORE_', names(score_weighted)), with = FALSE])
  
  return(score_combined)
}

# Adjust weights to correspond to PGS with SD of 1
adjust_weights <- function(weights, pgs_sd) {
  adjusted_weights <- weights * pgs_sd
  # Normalize weights so they sum to 1
  normalized_weights <- adjusted_weights * (1 / sum(adjusted_weights))
  return(normalized_weights)
}

# Create function to run LRT on allele frequencies
lrt_af_dual <- function(p1, n1, p0, n0) {
  # Convert allele frequencies to counts of alternate alleles
  k1 <- round(2 * n1 * p1)
  k0 <- round(2 * n0 * p0)
  N1 <- 2 * n1
  N0 <- 2 * n0
  
  # Estimate common allele frequency under null
  p_common <- (k1 + k0) / (N1 + N0)
  
  # Log-likelihood under null: same freq
  logL0 <- k1 * log(p_common) + (N1 - k1) * log(1 - p_common) +
    k0 * log(p_common) + (N0 - k0) * log(1 - p_common)
  
  # Log-likelihood under alternative: separate freqs
  logL1 <- k1 * log(p1) + (N1 - k1) * log(1 - p1) +
    k0 * log(p0) + (N0 - k0) * log(1 - p0)
  
  stat <- 2 * (logL1 - logL0)
  pval <- pchisq(stat, df = 1, lower.tail = FALSE)
  
  return(list(stat = stat, p = pval))
}

# New version of lassosum p2cor function that allows for increased precision
z2cor<-function(z, n){
  t <- qt(pnorm(abs(z), lower.tail = FALSE, log.p = TRUE), df = n,
          lower.tail = FALSE, log.p = TRUE) * sign(z)
  return(t/sqrt(n - 2 + t^2))
}

