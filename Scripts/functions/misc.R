#!/usr/bin/Rscript

# Make a data.frame giving labels to the 1KG reference populations
pop_1kg <- data.frame(
  pop = c('AFR','AMR','EAS','EUR','SAS'),
  label = c('African','Admixed American','East Asian','European','South Asian')
)

# Make a data.frame giving labels to the 1KG reference populations
pgs_method_labels <- data.frame(
  method = c('ptclump','dbslmm','ldpred2','sbayesr','lassosum','prscs','megaprs','external'),
  label = c('pT+clump','DBSLMM','LDpred2','SBayesR','lassosum','PRS-CS','MegaPRS','External')
)
pgs_method_labels[order(pgs_method_labels$method),]

# Calculate scores in target file
calc_score<-function(bfile, score, keep=NULL, extract=NULL, chr=1:22, frq=NULL, plink2='plink2', threads=1){
    # Create object indicating tmpdir
    tmp_folder<-tempdir()

    # Determine the number of scores
    score_small<-fread(score, nrows=5)
    n_scores<-ncol(score_small)-3

    # Assemble command and files for keep and extract
    cmd=NULL
    if(!is.null(keep)){
        cmd<-paste0(cmd, ' --keep ', keep)
    }
    if(!is.null(extract)){
        cmd<-paste0(cmd, ' --extract ', extract)
    }
    if(!is.null(frq)){
        cmd<-paste0(cmd, ' --read-freq ', frq,'CHROMOSOME_NUMBER.frq --score ',score,' header-read')
    } else {
        cmd<-paste0(cmd, ' --score ',score,' header-read no-mean-imputation')
    }
    if(n_scores > 1){
        cmd<-paste0(cmd, ' --score-col-nums 4-',3+n_scores)
    } else {
        cmd<-paste0(cmd, ' --score-col-nums 4')
    }
    # Calculate score in the target sample
    for(i in chr){
        cmd_i<-gsub('CHROMOSOME_NUMBER',i,cmd)
        cmd_full<-paste0(plink2, ' --bfile ',bfile,i,cmd_i,' --chr ',i,' --out ',tmp_folder,'/profiles.chr',i,' --threads ',threads)
        exit_status <- system(cmd_full, intern=FALSE)
        if (exit_status == 2) {
          stop()
        }
    }

    # Add up the scores across chromosomes
    # Read in score files IDs column from first non-missing chromosome
    for(i in chr){
        if(file.exists(paste0(tmp_folder,'/profiles.chr',i,'.sscore'))){
            scores_ids<-fread(paste0(tmp_folder,'/profiles.chr',i,'.sscore'))[,1:2, with=F]
            names(scores_ids)<-c('FID','IID')
            break
        }
    }

    # Read in the scores for each chromosome, adjust for the number of SNPs considered and add up
    scores<-list()
    for(i in chr){
        if(file.exists(paste0(tmp_folder,'/profiles.chr',i,'.sscore'))){
            sscore<-fread(paste0(tmp_folder,'/profiles.chr',i,'.sscore'))
            scores[[i]]<-sscore[,grepl('_AVG$', names(sscore)),with=F]
            # This allows for difference plink formats across plink versions
            if(any(names(sscore) == 'NMISS_ALLELE_CT')){
                scores[[i]]<-as.matrix(scores[[i]]*sscore$NMISS_ALLELE_CT[1]/2)
            } else {
                scores[[i]]<-as.matrix(scores[[i]]*sscore$ALLELE_CT[1]/2)
            }
        } else {
            cat0('No scores for chromosome ',i,'. Check plink logs file for reason.\n')
        }
    }

    # Remove NULL elements from list (these are inserted by R when list objects are numbered)
    scores[sapply(scores, is.null)] <- NULL

    # sum scores across chromosomes
    scores<-Reduce(`+`, scores)

    # Combine score with IDs
    scores<-data.table(scores_ids,
                       scores)

    # Rename columns
    names(scores)<-c('FID','IID',names(score_small)[-1:-3])

    return(scores)
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

# Make a version of cat with sep=''
cat0 <- function(...) {
  cat(..., sep = '')
}

# Read in SNP data from either plink1 binary, bgen or vcf
read_geno <- function(target, format) {
  if (!all(format %in% c('plink1', 'bgen', 'vcf'))) {
    stop('Specified format must be either plink1, bgen, or vcf.')
  }

  if (format == 'plink1') {
    target_snp <- fread(paste0(target, '.bim'))
    target_snp$V3 <- NULL
    names(target_snp) <- c('CHR', 'SNP', 'BP', 'A1', 'A2')
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
    target_snp <- fread(cmd = paste0("zcat ", target, ".vcf.gz | cut -f 1-5"))
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

# Detect wether object refers to a data.frame or a file
obj_or_file<-function(x, header = F, return_file = T){
  if(is.null(x)){
    output <- NULL
  } else {
    if(is.vector(x)){
      if(length(x) == 1 && file.exists(x)){
        if(return_file){
          output <- x
        } else {
          output <- fread(x, header = header)
          if(ncol(output) == 1){
            output <- output$V1
          }
        }
      } else {
        if(return_file){
          output<-tempfile()
          write.table(x, output, col.names = header, row.names=F, quote=F)
        } else {
          output <- x
        }
      }
    } else if(is.data.frame(x)){
        if(return_file){
          output<-tempfile()
          write.table(x, output, col.names = header, row.names=F, quote=F)
        } else {
          output <- x
        }
    } else {
      stop("Input must be a vector or data frame object, or a valid file path.")
    }
  }
  return(output)
}

# Make a subset of plink1 binaries
plink_subset<-function(plink='plink', chr = 1:22, keep = NULL, extract = NULL, bfile, out, memory = 4000, threads = 1){

  # Prepare plink options
  plink_opt<-NULL
  if(!is.null(keep)){
    keep <- obj_or_file(keep)
    plink_opt<-paste0(plink_opt, paste0('--keep ', keep, ' '))
  }
  if(!is.null(extract)){
    extract <- obj_or_file(extract)
    plink_opt<-paste0(plink_opt, paste0('--extract ', extract, ' '))
  }

  # Run plink
  for(chr_i in chr){
    cmd <- paste0(plink,' --bfile ',bfile,chr_i,' --threads ', threads,' ',plink_opt,'--make-bed --out ',out,chr_i,' --memory ',memory)
    exit_status <- system(cmd, intern=FALSE)
    if (exit_status == 2) {
      stop()
    }
  }
}

# Return a list of variants surviving QC using plink
plink_qc_snplist<-function(bfile, plink = 'plink', keep = NULL, chr = 1:22, threads = 1, memory = 4000, geno = NULL, maf = NULL, hwe = NULL){
  plink_opt<-NULL
  if(!is.null(geno)){
    plink_opt <- paste0(plink_opt, paste0('--geno ',geno,' '))
  }
  if(!is.null(maf)){
    plink_opt <- paste0(plink_opt, paste0('--maf ',maf,' '))
  }
  if(!is.null(hwe)){
    plink_opt <- paste0(plink_opt, paste0('--hwe ',hwe,' '))
  }
  if(!is.null(keep)){
    keep<-obj_or_file(keep)
    plink_opt <- paste0(plink_opt, paste0('--keep ',keep,' '))
  }

  temp_file <- tempfile()
  snplist <- NULL
  for(chr_i in chr){
    system(paste0(plink, ' --bfile ', bfile, chr_i, ' --threads ', threads, ' ', plink_opt, ' --write-snplist --out ', temp_file,' --memory ', memory))
    snplist<-c(snplist, fread(paste0(temp_file, '.snplist'), header=F)$V1)
  }

  return(snplist)
}

# Coordinates of high ld regions in build GRCh37
long_ld_coord<-do.call(rbind, list(
    data.frame(CHR = 1, P0 = 48e6, P1 = 52e6),
    data.frame(CHR = 2, P0 = 86e6, P1 = 100.5e6),
    data.frame(CHR = 2, P0 = 134.5e6, P1 = 138e6),
    data.frame(CHR = 2, P0 = 183e6, P1 = 190e6),
    data.frame(CHR = 3, P0 = 47.5e6, P1 = 50e6),
    data.frame(CHR = 3, P0 = 83.5e6, P1 = 87e6),
    data.frame(CHR = 3, P0 = 89e6, P1 = 97.5e6),
    data.frame(CHR = 5, P0 = 44.5e6, P1 = 50.5e6),
    data.frame(CHR = 5, P0 = 98e6, P1 = 100.5e6),
    data.frame(CHR = 5, P0 = 129e6, P1 = 132e6),
    data.frame(CHR = 5, P0 = 135.5e6, P1 = 138.5e6),
    data.frame(CHR = 6, P0 = 25.5e6, P1 = 33.5e6),
    data.frame(CHR = 6, P0 = 57e6, P1 = 64e6),
    data.frame(CHR = 6, P0 = 140e6, P1 = 142.5e6),
    data.frame(CHR = 7, P0 = 55e6, P1 = 66e6),
    data.frame(CHR = 8, P0 = 8e6, P1 = 12e6),
    data.frame(CHR = 8, P0 = 43e6, P1 = 50e6),
    data.frame(CHR = 8, P0 = 112e6, P1 = 115e6),
    data.frame(CHR = 10, P0 = 37e6, P1 = 43e6),
    data.frame(CHR = 11, P0 = 46e6, P1 = 57e6),
    data.frame(CHR = 11, P0 = 87.5e6, P1 = 90.5e6),
    data.frame(CHR = 12, P0 = 33e6, P1 = 40e6),
    data.frame(CHR = 12, P0 = 109.5e6, P1 = 112e6),
    data.frame(CHR = 20, P0 = 32e6, P1 = 34.5e6)
))

# Remove variants within genomic regions (REF: PMC2443852)
remove_regions<-function(bim, regions){
  exclude<-NULL
  for(i in 1:nrow(regions)){
    exclude<-c(exclude, bim$SNP[(   bim$CHR == regions$CHR[i]  &
                                    bim$BP >= regions$P0[i] &
                                    bim$BP <= regions$P1[i])])
  }

  return(bim[!(bim$SNP %in% exclude),])
}

# Merge plink files
plink_merge<-function(bfile, plink = 'plink', chr = 1:22, extract = NULL, keep = NULL, flip = NULL, memory = 4000, out){
  tmp_dir<-tempdir()

  # Create merge list
  ref_merge_list<-paste0(bfile, chr)
  write.table(ref_merge_list, paste0(tmp_dir,'/ref_mergelist.txt'), row.names=F, col.names=F, quote=F)

  # Prepare command
  plink_opt<-NULL
  if(!is.null(extract)){
    write.table(extract, paste0(tmp_dir,'/extract.snplist'), col.names = F, row.names = F, quote=F)
    plink_opt<-paste0(plink_opt, paste0('--extract ',tmp_dir,'/extract.snplist '))
  }
  if(!is.null(flip)){
    write.table(flip, paste0(tmp_dir,'/flip_list.txt'), col.names = F, row.names = F, quote=F)
    plink_opt<-paste0(plink_opt, paste0('--flip ',tmp_dir,'/flip_list.txt '))
  }
  if(!is.null(keep)){
    # Check extract file
    keep<-obj_or_file(keep)
    plink_opt<-paste0(plink_opt, paste0('--keep ', keep, ' '))
  }

  if(length(chr) > 1){
    # Merge
    cmd<-paste0(plink,' --merge-list ', tmp_dir,'/ref_mergelist.txt ', plink_opt,'--threads 1 --make-bed --out ', out,' --memory ', memory)
  } else {
    # Skip merge
    cmd<-paste0(plink,' --bfile ', bfile, chr, ' ', plink_opt,'--threads 1 --make-bed --out ', out,' --memory ', memory)
  }
  system(cmd)
}

# Perform PCA using plink files
plink_pca<-function(bfile, chr = 1:22, plink = 'plink', plink2 = 'plink2', extract = NULL, keep = NULL, flip = NULL, memory = 4000, n_pc = 6){
  tmp_dir<-tempdir()

  # Merge subset reference
  plink_merge(bfile = bfile, chr = chr, plink = plink, keep = keep, extract = extract, flip = flip, memory = memory, out = paste0(tmp_dir,'/ref_merge'))

  # Calculate SNP weights
  system(paste0(plink2,' --bfile ',tmp_dir,'/ref_merge --threads 1 --pca ',n_pc,' biallelic-var-wts  --out ',tmp_dir,'/ref_merge --memory ', memory))

  # Format the SNP-weights
  snp_weights<-fread(paste0(tmp_dir,'/ref_merge.eigenvec.var'))
  snp_weights<-snp_weights[, -1, with=F]
  names(snp_weights)[1:3]<-c('SNP','A1','A2')

  return(snp_weights)
}

# Performing LD pruning
plink_prune<-function(bfile, keep = NULL, plink = 'plink', chr = 1:22, extract = NULL, memory = 4000){
  # Create a temporary file path to store pruning output
  tmp_file<-tempfile()

  # Check extract file
  extract<-obj_or_file(extract)

  # Prepare plink options
  plink_opt<-NULL
  if(!is.null(extract)){
    extract <- obj_or_file(extract)
    plink_opt<-paste0(plink_opt, '--extract ', extract,' ')
  }
  if(!is.null(keep)){
    keep <- obj_or_file(keep)
    plink_opt<-paste0(plink_opt, '--keep ', keep,' ')
  }

  # Perfom pruning and read in SNP-list
  ld_indep<-NULL
  for(chr_i in chr){
    cmd <- paste0(plink,' --bfile ',bfile,chr_i,' --threads 1 ',plink_opt,'--indep-pairwise 1000 5 0.2 --out ',tmp_file,'.chr',chr_i,' --memory ',memory)
    exit_status <- system(cmd, intern=FALSE)
    if(file.exists(paste0(tmp_file,'.chr',chr_i,'.prune.in'))){
      ld_indep<-c(ld_indep, fread(paste0(tmp_file,'.chr',chr_i,'.prune.in'), header=F)$V1)
    }
    if (exit_status == 2) {
      stop()
    }
  }

  return(ld_indep)
}

# Set range of chromosome numbers to use by default
CHROMS <- 1:22

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

# Peforming LD-based clumping
plink_clump<-function(bfile, plink = 'plink', chr = 1:22, sumstats, keep = NULL, memory = 4000, log_file = NULL){
  log_add(log_file = log_file, message = 'Performing LD-based clumping.')
  tmp_file <- tempfile()

  sumstats <- obj_or_file(sumstats, header=T)

  opt_plink<-NULL
  if(!is.null(keep)){
    keep <- obj_or_file(keep)
    opt_plink<-paste0(opt_plink, '--keep ', keep, ' ')
  }

  clumped<-NULL
  for(chr_i in chr){
    cmd<-paste0(plink,' --bfile ', bfile, chr_i,' ',opt_plink, '--clump ', sumstats,' --clump-p1 1 --clump-p2 1 --clump-r2 0.1 --clump-kb 250 --out ',tmp_file,'.chr',chr_i,' --memory ', memory)
    exit_status <- system(cmd, intern=FALSE)
    if (file.exists(paste0(tmp_file,'.chr',chr_i,'.clumped'))) {
      clumped <- c(clumped, fread(paste0(tmp_file,'.chr',chr_i,'.clumped'))$SNP)
    }
    if(exit_status == 2){
      stop()
    }
  }

  log_add(log_file = log_file, message = paste0(length(clumped),' variants remain after clumping.'))

  return(clumped)
}

test_start<-function(log_file){
  test_start.time <- Sys.time()
  log_add(log_file = log_file, message = paste0('Test started at ',as.character(test_start.time)))
  return(test_start.time)
}

test_finish<-function(log_file, test_start.time){
  end.time <- Sys.time()
  time.taken <- end.time - test_start.time
  sink(file = log_file, append = T)
  cat('Test run finished at',as.character(end.time),'\n')
  cat('Test duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
  sink()
}

# Run LDSC
ldsc <- function(sumstats, ldsc, munge_sumstats, ldsc_ref, pop_prev = NULL, sample_prev = NULL, log_file = NULL){
  tmp_dir<-tempdir()

  # Munge the sumstats
  system(paste0(munge_sumstats, ' --sumstats ', sumstats,' --merge-alleles ', ldsc_ref, '/w_hm3.snplist --out ', tmp_dir,'/munged'))

  # Define the file paths
  sumstats_path <- file.path(tmp_dir,'/munged.sumstats.gz')
  output_path <- file.path(tmp_dir, '/ldsc')

  # Run the appropriate LDSC command based on the availability of prevalence data
  if(!is.null(pop_prev) && !is.null(sample_prev)) {
    system(paste0(ldsc, ' --h2 ', sumstats_path, ' --ref-ld-chr ', ldsc_ref, '/ --w-ld-chr ', ldsc_ref, '/ --out ', output_path, ' --samp-prev ', sample_prev, ' --pop-prev ', pop_prev))
    scale_type <- "Liability"
  } else {
    system(paste0(ldsc, ' --h2 ', sumstats_path, ' --ref-ld-chr ', ldsc_ref, '/ --w-ld-chr ', ldsc_ref, '/ --out ', output_path))
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
dbslmm <- function(dbslmm, plink = plink, ld_blocks, chr = 1:22, bfile, sumstats, h2, nsnp, nindiv, log_file = NULL){
  # Create temp directory
  tmp_dir<-tempdir()

  # Make sure there is permission to run dbslmm
  system(paste0('chmod 777 ', dbslmm))

  # Run dbslmm for each chromosome
  for(chr_i in chr){
    cmd <-paste0('Rscript ', dbslmm,'/DBSLMM.R --plink ', plink,' --block ', ld_blocks,'/fourier_ls-chr', chr_i, '.bed --dbslmm ', dbslmm, '/dbslmm --h2 ', ldsc_h2,' --ref ', bfile, chr_i, ' --summary ', sumstats, chr_i, '.assoc.txt --n ', nindiv,' --nsnp ', nsnp, ' --outPath ', tmp_dir, '/ --thread 1')
    exit_status <- system(cmd, intern = F)
  }

  # Read in DBSLMM output
  dbslmm_output <- list.files(path=tmp_dir, pattern='.dbslmm.txt')
  if(length(dbslmm_output) != 22){
    log_add(log_file = log_file, message = 'At least one chromosome did not complete.')
  }

  dbslmm_all<-NULL
  for(i in dbslmm_output){
    dbslmm_all<-rbind(dbslmm_all, fread(paste0(tmp_dir,'/',i)))
  }

  dbslmm_all<-dbslmm_all[,c(1,2,4), with=T]
  names(dbslmm_all)<-c('SNP', 'A1', 'SCORE_DBSLMM')

  # Insert A2 information
  bim<-read_bim(bfile, chr = chr)
  dbslmm_all<-merge(dbslmm_all, bim[, c('SNP','A1','A2'), with = F], by = c('SNP','A1'), all.x = T)
  if(any(is.na(dbslmm_all$A2))){
    stop('Insertion of A2 in dbslmm score file failed.')
  }

  dbslmm_all <- dbslmm_all[, c('SNP','A1','A2','SCORE_DBSLMM')]

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

# Generate kinship matrix and identify unrelated individuals (>2nd degree)
plink_king<-function(bfile, extract = NULL, chr = 1:22, plink = 'plink', plink2 = 'plink2', out, keep=NA){
  # Create object indicating tmpdir
  tmp_dir<-tempdir()

  # Identify variants not in high LD regions
  bim <- read_bim(dat = bfile, chr = chr)
  bim <- remove_regions(bim = bim, regions = long_ld_coord)

  # Identify intersect of snplists to extract
  if(!is.null(extract)){
    extract_snplist<-intersect(bim$SNP, extract)
  }

  # Merge per chromosome files extracting selected variants
  plink_merge(bfile = bfile, chr = chr, plink = plink, extract = extract_snplist, out = paste0(tmp_dir,'/merged'))

  # Run KING estimator
  system(paste0(plink2, ' --bfile ', tmp_dir, '/merged --threads 1 --make-king triangle bin --out ', tmp_dir, '/merged'))

  # Identify unrelated individuals (remove 2nd degree relatives)
  system(paste0(plink2, ' --bfile ', tmp_dir, '/merged --threads 1 --king-cutoff ',tmp_dir, '/merged 0.0884 --out ', tmp_dir, '/merged'))

  # Move the kinship matrix
  system(paste0('mv ', tmp_dir, '/merged.king.bin ', out, '.king.bin'))
  system(paste0('mv ', tmp_dir, '/merged.king.id ', out, '.king.id'))

  # Format and save the list of unrelated individuals
  system(paste0('tail -n +2 ', tmp_dir, '/merged.king.cutoff.in.id > ', out, '.unrelated.keep'))
}

# Read in PGS
read_pgs <- function(config, name = NULL, pgs_methods = NULL, gwas = NULL, pop = NULL){

  # Read in target_list
  target_list <- read_param(config = config, param = 'target_list')
  if(!is.null(name)){
    if(any(!(name %in% target_list$name))){
      stop('Requested target samples are not present in target_list')
    }
    target_list <- target_list[target_list$name %in% name,]
  }

  # Read in gwas_list
  gwas_list <- read_param(config = config, param = 'gwas_list')

  # Read in score_list
  score_list <- read_param(config = config, param = 'score_list')

  if(!is.null(gwas)){
    if(!is.null(score_list)){
      full_gwas_list <- c(gwas_list$name, score_list$name)
    } else {
      full_gwas_list <- gwas_list$name
    }

    if(any(!(gwas %in% full_gwas_list))){
      stop('Requested GWAS are not present in gwas_list/score_list')
    }
    gwas_list <- gwas_list[target_list$name %in% gwas,]

    if(!is.null(score_list)){
      score_list <- score_list[score_list$name %in% gwas,]
    }
  }

  # Identify PGS methods to be included
  pgs_methods_list <- read_param(config = config, param = 'pgs_methods', return_obj = F)

  if(!is.null(pgs_methods)){
    if(!is.null(score_list)){
      if(any(!(pgs_methods %in% c(pgs_methods_list, 'external')))){
        stop('Requested pgs_methods are not present in pgs_methods in config')
      }
    } else {
      if(any(!(pgs_methods %in% pgs_methods_list))){
        stop('Requested pgs_methods are not present in pgs_methods in config')
      }
    }
    pgs_methods_list <- pgs_methods_list[pgs_methods_list %in% pgs_methods]
  }

  # Define PGS methods applied to non-EUR GWAS
  pgs_methods_noneur <- c('ptclump','lassosum','megaprs')

  # Identify outdir parameter
  outdir <- read_param(config = config, param = 'outdir', return_obj = F)

  pgs <- list()
  for (name_i in target_list$name) {
    # Read in keep_list to determine populations available
    keep_list_i <- fread(paste0(outdir,'/',name_i,'/ancestry/keep_list.txt'))
    pgs[[name_i]] <- list()
    for (pop_i in keep_list_i$POP) {
      pgs[[name_i]][[pop_i]] <- list()
      for (gwas_i in gwas_list$name) {
        pgs[[name_i]][[pop_i]][[gwas_i]] <- list()

        for (pgs_method_i in pgs_methods_list) {
          if (gwas_list$population[gwas_list$name == gwas_i] == 'EUR' | (gwas_list$population[gwas_list$name == gwas_i] != 'EUR' & (pgs_method_i %in% pgs_methods_noneur))) {
            pgs[[name_i]][[pop_i]][[gwas_i]][[pgs_method_i]] <-
              fread(
                paste0(
                  outdir, '/', name_i, '/pgs/', pop_i, '/', pgs_method_i, '/',  gwas_i, '/', name_i, '-', gwas_i, '-', pop_i, '.profiles'
                )
              )
          }
        }
      }
      if(!is.null(score_list)){
        for (score_i in score_list$name) {
          pgs[[name_i]][[pop_i]][[score_i]] <- list()
          pgs_method_i <- 'external'
          pgs[[name_i]][[pop_i]][[score_i]][[pgs_method_i]] <-
            fread(
              paste0(
                outdir, '/', name_i, '/pgs/', pop_i, '/', pgs_method_i, '/',  score_i, '/', name_i, '-', score_i, '-', pop_i, '.profiles'
              )
            )
        }
      }
    }
  }

  return(pgs)
}

# Create function to read in parameters in the config file
read_param <- function(config, param, return_obj = T){

  # Read in the config file
  config_file <- readLines(config)

  if(all(grepl(paste0('^',param,':'), config_file) == F)){
    # Check default config file
    config_file <- readLines(config)

    if(all(grepl(paste0('^',param,':'), config_file) == F)){
      cat('Requested parameter is not present in user specified config file or default config file.')
      return(NULL)
    } else {
      cat('Parameter is not present in user specified config file, so will use value in default config file.')
    }
  }

  # Identify value for param
  file <- gsub(paste0(param,': '), '', config_file[grepl(paste0('^',param,':'), config_file)])
  file[file == 'NA']<-NA

  if(return_obj){
    if(!is.na(file)){
      obj <- fread(file)
    } else {
      obj <- NULL
    }
    return(obj)
  } else {

    file <- unlist(strsplit(gsub("'", '', gsub(']', '', gsub('\\[', '', file))),','))
    file <- file[order(file)]

    return(file)
  }
}

# Read in ancestry classifications
read_ancestry <- function(config, name){

  # Read in the config file
  config_file <- readLines(config)

  # Identify outdir
  outdir <- read_param(config = params$config, param = 'outdir', return_obj = F)

  # Read keep_list
  keep_list <- fread(paste0(outdir,'/',name,'/ancestry/keep_list.txt'))

  # Read in keep lists
  keep_files <- list()
  for(pop in keep_list$POP){
    keep_files[[pop]] <- fread(keep_list$file[keep_list$POP == pop], header = F)
  }

  # Read in model predictions
  model_pred <- fread(paste0(outdir,'/',name,'/ancestry/',name,'.Ancestry.model_pred'))

  # Read in ancestry log file
  log <- readLines(paste0(outdir,'/',name,'/ancestry/',name,'.Ancestry.log'))

  output <- list(
    keep_list = keep_list,
    keep_files = keep_files,
    model_pred = model_pred,
    log = log
  )

  return(output)
}

# Return score corresponding to pseudovalidation
find_pseudo <- function(config, gwas, pgs_method){

  if(length(pgs_method) > 1){
    stop('Only one pgs_method can be specified at a time')
  }
  if(length(gwas) > 1){
    stop('Only one gwas can be specified at a time')
  }

  # Read in gwas_list
  gwas_list <- read_param(config = config, param = 'gwas_list')

  # Read in score_list
  score_list <- read_param(config = config, param = 'score_list')

  if(!is.null(gwas)){
    if(!is.null(score_list)){
      full_gwas_list <- c(gwas_list$name, score_list$name)
    } else {
      full_gwas_list <- gwas_list$name
    }

    if(any(!(gwas %in% full_gwas_list))){
      stop('Requested GWAS are not present in gwas_list/score_list')
    }
    gwas_list <- gwas_list[target_list$name %in% gwas,]

    if(!is.null(score_list)){
      score_list <- score_list[score_list$name %in% gwas,]
    }
  }

  # Identify PGS methods to be included
  pgs_methods_list <- read_param(config = config, param = 'pgs_methods', return_obj = F)

  if(!is.null(pgs_method)){
    if(!is.null(score_list)){
      if(any(!(pgs_method %in% c(pgs_methods_list, 'external')))){
        stop('Requested pgs_method are not present in pgs_methods in config')
      }
    } else {
      if(any(!(pgs_method %in% pgs_methods_list))){
        stop('Requested pgs_method are not present in pgs_methods in config')
      }
    }
    pgs_methods_list <- pgs_methods_list[pgs_methods_list %in% pgs_method]
  }

  # Identify outdir parameter
  outdir <- read_param(config = config, param = 'outdir', return_obj = F)

  # Use most stringent p-value threshold of 0.05 as pseudo
  if(pgs_method == 'ptclump'){
    return('0_1')
  }

  # Pseudoval only methods
  if(pgs_method == 'sbayesr'){
    return('SBayesR')
  }
  if(pgs_method == 'dbslmm'){
    return('DBSLMM')
  }

  # Retrieve pseudoval param
  if(pgs_method == 'ldpred2'){
    return('beta_auto')
  }
  if(pgs_method == 'prscs'){
    return('phi_auto')
  }
  if(pgs_method == 'megaprs'){
    # Read in megaprs log file
    log <- readLines(paste0(outdir,'/resources/data/ref/pgs_score_files/',pgs_method,'/',gwas,'/ref-',gwas,'.log'))
    log <- log[grepl('identified as the best with correlation', log)]
    pseudoval <- gsub(' .*','', gsub('Model ', '', log))
    return(paste0('ldak_Model', pseudoval))
  }
  if(pgs_method == 'lassosum'){
    # Read in megaprs log file
    log <- readLines(paste0(outdir,'/resources/data/ref/pgs_score_files/',pgs_method,'/',gwas,'/ref-',gwas,'.log'))
    s_val <- gsub('.* ', '', log[grepl('^s = ', log)])
    lambda_val <- gsub('.* ', '', log[grepl('^lambda = ', log)])
    return(paste0('s', s_val, '_lambda', lambda_val))
  }

  # If pgs_method is external, return the only score
  if(pgs_method == 'external'){
    return('external')
  }
}

# Read in lassosum pseudoval results
read_pseudo_r <- function(config, gwas){

  if(length(gwas) > 1){
    stop('Only one gwas can be specified at a time')
  }

  # Find outdir param
  outdir <- read_param(config = config, param = 'outdir', return_obj = F)

  # Read in lassosum log file
  log <- readLines(paste0(outdir,'/resources/data/ref/pgs_score_files/lassosum/',gwas,'/ref-',gwas,'.log'))
  r <- as.numeric(gsub('value = ','',log[grepl('value = ', log)]))

  return(r)
}

# Read in external score file
read_score <- function(score, log_file = NULL){
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

	return(score)
}
