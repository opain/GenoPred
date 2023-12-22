#!/usr/bin/Rscript

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
read_geno<-function(target, format){

  if(format == 'samp_imp_plink1'){
    target_snp<-fread(paste0(target,'.bim'))
    target_snp$V3<-NULL
    names(target_snp)<-c('CHR','SNP','BP','A1','A2')
  }

  if(format == 'samp_imp_bgen'){
    connection = dbConnect( RSQLite::SQLite(), paste0(target,'.bgen.bgi'))
    target_snp = dbGetQuery( connection, "SELECT * FROM Variant" )
    target_snp<-target_snp[,c('chromosome','rsid','position','allele1','allele2')]
    names(target_snp)<-c('CHR','SNP','BP','A1','A2')
    dbDisconnect(connection)
    target_snp<-data.table(target_snp)
    target_snp$CHR<-as.numeric(gsub('chr','',target_snp$CHR))
  }

  if(format == 'samp_imp_vcf'){
    target_snp<-fread(cmd=paste0("zcat ",target,".vcf.gz | cut -f 1-5"))
    names(target_snp)<-c('CHR','BP','SNP','A1','A2')
    target_snp$CHR<-as.numeric(gsub('chr','',target_snp$CHR))
  }

  target_snp<-target_snp[, c('CHR','BP','SNP','A1','A2'), with=F]

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
plink_subset<-function(plink='plink', chr = 1:22, keep = NA, bfile, out, memory = 4000, threads = 1){

  # If object, create file
  keep <- obj_or_file(keep)

  # Prepare plink options
  plink_opt<-NULL
  if(!is.null(keep)){
    plink_opt<-paste0(plink_opt, paste0('--keep ',keep,' '))
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
plink_qc_snplist<-function(bfile, plink = 'plink', chr = 1:22, threads = 1, memory = 4000, geno = NULL, maf = NULL, hwe = NULL){
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

  temp_file<-tempfile()
  snplist<-NULL
  for(chr_i in chr){
    system(paste0(plink,' --bfile ',bfile,chr_i,' --threads ',threads,' ',plink_opt,' --write-snplist --out ', temp_file,' --memory ', memory))
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

# Perform PCA using plink files
plink_pca<-function(bfile, plink = 'plink', plink2 = 'plink2', extract = NULL, flip = NULL, memory = 4000, n_pc = 6){
  ###########
  # Merge subset reference
  ###########
  tmp_dir<-tempdir()

  # Create merge list
  ref_merge_list<-paste0(bfile,1:22)
  write.table(ref_merge_list, paste0(tmp_dir,'/ref_mergelist.txt'), row.names=F, col.names=F, quote=F)

  # Merge
  plink_opt<-NULL
  if(!is.null(extract)){
    write.table(extract, paste0(tmp_dir,'/extract.snplist'), col.names = F, row.names = F, quote=F)
    plink_opt<-c(plink_opt, paste0('--extract ',tmp_dir,'/extract.snplist '))
  }
  if(!is.null(flip)){
    write.table(flip, paste0(tmp_dir,'/flip_list.txt'), col.names = F, row.names = F, quote=F)
    plink_opt<-c(plink_opt, paste0('--flip ',tmp_dir,'/flip_list.txt '))
  }
  cmd<-paste0(plink,' --merge-list ',tmp_dir,'/ref_mergelist.txt ',plink_opt,'--threads 1 --make-bed --out ',tmp_dir,'/ref_merge --memory ',memory)
  system(cmd)

  # Calculate SNP weights
  system(paste0(plink2,' --bfile ',tmp_dir,'/ref_merge --threads 1 --pca ',n_pc,' biallelic-var-wts  --out ',tmp_dir,'/ref_merge --memory ', memory))

  # Format the SNP-weights
  snp_weights<-fread(paste0(tmp_dir,'/ref_merge.eigenvec.var'))
  snp_weights<-snp_weights[, -1, with=F]
  names(snp_weights)[1:3]<-c('SNP','A1','A2')

  return(snp_weights)
}

# Performing LD pruning
plink_prune<-function(bfile, plink = 'plink', chr =1:22, extract = NULL, memory = 4000){
  # Create a temporary file path to store pruning output
  tmp_file<-tempfile()

  # Check extract file
  extract<-obj_or_file(extract)

  # Prepare plink options
  plink_opt<-NULL
  if(!is.null(extract)){
    plink_opt<-paste0(plink_opt, '--extract ', extract,' ')
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
  if(all(!(chr %in% 1:22))){
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
  if(!(all(c('CHR','BP','A1','A2','BETA','SE','P','FREQ') %in% names(gwas)))){
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
  swap_index <- sumstats_ref$A1.x == sumstats_ref$A1.y & sumstats_ref$A2.x == sumstats_ref$A2.y
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
