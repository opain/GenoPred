#!/usr/bin/Rscript

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

# Calculate scores in target file
plink_score<-function(bfile, score, keep=NULL, extract=NULL, chr=1:22, frq=NULL, plink2='plink2', threads=1){
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
