#!/usr/bin/Rscript

# Make a subset of plink1 binaries
plink_subset<-function(plink=NULL, plink2=NULL, chr = 1:22, keep = NULL, extract = NULL, bfile=NULL, pfile=NULL, out, memory = 4000, threads = 1){
  if(is.null(bfile) & is.null(pfile)){
    stop("bfile or pfile must be specified.")
  }
  if(!is.null(bfile) & !is.null(pfile)){
    stop("Both bfile and pfile cannot be specified.")
  }
  if(is.null(plink) & is.null(plink2)){
    stop("plink or plink2 must be specified.")
  }
  if(!is.null(plink) & !is.null(plink2)){
    stop("Both plink and plink2 cannot be specified.")
  }
  if(!is.null(pfile) & is.null(plink2)){
    stop("plink2 must be specified when using pfile.")
  }

  # Prepare plink options
  plink_opt<-NULL
  if(!is.null(plink)){
    plink_opt<-paste0(plink_opt, paste0(plink, ' '))
  } else {
    plink_opt<-paste0(plink_opt, paste0(plink2, ' '))
  }
  if(!is.null(bfile)){
    plink_opt<-paste0(plink_opt, paste0('--bfile ',bfile,'CHROMOSOME_NUMBER --make-bed '))
  } else {
    plink_opt<-paste0(plink_opt, paste0('--pfile ',pfile,'CHROMOSOME_NUMBER --make-pgen '))
  }
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
    plink_opt <- gsub('CHROMOSOME_NUMBER', chr_i, plink_opt)
    cmd <- paste0(plink_opt, '--threads ', threads,' --out ',out,chr_i,' --memory ',memory)
    exit_status <- system(cmd, intern=FALSE)
    if (exit_status == 2) {
      stop()
    }
  }
}

# Return a list of variants surviving QC using plink
plink_qc_snplist<-function(bfile=NULL, pfile=NULL, plink=NULL, plink2=NULL, keep = NULL, chr = 1:22, threads = 1, memory = 4000, geno = NULL, maf = NULL, hwe = NULL){
  if(is.null(bfile) & is.null(pfile)){
    stop("bfile or pfile must be specified.")
  }
  if(!is.null(bfile) & !is.null(pfile)){
    stop("Both bfile and pfile cannot be specified.")
  }
  if(is.null(plink) & is.null(plink2)){
    stop("plink or plink2 must be specified.")
  }
  if(!is.null(plink) & !is.null(plink2)){
    stop("Both plink and plink2 cannot be specified.")
  }
  if(!is.null(pfile) & is.null(plink2)){
    stop("plink2 must be specified when using pfile.")
  }

  plink_opt<-NULL
  if(!is.null(plink)){
    plink_opt<-paste0(plink_opt, paste0(plink, ' '))
  } else {
    plink_opt<-paste0(plink_opt, paste0(plink2, ' '))
  }
  if(!is.null(bfile)){
    plink_opt<-paste0(plink_opt, paste0('--bfile ',bfile,'CHROMOSOME_NUMBER '))
  } else {
    plink_opt<-paste0(plink_opt, paste0('--pfile ',pfile,'CHROMOSOME_NUMBER '))
  }
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
    plink_opt <- gsub('CHROMOSOME_NUMBER', chr_i, plink_opt)
    cmd <- paste0(plink_opt,'--threads ', threads,' --write-snplist --out ', temp_file,' --memory ', memory)
    exit_status <- system(cmd, intern=FALSE)
    if (exit_status == 2) {
      stop()
    }
    snplist<-c(snplist, fread(paste0(temp_file, '.snplist'), header=F)$V1)
  }

  return(snplist)
}

# Merge plink files
plink_merge<-function(bfile=NULL, pfile=NULL, plink=NULL, plink2=NULL, chr = 1:22, extract = NULL, keep = NULL, flip = NULL, memory = 4000, out){
  if(is.null(bfile) & is.null(pfile)){
    stop("bfile or pfile must be specified.")
  }
  if(!is.null(bfile) & !is.null(pfile)){
    stop("Both bfile and pfile cannot be specified.")
  }
  if(is.null(plink) & is.null(plink2)){
    stop("plink or plink2 must be specified.")
  }
  if(!is.null(plink) & !is.null(plink2)){
    stop("Both plink and plink2 cannot be specified.")
  }
  if(!is.null(pfile) & is.null(plink2)){
    stop("plink2 must be specified when using pfile.")
  }

  tmp_dir<-tempdir()

  # Create merge list
  if(!is.null(bfile)){
    ref_merge_list<-paste0(bfile, chr)
  } else {
    ref_merge_list<-paste0(pfile, chr)
  }
  write.table(ref_merge_list, paste0(tmp_dir,'/ref_mergelist.txt'), row.names=F, col.names=F, quote=F)

  # Prepare command
  plink_opt<-NULL
  if(!is.null(plink)){
    plink_opt<-paste0(plink_opt, paste0(plink, ' '))
  } else {
    plink_opt<-paste0(plink_opt, paste0(plink2, ' '))
  }
  if(!is.null(bfile)){
    plink_opt<-paste0(plink_opt, paste0('--make-bed '))
  } else {
    plink_opt<-paste0(plink_opt, paste0('--make-pgen '))
  }
  if(length(chr) > 1){
    if(!is.null(bfile)){
      plink_opt<-paste0(plink_opt, paste0('--merge-list ', tmp_dir,'/ref_mergelist.txt '))
    } else {
      plink_opt<-paste0(plink_opt, paste0('--pmerge-list ', tmp_dir,'/ref_mergelist.txt '))
    }
  } else {
    if(!is.null(bfile)){
      plink_opt<-paste0(plink_opt, paste0('--bfile ', bfile, chr, ' '))
    } else {
      plink_opt<-paste0(plink_opt, paste0('--pfile ', pfile, chr, ' '))
    }
  }
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

  cmd<-paste0(plink_opt,'--threads 1 --out ', out,' --memory ', memory)
  exit_status <- system(cmd, intern=FALSE)
  if (exit_status == 2) {
    stop()
  }
}

# Perform PCA using plink files
plink_pca<-function(bfile=NULL, pfile=NULL, chr = 1:22, plink2, extract = NULL, keep = NULL, flip = NULL, memory = 4000, n_pc = 6){
  if(is.null(bfile) & is.null(pfile)){
    stop("bfile or pfile must be specified.")
  }
  if(!is.null(bfile) & !is.null(pfile)){
    stop("Both bfile and pfile cannot be specified.")
  }

  tmp_dir<-tempdir()

  # Merge subset reference
  if(!is.null(bfile)){
    plink_merge(bfile = bfile, chr = chr, plink2 = plink2, keep = keep, extract = extract, flip = flip, memory = memory, out = paste0(tmp_dir,'/ref_merge'))
  } else {
    plink_merge(pfile = pfile, chr = chr, plink2 = plink2, keep = keep, extract = extract, flip = flip, memory = memory, out = paste0(tmp_dir,'/ref_merge'))
  }

  plink_opt<-paste0(plink2, ' ')
  if(!is.null(bfile)){
    plink_opt<-paste0(plink_opt, paste0('--bfile ',tmp_dir,'/ref_merge '))
  } else {
    plink_opt<-paste0(plink_opt, paste0('--pfile ',tmp_dir,'/ref_merge '))
  }

  # Calculate SNP weights
  system(paste0(plink_opt,' --threads 1 --pca ',n_pc,' biallelic-var-wts  --out ',tmp_dir,'/ref_merge --memory ', memory))

  # Format the SNP-weights
  snp_weights<-fread(paste0(tmp_dir,'/ref_merge.eigenvec.var'))
  snp_weights<-snp_weights[, -1, with=F]
  names(snp_weights)[1:3]<-c('SNP','A1','A2')

  return(snp_weights)
}

# Performing LD pruning
plink_prune<-function(bfile=NULL, pfile=NULL, keep = NULL, plink=NULL, plink2=NULL, chr = 1:22, extract = NULL, memory = 4000){
  if(is.null(bfile) & is.null(pfile)){
    stop("bfile or pfile must be specified.")
  }
  if(!is.null(bfile) & !is.null(pfile)){
    stop("Both bfile and pfile cannot be specified.")
  }
  if(is.null(plink) & is.null(plink2)){
    stop("plink or plink2 must be specified.")
  }
  if(!is.null(plink) & !is.null(plink2)){
    stop("Both plink and plink2 cannot be specified.")
  }
  if(!is.null(pfile) & is.null(plink2)){
    stop("plink2 must be specified when using pfile.")
  }

  # Create a temporary file path to store pruning output
  tmp_file<-tempfile()

  # Check extract file
  extract<-obj_or_file(extract)

  # Prepare plink options
  plink_opt<-NULL
  if(!is.null(plink)){
    plink_opt<-paste0(plink_opt, paste0(plink, ' '))
  } else {
    plink_opt<-paste0(plink_opt, paste0(plink2, ' '))
  }
  if(!is.null(bfile)){
    plink_opt<-paste0(plink_opt, paste0('--bfile ',bfile,'CHROMOSOME_NUMBER '))
  } else {
    plink_opt<-paste0(plink_opt, paste0('--pfile ',pfile,'CHROMOSOME_NUMBER '))
  }
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
    plink_opt <- gsub('CHROMOSOME_NUMBER', chr_i, plink_opt)
    cmd <- paste0(plink_opt, '--threads 1 --indep-pairwise 1000 5 0.2 --out ',tmp_file,'.chr',chr_i,' --memory ',memory)
    exit_status <- system(cmd, intern=FALSE)
    if (exit_status == 2) {
      stop()
    }
    if(file.exists(paste0(tmp_file,'.chr',chr_i,'.prune.in'))){
      ld_indep<-c(ld_indep, fread(paste0(tmp_file,'.chr',chr_i,'.prune.in'), header=F)$V1)
    }
  }
  return(ld_indep)
}

# Peforming LD-based clumping
plink_clump<-function(bfile=NULL, pfile=NULL, plink=NULL, plink2=NULL, chr = 1:22, sumstats, keep = NULL, memory = 4000, log_file = NULL){
  if(is.null(bfile) & is.null(pfile)){
    stop("bfile or pfile must be specified.")
  }
  if(!is.null(bfile) & !is.null(pfile)){
    stop("Both bfile and pfile cannot be specified.")
  }
  if(is.null(plink) & is.null(plink2)){
    stop("plink or plink2 must be specified.")
  }
  if(!is.null(plink) & !is.null(plink2)){
    stop("Both plink and plink2 cannot be specified.")
  }
  if(!is.null(pfile) & is.null(plink2)){
    stop("plink2 must be specified when using pfile.")
  }

  log_add(log_file = log_file, message = 'Performing LD-based clumping.')
  tmp_file <- tempfile()

  sumstats <- obj_or_file(sumstats, header=T)

  plink_opt<-NULL
  if(!is.null(plink)){
    plink_opt<-paste0(plink_opt, paste0(plink, ' '))
  } else {
    plink_opt<-paste0(plink_opt, paste0(plink2, ' '))
  }
  if(!is.null(bfile)){
    plink_opt<-paste0(plink_opt, paste0('--bfile ',bfile,'CHROMOSOME_NUMBER '))
  } else {
    plink_opt<-paste0(plink_opt, paste0('--pfile ',pfile,'CHROMOSOME_NUMBER '))
  }
  if(!is.null(keep)){
    keep <- obj_or_file(keep)
    opt_plink<-paste0(opt_plink, '--keep ', keep, ' ')
  }

  clumped<-NULL
  for(chr_i in chr){
    cmd<-paste0(plink,' --bfile ', bfile, chr_i,' ',opt_plink, '--clump ', sumstats,' --clump-p1 1 --clump-p2 1 --clump-r2 0.1 --clump-kb 250 --out ',tmp_file,'.chr',chr_i,' --threads 1 --memory ', memory)
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
plink_king<-function(bfile=NULL, pfile=NULL, extract = NULL, chr = 1:22, plink=NULL, plink2=NULL, out, keep=NA){
  if(is.null(bfile) & is.null(pfile)){
    stop("bfile or pfile must be specified.")
  }
  if(!is.null(bfile) & !is.null(pfile)){
    stop("Both bfile and pfile cannot be specified.")
  }
  if(is.null(plink) & is.null(plink2)){
    stop("plink or plink2 must be specified.")
  }
  if(!is.null(plink) & !is.null(plink2)){
    stop("Both plink and plink2 cannot be specified.")
  }
  if(!is.null(pfile) & is.null(plink2)){
    stop("plink2 must be specified when using pfile.")
  }

  # Create object indicating tmpdir
  tmp_dir<-tempdir()

  # Identify variants not in high LD regions
  bim <- read_bim(dat = bfile, chr = chr)
  bim <- remove_regions(bim = bim, regions = long_ld_coord)

  # Identify intersect of snplists to extract
  if(!is.null(extract)){
    extract_snplist<-intersect(bim$SNP, extract)
  }

  plink_opt<-NULL
  if(!is.null(plink)){
    plink_opt<-paste0(plink_opt, paste0(plink, ' '))
  } else {
    plink_opt<-paste0(plink_opt, paste0(plink2, ' '))
  }
  if(!is.null(bfile)){
    plink_opt<-paste0(plink_opt, paste0('--bfile ',bfile,'CHROMOSOME_NUMBER '))
  } else {
    plink_opt<-paste0(plink_opt, paste0('--pfile ',pfile,'CHROMOSOME_NUMBER '))
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
plink_score<-function(bfile=NULL, pfile=NULL, score, keep=NULL, extract=NULL, chr=1:22, frq=NULL, plink2=NULL, threads=1){
  if(is.null(bfile) & is.null(pfile)){
    stop("bfile or pfile must be specified.")
  }
  if(!is.null(bfile) & !is.null(pfile)){
    stop("Both bfile and pfile cannot be specified.")
  }

  # Create object indicating tmpdir
  tmp_folder<-tempdir()

  # Determine the number of scores
  score_small<-fread(score, nrows=5)
  n_scores<-ncol(score_small)-3

  # Assemble command and files for keep and extract
  plink_opt<-paste0(plink2, ' ')
  if(!is.null(bfile)){
    plink_opt<-paste0(plink_opt, paste0('--bfile ',bfile,'CHROMOSOME_NUMBER '))
  } else {
    plink_opt<-paste0(plink_opt, paste0('--pfile ',pfile,'CHROMOSOME_NUMBER '))
  }
  if(!is.null(keep)){
    plink_opt<-paste0(plink_opt, '--keep ', keep, ' ')
  }
  if(!is.null(extract)){
    plink_opt<-paste0(plink_opt, '--extract ', extract, ' ')
  }
  if(!is.null(frq)){
    plink_opt<-paste0(plink_opt, '--read-freq ', frq,'CHROMOSOME_NUMBER.frq --score ', score,' header-read ')
  } else {
    plink_opt<-paste0(plink_opt, '--score ', score,' header-read no-mean-imputation ')
  }
  if(n_scores > 1){
    plink_opt<-paste0(plink_opt, '--score-col-nums 4-',3+n_scores, ' ')
  } else {
    plink_opt<-paste0(plink_opt, '--score-col-nums 4 ')
  }
  # Calculate score in the target sample
  for(i in chr){
    plink_opt<-gsub('CHROMOSOME_NUMBER',i,plink_opt)
    cmd_full<-paste0(plink_opt,'--chr ',i,' --out ',tmp_folder,'/profiles.chr',i,' --threads ',threads)
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
