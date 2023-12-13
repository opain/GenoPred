#!/usr/bin/Rscript

# Calculate scores in target file
calc_score<-function(bfile, score, keep=NULL, extract=NULL, chr=1:22, frq=NULL, plink2='plink2', threads=1){
    # Create object indicating tmpdir
    tmp_folder<-tempdir()

    # Determine the number of scores
    score_small<-fread(score, nrows=5)
    n_scores<-ncol(score_small)-2
    
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
        cmd<-paste0(cmd, ' --score-col-nums 3-',2+n_scores)
    }
    # Calculate score in the target sample
    for(i in chr){
        cmd_i<-gsub('CHROMOSOME_NUMBER',i,cmd)
        system(paste0(plink2, ' --bfile ',bfile,i,cmd_i,' --chr ',i,' --out ',tmp_folder,'/profiles.chr',i,' --threads ',threads))
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
    nsnp<-NULL
    for(i in chr){
        if(file.exists(paste0(tmp_folder,'/profiles.chr',i,'.sscore'))){
            sscore<-fread(paste0(tmp_folder,'/profiles.chr',i,'.sscore'))
            scores[[i]]<-sscore[,grepl('_AVG$', names(sscore)),with=F]
            # This allows for difference plink formats across plink versions
            if(any(names(sscore) == 'NMISS_ALLELE_CT')){
                scores[[i]]<-as.matrix(scores[[i]]*sscore$NMISS_ALLELE_CT[1]/2)
                nsnp<-c(nsnp, sscore$NMISS_ALLELE_CT[1]/2)
            } else {
                scores[[i]]<-as.matrix(scores[[i]]*sscore$ALLELE_CT[1]/2)
                nsnp<-c(nsnp, sscore$ALLELE_CT[1]/2)
            }
        } else {
            cat0('No scores for chromosome ',i,'. Check plink logs file for reason.\n')
        }
    }
    nsnp<-sum(nsnp)
    
    # Remove NULL elements from list (these are inserted by R when list objects are numbered)
    scores[sapply(scores, is.null)] <- NULL
    
    # sum scores across chromosomes
    scores<-Reduce(`+`, scores)
    
    # Combine score with IDs
    scores<-data.table(scores_ids,
                       scores)
    
    # Rename columns
    names(scores)<-c('FID','IID',names(score_small)[-1:-2])
    
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
dat_or_file<-function(x, header = F){
  if(is.null(x)){
    file_path <- NULL
  } else {
    if (is.data.frame(x)) {
      # Create file in temp directory
      file_path<-tempfile()
      write.table(x, file_path, col.names = header, row.names=F, quote=F)
    } else if (is.character(x) && file.exists(x)) {
      file_path <- x
    } else {
      stop("Input must be a data frame or a valid file path.")
    }
  }
  return(file_path)
}

# Make a subset of plink1 binaries
plink_subset<-function(plink='plink', chr = 1:22, keep = NA, bfile, out, memory = 4000, threads = 1){

  # If object, create file
  keep <- dat_or_file(keep)

  # Prepare plink options
  plink_opt<-NULL
  if(!is.null(keep)){
    plink_opt<-paste0(plink_opt, paste0('--keep ',keep,' '))
  }

  # Run plink
  for(chr_i in chr){
    cmd <- paste0(opt$plink,' --bfile ',bfile,chr_i,' --threads ', threads,' ',plink_opt,'--make-bed --out ',out,chr_i,' --memory ',memory)
    exit_status <- system(cmd, intern=FALSE)
    if (exit_status != 0) {
      cat("Error occurred in running plink command for chromosome",chr_i,"\n")
    }
  }
}

# Return a list of variants surviving QC using plink
plink_qc_snplist<-function(bfile, chr = 1:22, threads = 1, memory = 4000, geno = NULL, maf = NULL, hwe = NULL){
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
    system(paste0(opt$plink,' --bfile ',bfile,chr_i,' --threads ',threads,' ',plink_opt,' --write-snplist --out ', temp_file,' --memory ', memory))
    snplist<-c(snplist, fread(paste0(temp_file, '.snplist'), header=F)$V1)
  }

  return(snplist)
}

# Coordinates of high ld regions
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

# Read in per chromosome bim files
read_bim<-function(x){
  bim<-NULL
  for(chr_i in 1:22){
      bim <- rbind(bim, fread(paste0(x, chr_i,'.bim')))
  }
  names(bim)<-c('CHR','SNP','POS','BP','A1','A2')
  return(bim)
}

# Perform PCA using plink files
plink_pca<-function(bfile, extract = NULL, flip = NULL, memory = 4000, n_pc = 6){
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
  cmd<-paste0(opt$plink,' --merge-list ',tmp_dir,'/ref_mergelist.txt ',plink_opt,'--threads 1 --make-bed --out ',tmp_dir,'/ref_merge --memory ',memory)
  system(cmd)

  # Calculate SNP weights
  system(paste0(opt$plink2,' --bfile ',tmp_dir,'/ref_merge --threads 1 --pca ',n_pc,' biallelic-var-wts  --out ',tmp_dir,'/ref_merge --memory ', memory))

  # Format the SNP-weights
  snp_weights<-fread(paste0(tmp_dir,'/ref_merge.eigenvec.var'))
  snp_weights<-snp_weights[,c(-1,-4),with=F]
  names(snp_weights)[1:2]<-c('SNP','A1')

  return(snp_weights)
}