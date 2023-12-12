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