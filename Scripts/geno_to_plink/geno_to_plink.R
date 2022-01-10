#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
suppressMessages(library("optparse"))

option_list = list(
make_option("--target", action="store", default=NA, type='character',
    help="Path to per chromosome target sample plink files [required]"),
make_option("--ref", action="store", default=NA, type='character',
    help="Path to per chromosome target sample plink files [required]"),
    make_option("--keep", action="store", default=NA, type='character',
    help="Path to a file with FID/IIDs of individuals to keep (currently only works with plink/bgen input) [optional]"),
make_option("--format", action="store", default=NA, type='character',
    help="Format of target files [required]"),
make_option("--plink2", action="store", default=NA, type='character',
    help="Path to plink2 [required]"),
make_option("--qctool2", action="store", default=NA, type='character',
    help="Path to qctool v2 [required]"),
make_option("--liftover", action="store", default=NA, type='character',
    help="Path to liftover [required]"),
make_option("--liftover_track", action="store", default=NA, type='character',
    help="Path to liftover track [required]"),
make_option("--bgen_ref", action="store", default='ref-unknown', type='character',
    help="bgen REF/ALT mode. One of [ref-first, ref-last, ref-unknown]. See plink2 documentation for details. default: 'ref-unknown'"),
    make_option("--threads", action="store", default='8', type='character',
    help="threads (used in call to plink2)"),
    make_option("--mem_mb", action="store", default='16000', type='character',
    help="memory in MB (used in call to plink2)"),
make_option("--out", action="store", default=NA, type='character',
		help="Path for output files [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

opt$out_dir<-paste0(dirname(opt$out),'/')
system(paste0('mkdir -p ',opt$out_dir))

sink(file = paste(opt$out,'.geno_to_plink.log',sep=''), append = F)
cat(
"#################################################################
# geno_to_plink.R
# For questions contact Oliver Pain (oliver.pain@kcl.ac.uk)
#################################################################
Analysis started at",as.character(start.time),'
Options are:\n')
print(opt)
cat('Analysis started at',as.character(start.time),'\n')
sink()

library(data.table)
library(bigsnpr)

###########
# Read in reference SNP data
###########

sink(file = paste(opt$out,'.geno_to_plink.log',sep=''), append = T)
cat('Reading reference SNP data...')
sink()

# Read in reference SNP data
ref<-list()
ref[['GRCh37']]<-fread(paste0(opt$ref,'.bim'))
ref[['GRCh37']]$V3<-NULL
names(ref[['GRCh37']])<-c('chr','snp','pos','a1','a2')
CHR<-ref[['GRCh37']]$chr[1]

# Create snp_modifyBuild_offline
make_executable <- function(exe) {
  Sys.chmod(exe, mode = (file.info(exe)$mode | "111"))
}

snp_modifyBuild_offline<-function (info_snp, liftOver, chain){
  if (!all(c("chr", "pos") %in% names(info_snp)))
    stop2("Please use proper names for variables in 'info_snp'. Expected %s.",
          "'chr' and 'pos'")
  liftOver <- normalizePath(liftOver)
  make_executable(liftOver)
  BED <- tempfile(fileext = ".BED")
  info_BED <- with(info_snp, data.frame(paste0("chr", chr),
                                        pos0 = pos - 1L, pos, id = rows_along(info_snp)))
  bigreadr::fwrite2(info_BED, BED, col.names = FALSE, sep = " ")
  lifted<-paste0(opt$out,'.lifted')
  unmapped<-paste0(opt$out,'.unmapped')
  system(paste(liftOver, BED, chain, lifted, unmapped))
  new_pos <- bigreadr::fread2(lifted)
  bad <- grep("^#", readLines(unmapped), value = TRUE, invert = TRUE)
  print(paste0(length(bad)," variants have not been mapped."))
  info_snp$pos <- NA
  info_snp$pos[new_pos$V4] <- new_pos$V3
  info_snp
}

# Liftover BP to GRCh38
ref[['GRCh38']]<-snp_modifyBuild_offline(ref[['GRCh37']], liftOver=opt$liftover, chain=opt$liftover_track)

names(ref[['GRCh37']])<-c('CHR','SNP','BP','A1','A2')
names(ref[['GRCh38']])<-c('CHR','SNP','BP','A1','A2')

sink(file = paste(opt$out,'.geno_to_plink.log',sep=''), append = T)
cat('Done!\n')
sink()

###################
# Read in target SNP data
###################

sink(file = paste(opt$out,'.geno_to_plink.log',sep=''), append = T)
cat('Reading target SNP data...\n')
sink()

if(opt$format == 'samp_imp_plink1'){
  target_snp<-fread(paste0(opt$target,'.bim'))
  target_snp$V3<-NULL
  names(target_snp)<-c('CHR','SNP','BP','A1','A2')
  if(!is.na(opt$keep)){
      stopifnot(file.exists(opt$keep))
  }
}

if(opt$format == 'samp_imp_bgen'){

  samplefile <- paste0(gsub('.chr.*','',opt$target),'.sample')
  if (!file.exists(samplefile)){
    samplefile <- paste0(opt$target,'.sample')
    if (!file.exists(samplefile)){
      stop('cannot determine input .sample file from prefix')
    }
  }

  if(!is.na(opt$keep)){
    stopifnot(file.exists(opt$keep))
  }
    
  sample_file<-fread(samplefile)
  sample_file<-sample_file[2,1:2]
  write.table(sample_file, paste0(opt$out,'_tmp_keep'), col.names = F, row.names = F, quote = F)
  
  # Convert bgen file to plink format containg data for one individual
  # system(paste0(opt$qctool2,' -g ',opt$target,'.bgen -s ',samplefile,' -ofiletype binary_ped -og ',opt$out,'_tmp -incl-samples ',opt$out,'_tmp_keep'))
  system(paste0(opt$plink2,' --threads ',opt$threads,' --memory ',opt$mem_mb,' --bgen ',opt$target,'.bgen ',opt$bgen_ref,' --sample ',samplefile, ' --make-bed --keep ', opt$out, '_tmp_keep --out ',opt$out,'_tmp'))
    
  target_snp<-fread(paste0(opt$out,'_tmp.bim'))
  target_snp$V3<-NULL
  names(target_snp)<-c('CHR','SNP','BP','A1','A2')
  
  system(paste0('rm ',opt$out,'_tmp.bim'))
  system(paste0('rm ',opt$out,'_tmp.bed'))
  system(paste0('rm ',opt$out,'_tmp.fam'))
  system(paste0('rm ',opt$out,'_tmp_keep'))
}


if(opt$format == 'samp_imp_vcf'){
  target_snp<-fread(cmd=paste0("zcat ",opt$target,".vcf.gz | cut -f 1-5"))
  names(target_snp)<-c('CHR','BP','SNP','A1','A2')
}

sink(file = paste(opt$out,'.geno_to_plink.log',sep=''), append = T)
cat('Done!\n')
sink()

###################
# Determine target genome build
###################

sink(file = paste(opt$out,'.geno_to_plink.log',sep=''), append = T)
cat('Determining build of target genome...\n')
sink()

# Insert IUPAC codes in ref (GRCh38)
ref[['GRCh38']]$IUPAC[ref[['GRCh38']]$A1 == 'A' & ref[['GRCh38']]$A2 =='T' | ref[['GRCh38']]$A1 == 'T' & ref[['GRCh38']]$A2 =='A']<-'W'
ref[['GRCh38']]$IUPAC[ref[['GRCh38']]$A1 == 'C' & ref[['GRCh38']]$A2 =='G' | ref[['GRCh38']]$A1 == 'G' & ref[['GRCh38']]$A2 =='C']<-'S'
ref[['GRCh38']]$IUPAC[ref[['GRCh38']]$A1 == 'A' & ref[['GRCh38']]$A2 =='G' | ref[['GRCh38']]$A1 == 'G' & ref[['GRCh38']]$A2 =='A']<-'R'
ref[['GRCh38']]$IUPAC[ref[['GRCh38']]$A1 == 'C' & ref[['GRCh38']]$A2 =='T' | ref[['GRCh38']]$A1 == 'T' & ref[['GRCh38']]$A2 =='C']<-'Y'
ref[['GRCh38']]$IUPAC[ref[['GRCh38']]$A1 == 'G' & ref[['GRCh38']]$A2 =='T' | ref[['GRCh38']]$A1 == 'T' & ref[['GRCh38']]$A2 =='G']<-'K'
ref[['GRCh38']]$IUPAC[ref[['GRCh38']]$A1 == 'A' & ref[['GRCh38']]$A2 =='C' | ref[['GRCh38']]$A1 == 'C' & ref[['GRCh38']]$A2 =='A']<-'M'

# Insert IUPAC codes in ref (ref[['GRCh37']])
ref[['GRCh37']]$IUPAC[ref[['GRCh37']]$A1 == 'A' & ref[['GRCh37']]$A2 =='T' | ref[['GRCh37']]$A1 == 'T' & ref[['GRCh37']]$A2 =='A']<-'W'
ref[['GRCh37']]$IUPAC[ref[['GRCh37']]$A1 == 'C' & ref[['GRCh37']]$A2 =='G' | ref[['GRCh37']]$A1 == 'G' & ref[['GRCh37']]$A2 =='C']<-'S'
ref[['GRCh37']]$IUPAC[ref[['GRCh37']]$A1 == 'A' & ref[['GRCh37']]$A2 =='G' | ref[['GRCh37']]$A1 == 'G' & ref[['GRCh37']]$A2 =='A']<-'R'
ref[['GRCh37']]$IUPAC[ref[['GRCh37']]$A1 == 'C' & ref[['GRCh37']]$A2 =='T' | ref[['GRCh37']]$A1 == 'T' & ref[['GRCh37']]$A2 =='C']<-'Y'
ref[['GRCh37']]$IUPAC[ref[['GRCh37']]$A1 == 'G' & ref[['GRCh37']]$A2 =='T' | ref[['GRCh37']]$A1 == 'T' & ref[['GRCh37']]$A2 =='G']<-'K'
ref[['GRCh37']]$IUPAC[ref[['GRCh37']]$A1 == 'A' & ref[['GRCh37']]$A2 =='C' | ref[['GRCh37']]$A1 == 'C' & ref[['GRCh37']]$A2 =='A']<-'M'

# Insert IUPAC codes in target
target_snp$IUPAC[target_snp$A1 == 'A' & target_snp$A2 =='T' | target_snp$A1 == 'T' & target_snp$A2 =='A']<-'W'
target_snp$IUPAC[target_snp$A1 == 'C' & target_snp$A2 =='G' | target_snp$A1 == 'G' & target_snp$A2 =='C']<-'S'
target_snp$IUPAC[target_snp$A1 == 'A' & target_snp$A2 =='G' | target_snp$A1 == 'G' & target_snp$A2 =='A']<-'R'
target_snp$IUPAC[target_snp$A1 == 'C' & target_snp$A2 =='T' | target_snp$A1 == 'T' & target_snp$A2 =='C']<-'Y'
target_snp$IUPAC[target_snp$A1 == 'G' & target_snp$A2 =='T' | target_snp$A1 == 'T' & target_snp$A2 =='G']<-'K'
target_snp$IUPAC[target_snp$A1 == 'A' & target_snp$A2 =='C' | target_snp$A1 == 'C' & target_snp$A2 =='A']<-'M'

# Check condordance of BP across builds
matched<-list()
matched[['GRCh37']]<-merge(target_snp, ref[['GRCh37']], by=c('CHR','BP','IUPAC'))
matched[['GRCh38']]<-merge(target_snp, ref[['GRCh38']], by=c('CHR','BP','IUPAC'))

sink(file = paste(opt$out,'.geno_to_plink.log',sep=''), append = T)
cat('GRCh37 match: ',round(nrow(matched[['GRCh37']])/nrow(ref[['GRCh37']])*100, 2),'%\n',sep='')
cat('GRCh38 match: ',round(nrow(matched[['GRCh38']])/nrow(ref[['GRCh38']])*100,2),'%\n',sep='')
sink()

if((nrow(matched[['GRCh37']])/nrow(ref[['GRCh37']])) > 0.7 & (nrow(matched[['GRCh37']])/nrow(ref[['GRCh37']])) > (nrow(matched[['GRCh38']])/nrow(ref[['GRCh38']]))){
  target_build<-'GRCh37'
}

if((nrow(matched[['GRCh38']])/nrow(ref[['GRCh38']])) > 0.7 & (nrow(matched[['GRCh38']])/nrow(ref[['GRCh38']])) > (nrow(matched[['GRCh37']])/nrow(ref[['GRCh37']]))){
  target_build<-'GRCh38'
}

###################
# Extract overlapping variants in plink format and insert RSIDs
###################

# To avoid issues due to duplicate IDs, we must extract variants based on original ID, update IDs manually to the reference RSID, and then extract those SNPs from the PLINK files.
extract_list_1<-matched[[target_build]]$SNP.x
extract_list_2<-matched[[target_build]]$SNP.y

write.table(extract_list_1, paste0(opt$out,'_extract_list_1.txt'), col.names = F, row.names = F, quote=F)
write.table(extract_list_2, paste0(opt$out,'_extract_list_2.txt'), col.names = F, row.names = F, quote=F)

# First extract variants based on original ID
if(opt$format == 'samp_imp_plink1'){
    
  plink_call <- paste0(opt$plink2,' --bfile ',opt$target, ' --extract ', opt$out,'_extract_list_1.txt --make-bed --memory 5000 --threads 1 --out ', opt$out,'_tmp')
    
  if (!is.na(opt$keep)){
    plink_call <- paste0(plink_call,' --keep ',opt$keep)
  }
    
  system(plink_call)
}

if(opt$format == 'samp_imp_bgen'){
    
  # --hard-call-threshold
  # system(paste0(opt$qctool2,' -g ',opt$target,'.bgen -s ',samplefile,' -ofiletype binary_ped -og ',opt$out,'_tmp -incl-rsids ',opt$out,'_extract_list_1.txt -threshold 0.9'))
    
  plink_call <- paste0(opt$plink2,' --threads ',opt$threads, ' --extract ',opt$out,'_extract_list_1.txt --memory ',opt$mem_mb,' --bgen ',opt$target,'.bgen ',opt$bgen_ref,' --sample ',samplefile, ' --hard-call-threshold 0.1 --make-bed --out ',opt$out,'_tmp')
    
  if (!is.na(opt$keep)){
    plink_call <- paste0(plink_call,' --keep ',opt$keep)
  }
    
  system(plink_call)
}

if(opt$format == 'samp_imp_vcf'){
  system(paste0(opt$plink2,' --vcf ',opt$target,'.vcf.gz --vcf-min-gq 10 --extract ', opt$out,'_extract_list_1.txt --make-bed --memory ',opt$mem_mb,' --threads ',opt$threads,' --out ', opt$out,'_tmp'))
}

# Now edit bim file to update IDs to reference IDs
targ_bim<-fread(paste0(opt$out,'_tmp.bim'))
names(targ_bim)<-c('CHR','SNP','POS','BP','A1','A2')
targ_bim$ID<-paste(targ_bim$CHR, targ_bim$BP, targ_bim$A1, targ_bim$A2, sep=':')

targ_bim_update<-targ_bim

targ_bim_update$IUPAC<-NA
targ_bim_update$IUPAC[targ_bim_update$A1 == 'A' & targ_bim_update$A2 =='T' | targ_bim_update$A1 == 'T' & targ_bim_update$A2 =='A']<-'W'
targ_bim_update$IUPAC[targ_bim_update$A1 == 'C' & targ_bim_update$A2 =='G' | targ_bim_update$A1 == 'G' & targ_bim_update$A2 =='C']<-'S'
targ_bim_update$IUPAC[targ_bim_update$A1 == 'A' & targ_bim_update$A2 =='G' | targ_bim_update$A1 == 'G' & targ_bim_update$A2 =='A']<-'R'
targ_bim_update$IUPAC[targ_bim_update$A1 == 'C' & targ_bim_update$A2 =='T' | targ_bim_update$A1 == 'T' & targ_bim_update$A2 =='C']<-'Y'
targ_bim_update$IUPAC[targ_bim_update$A1 == 'G' & targ_bim_update$A2 =='T' | targ_bim_update$A1 == 'T' & targ_bim_update$A2 =='G']<-'K'
targ_bim_update$IUPAC[targ_bim_update$A1 == 'A' & targ_bim_update$A2 =='C' | targ_bim_update$A1 == 'C' & targ_bim_update$A2 =='A']<-'M'

targ_bim_update<-merge(targ_bim_update, ref[[target_build]], by.x=c('CHR','BP','IUPAC'), by.y=c('CHR','BP','IUPAC'), all.x=T)

# Sort bim to be in the original order
targ_bim_update<-targ_bim_update[match(targ_bim$ID, targ_bim_update$ID),]

# Set ID for missing variants to the original ID
targ_bim_update$SNP.y[is.na(targ_bim_update$SNP.y)]<-targ_bim_update$SNP.x[is.na(targ_bim_update$SNP.y)]
# Set ID to containg _dup if has a duplicate RSID
targ_bim_update$SNP.y[duplicated(targ_bim_update$SNP.y)]<-paste0(targ_bim_update$SNP.y[duplicated(targ_bim_update$SNP.y)],'_dup')

targ_bim_update_clean<-targ_bim_update[,c('CHR','SNP.y','POS','BP','A1.x','A2.x'),with=F]
setnames(targ_bim_update_clean, c('CHR','SNP','POS','BP','A1','A2'))

fwrite(targ_bim_update_clean, paste0(opt$out,'_tmp.bim'), col.names=F, row.names=F, quote=F, sep=' ')

# Extract variants based on new reference RSIDs
system(paste0(opt$plink2,' --bfile ',opt$out,'_tmp --extract ', opt$out,'_extract_list_2.txt --make-bed --memory ',opt$mem_mb,' --threads ',opt$threads,' --out ', opt$out))

system(paste0('rm ', opt$out,'.log'))
system(paste0('rm ', opt$out,'_extract*'))
system(paste0('rm ', opt$out,'_tmp*'))
system(paste0('rm ', opt$out,'.unmapped'))
system(paste0('rm ', opt$out,'.lifted'))

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = paste(opt$out,'.geno_to_plink.log',sep=''), append = T)
cat('Analysis finished at',as.character(end.time),'\n')
cat('Analysis duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
sink()
