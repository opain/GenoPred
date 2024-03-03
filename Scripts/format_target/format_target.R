#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
library("optparse")

option_list = list(
make_option("--target", action="store", default=NA, type='character',
    help="Prefix to a target sample plink files for a given chromosome [required]"),
make_option("--ref", action="store", default=NA, type='character',
    help="Prefix to a reference sample plink and .rds files for a given chromosome [required]"),
make_option("--format", action="store", default=NA, type='character',
    help="Format of target files. Either plink1, plink2, bgen, or vcf. [required]"),
make_option("--plink", action="store", default='plink', type='character',
    help="Path to plink1.9 [optional]"),
make_option("--plink2", action="store", default='plink2', type='character',
    help="Path to plink2 [optional]"),
make_option("--output", action="store", default=NA, type='character',
		help="Path for output files [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

# Load dependencies
library(GenoUtils)
library(data.table)
library(RSQLite)
source('../functions/misc.R')
source_all('../functions')

# Create output directory
opt$out_dir<-paste0(dirname(opt$output),'/')
system(paste0('mkdir -p ',opt$out_dir))

# Create temp directory
tmp_dir<-tempdir()

# Initiate log file
log_file <- paste0(opt$output,'.format_target.log')
log_header(log_file = log_file, opt = opt, script = 'format_target.R', start.time = start.time)

###########
# Read in reference SNP data
###########

log_add(log_file = log_file, message = 'Reading in reference SNP data.')

# Read in reference SNP data
ref<-readRDS(paste0(opt$ref,'.rds'))

log_add(log_file = log_file, message = paste0('Reference data contains ', nrow(ref),' variants.'))

###################
# Read in target SNP data
###################

log_add(log_file = log_file, message = 'Reading in target SNP data.')

target_snp<-read_geno(target = opt$target, format = opt$format)

log_add(log_file = log_file, message = paste0('Target data contains ', nrow(target_snp),' variants.'))

# Throw an error if there are not enough SNPs in the target data
if(nrow(target_snp) < nrow(ref)*0.9){
  log_add(log_file = log_file, message = 'Error: Check your target data has been imputed already.')
  stop("Check your target data has been imputed already.")
}

###################
# Determine target genome build
###################

target_build <- detect_build( ref = ref,
                              targ = target_snp,
                              log_file = log_file)

###################
# Extract overlapping variants in plink format and insert RSIDs
###################

# Subset relevent build from ref
ref_subset<-ref[, c("CHR","SNP",paste0("BP_",target_build),"A1","A2","IUPAC"), with=F]
names(ref_subset)[names(ref_subset) == paste0("BP_",target_build)]<-'BP'

# Merge target and ref by CHR and BP
ref_target<-merge(target_snp, ref_subset, by=c('CHR','BP'))

# Insert IUPAC codes
names(ref_target)[names(ref_target) == 'IUPAC']<-'IUPAC.y'
ref_target$IUPAC.x<-snp_iupac(ref_target$A1.x, ref_target$A2.x)

# Only retain variants that are non-ambiguous SNPs in the target
# Some might have been retained due to matching CHR and BP position with SNPs in reference
ref_target<-ref_target[!is.na(ref_target$IUPAC.x),]

# Identify variants that need to be flipped
flip <- detect_strand_flip(ref_target$IUPAC.x, ref_target$IUPAC.y)
if(any(flip)){
  write.table(ref_target$SNP.y[flip], paste0(tmp_dir,'/flip_list.txt'), col.names = F, row.names = F, quote=F)
}

# Remove variants where IUPAC codes do not match (allowing for strand flips)
matched <- which((ref_target$IUPAC.x == ref_target$IUPAC.y) | flip)
ref_target<-ref_target[matched,]

log_add(log_file = log_file, message = paste0('Target contains ', nrow(ref_target),' reference variants.'))

if(nrow(ref_target) < nrow(ref)*0.7){
  log_add(log_file = log_file, message = 'Error: Less than 70% of reference variants are present in the target.')
  stop("Less than 70% of reference variants are present in the target")
}

# To avoid issues due to duplicate IDs, we must extract variants based on original ID, update IDs manually to the reference RSID, and then extract those SNPs from the PLINK files.
write.table(ref_target$SNP.x, paste0(tmp_dir,'/extract_list_1.txt'), col.names = F, row.names = F, quote=F)
write.table(ref_target$SNP.y, paste0(tmp_dir,'/extract_list_2.txt'), col.names = F, row.names = F, quote=F)

# First extract variants based on original ID
if(opt$format == 'plink1'){
  system(paste0(opt$plink2,' --bfile ',opt$target, ' --extract ', tmp_dir,'/extract_list_1.txt --make-pgen --memory 5000 --threads 1 --out ', tmp_dir,'/subset'))
}
if(opt$format == 'plink2'){
  system(paste0(opt$plink2,' --pfile ',opt$target, ' --extract ', tmp_dir,'/extract_list_1.txt --make-pgen --memory 5000 --threads 1 --out ', tmp_dir,'/subset'))
}
if(opt$format == 'bgen'){
  system(paste0(opt$plink2,' --bgen ',opt$target,'.bgen ref-last --sample ',gsub('.chr.*','',opt$target),'.sample --extract ', tmp_dir,'/extract_list_1.txt --make-pgen --memory 5000 --threads 1 --out ', tmp_dir,'/subset'))
}
if(opt$format == 'vcf'){
  system(paste0(opt$plink2,' --vcf ',opt$target,'.vcf.gz --extract ', tmp_dir,'/extract_list_1.txt --make-pgen --memory 5000 --threads 1 --out ', tmp_dir,'/subset'))
}

# Ensure both FID and IID are present in the .psam file
targ_psam <- fread(paste0(tmp_dir,'/subset.psam'))
names(targ_psam)<-gsub('\\#', '', names(targ_psam))
if(sum(names(targ_psam) %in% c('FID', 'IID')) == 1){
  targ_psam$FID <- targ_psam$IID
  targ_psam <- targ_psam[, c('FID','IID', names(targ_psam)[!(names(targ_psam) %in% c('FID','IID'))]), with=F]
  names(targ_psam)[1]<-paste0('#FID')
  fwrite(targ_psam, paste0(tmp_dir,'/subset.psam'), col.names=T, row.names=F, quote=F, na='NA', sep=' ')
}

# Now edit bim file to update IDs to reference IDs
targ_pvar<-fread(paste0(tmp_dir,'/subset.pvar'))
names(targ_pvar)<-c('CHR','BP','SNP','A2','A1')

# Update SNP with reference SNP value based on CHR:BP:IUPAC in the previously matched ref and target data
targ_pvar$IUPAC<-snp_iupac(targ_pvar$A1, targ_pvar$A2)
targ_pvar$ID<-paste0(targ_pvar$CHR,':',targ_pvar$BP,':',targ_pvar$IUPAC)
targ_pvar$SNP<-targ_pvar$ID # Give SNP column a unique value before updating to reference value
ref_target$ID<-paste0(ref_target$CHR,':',ref_target$BP,':',ref_target$IUPAC.x)
targ_pvar[ref_target, on=.(ID), SNP := i.SNP.y]
targ_pvar<-targ_pvar[,c('CHR','BP','SNP','A2','A1'),with=F]

# Label SNP with _dup if the RSID is duplicated, so these variants are removed.
dup_snp<-duplicated(targ_pvar$SNP)
log_add(log_file = log_file, message = paste0('Removing ', sum(dup_snp),' duplicate variants - May have IUPAC NA.'))
targ_pvar$SNP[dup_snp]<-paste0(targ_pvar$SNP[dup_snp],'_dup')

# Write out new bim file
names(targ_pvar)<-c('#CHROM','POS','ID','REF','ALT')
fwrite(targ_pvar, paste0(tmp_dir,'/subset.pvar'), col.names=T, row.names=F, quote=F, na='NA', sep=' ')

# Extract variants based on new reference RSIDs
# and flip variants if there are any to be flipped
plink_opt<-NULL
if(sum(flip) > 0){
  plink_opt<-paste0(plink_opt, paste0('--flip ',tmp_dir,'/flip_list.txt '))
}

system(paste0(opt$plink2,' --pfile ',tmp_dir,'/subset ',plink_opt,'--extract ', tmp_dir,'/extract_list_2.txt --make-pgen --memory 5000 --threads 1 --out ', tmp_dir,'/subset'))

##################
# Insert missing SNPs into the reference data
##################

log_add(log_file = log_file, message = 'Inserting missing reference variants.')

# Update IDs in reference to avoid conflict with the target
ref_psam<-fread(paste0(opt$ref,'.psam'))
names(ref_psam)<-gsub('\\#', '', names(ref_psam))
ref_psam <- ref_psam[, names(ref_psam) %in% c('FID', 'IID'), with = F]
if(ncol(ref_psam) == 1){
  ref_ID_update<-data.frame(ref_psam$`IID`, paste0(ref_psam$`IID`,'_REF'))
} else {
  ref_ID_update<-data.frame(ref_psam$`FID`, ref_psam$`IID`, paste0(ref_psam$`FID`,'_REF'), paste0(ref_psam$`IID`,'_REF'))
}
fwrite(ref_ID_update, paste0(tmp_dir,'/ref_ID_update.txt'), sep=' ', col.names=F)
system(paste0(opt$plink2,' --pfile ',opt$ref,' --make-pgen --update-ids ',tmp_dir,'/ref_ID_update.txt --out ',tmp_dir,'/REF --memory 5000 --threads 1'))

# Merge target and reference plink files to insert missing SNPs
# plink2's pmerge only handles concatenation for the time being
# In the meantime, convert the ref and target into plink1 binaries, merge, and then convert back to plink2 binaries
system(paste0(opt$plink2,' --pfile ',tmp_dir,'/subset --make-bed --memory 5000 --threads 1 --out ',tmp_dir,'/subset'))
system(paste0(opt$plink2,' --pfile ',tmp_dir,'/REF --make-bed --out ',tmp_dir,'/REF --memory 5000 --threads 1'))
system(paste0(opt$plink,' --bfile ',tmp_dir,'/subset --bmerge ',tmp_dir,'/REF --make-bed --allow-no-sex --out ',tmp_dir,'/ref_targ'))

# Extract only target individuals
system(paste0(opt$plink2,' --bfile ',tmp_dir,'/ref_targ --remove ',tmp_dir,'/REF.psam --make-pgen --memory 5000 --threads 1 --out ',opt$output))

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = log_file, append = T)
cat('Analysis finished at',as.character(end.time),'\n')
cat('Analysis duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
sink()
