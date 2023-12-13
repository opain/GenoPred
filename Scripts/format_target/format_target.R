#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
library("optparse")

option_list = list(
make_option("--target", action="store", default=NA, type='character',
    help="Prefix to a target sample plink file for a given chromosome [required]"),
make_option("--ref", action="store", default=NA, type='character',
    help="Prefix to a reference sample plink/.rds file [required]"),
make_option("--format", action="store", default=NA, type='character',
    help="Format of target files. Either samp_imp_plink1, samp_imp_bgen, or samp_imp_vcf. [required]"),
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
source('../Scripts/functions/misc.R')

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

# Identify variants that need to be flipped
names(ref_target)[names(ref_target) == 'IUPAC']<-'IUPAC.y'
ref_target$IUPAC.x<-snp_iupac(ref_target$A1.x, ref_target$A2.x)
flip <- detect_strand_flip(ref_target$IUPAC.x, ref_target$IUPAC.y)
if(sum(flip) > 0){
  write.table(ref_target$SNP.y[flip], paste0(tmp_dir,'/flip_list.txt'), col.names = F, row.names = F, quote=F)
}

# Remove variants where IUPAC codes do not match (allowing for strand flips)
matched <- which((ref_target$IUPAC.x == ref_target$IUPAC.y) | flip)
ref_target<-ref_target[matched,]

log_add(log_file = log_file, message = paste0('Target contains ', nrow(ref_target),' reference variants.'))

# To avoid issues due to duplicate IDs, we must extract variants based on original ID, update IDs manually to the reference RSID, and then extract those SNPs from the PLINK files.
write.table(ref_target$SNP.x, paste0(tmp_dir,'/extract_list_1.txt'), col.names = F, row.names = F, quote=F)
write.table(ref_target$SNP.y, paste0(tmp_dir,'/extract_list_2.txt'), col.names = F, row.names = F, quote=F)

# First extract variants based on original ID
if(opt$format == 'samp_imp_plink1'){
  system(paste0(opt$plink2,' --bfile ',opt$target, ' --extract ', tmp_dir,'/extract_list_1.txt --make-bed --memory 5000 --threads 1 --out ', tmp_dir,'/subset'))
}
if(opt$format == 'samp_imp_bgen'){
  system(paste0(opt$plink2,' --bgen ',opt$target,'.bgen ref-last --sample ',gsub('.chr.*','',opt$target),'.sample --import-dosage-certainty 0.9 --extract ', tmp_dir,'/extract_list_1.txt --make-bed --memory 5000 --threads 1 --out ', tmp_dir,'/subset'))
}
if(opt$format == 'samp_imp_vcf'){
  system(paste0(opt$plink2,' --vcf ',opt$target,'.vcf.gz --vcf-min-gq 10 --import-dosage-certainty 0.9 --extract ', tmp_dir,'/extract_list_1.txt --make-bed --memory 5000 --threads 1 --out ', tmp_dir,'/subset'))
}

# Now edit bim file to update IDs to reference IDs
targ_bim<-fread(paste0(tmp_dir,'/subset.bim'))
names(targ_bim)<-c('CHR','SNP','POS','BP','A1','A2')

# Update SNP if original SNP value is present in the previously matched ref and target data
targ_bim$SNP.x<-targ_bim$SNP
targ_bim[ref_target, on=.(SNP.x), SNP := i.SNP.y]
targ_bim$SNP.x<-NULL

# Label SNP with _dup if the RSID is duplicated, so these variants are removed.
dup_snp<-duplicated(targ_bim$SNP)
log_add(log_file = log_file, message = paste0('Removing ', sum(dup_snp),' duplicate variants.'))
targ_bim$SNP[dup_snp]<-paste0(targ_bim$SNP[dup_snp],'_dup')

log_add(log_file = log_file, message = paste0(sum(!dup_snp)," of ", nrow(ref)," reference variants are in the target."))

# Write out new bim file
fwrite(targ_bim, paste0(tmp_dir,'/subset.bim'), col.names=F, row.names=F, quote=F, na='NA', sep=' ')

# Extract variants based on new reference RSIDs
# and flip variants if there are any to be flipped
plink_opt<-NULL
if(sum(flip) > 0){
  plink_opt<-c(plink_opt, paste0('--flip ',tmp_dir,'/flip_list.txt '))
}

system(paste0(opt$plink,' --bfile ',tmp_dir,'/subset ',plink_opt,'--extract ', tmp_dir,'/extract_list_2.txt --make-bed --memory 5000 --threads 1 --out ', tmp_dir,'/subset'))

##################
# Insert missing SNPs into the reference data
##################

log_add(log_file = log_file, message = 'Inserting missing reference variants.')

# Update IDs in reference to avoid conflict with the target
ref_fam<-fread(paste0(opt$ref,'.fam'))
ref_ID_update<-data.frame(ref_fam$V1, ref_fam$V2, paste0(ref_fam$V1,'_REF'),paste0(ref_fam$V2,'_REF'))
fwrite(ref_ID_update, paste0(tmp_dir,'/ref_ID_update.txt'), sep=' ', col.names=F)
system(paste0(opt$plink,' --bfile ',opt$ref,' --make-bed --update-ids ',tmp_dir,'/ref_ID_update.txt --out ',tmp_dir,'/REF --memory 7000'))

# Merge target and reference plink files to insert missing SNPs
system(paste0(opt$plink,' --bfile ',tmp_dir,'/subset --bmerge ',tmp_dir,'/REF --make-bed --out ',tmp_dir,'/subset'))

# Extract only target individuals
system(paste0(opt$plink,' --bfile ',tmp_dir,'/subset --remove ',tmp_dir,'/REF.fam --make-bed --out ',opt$output))

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = log_file, append = T)
cat('Analysis finished at',as.character(end.time),'\n')
cat('Analysis duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
sink()
