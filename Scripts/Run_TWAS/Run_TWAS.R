#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
suppressMessages(library("optparse"))

option_list = list(
make_option("--gwas", action="store", default=NA, type='character',
	help="Name of GWAS summary stats [required]"),
make_option("--pos", action="store", default=NA, type='character',
		help="File listing SNP-weights in .pos format [required]"),
make_option("--out", action="store", default=NA, type='character',
			help="Name of output files [required]"),
make_option("--fusion_dir", action="store", default=NA, type='character',
      help="Directory containing fusion software and reference data [required]"),
make_option("--coloc_P", action="store", default=NA, type='numeric',
      help="Specify as T to perform colocalisation [optional]"),
make_option("--ncores", action="store", default=5, type='numeric',
			help="Number of cores for parallel analysis [optional]")
)

opt = parse_args(OptionParser(option_list=option_list))

suppressMessages(library(foreach))
suppressMessages(library(doMC))
suppressMessages(library(data.table))
registerDoMC(opt$ncores)

# Make the folder for the output
system(paste0('mkdir -p ',opt$out,'_logs'))

sink(file = paste0(opt$out,'.log'), append = F)
cat(
'#################################################################
# Run_TWAS.R
# For questions contact Oliver Pain (oliver.pain@kcl.ac.uk)
#################################################################
Analysis started at',as.character(start.time),'
Options are:\n')
print(opt)
sink()

# Remove blank lines from sumstats
sink(file = paste0(opt$out,'.log'), append = T)
cat('Removing rows containing blanks from gwas file...')
sink()
sumstats<-data.frame(fread(cmd=paste0('zcat ',opt$gwas)))
sumstats<-sumstats[complete.cases(sumstats),]
sumstat_out <- paste0(opt$out,'_noNA.sumstats')
GWASN<-round(median(sumstats$N),0)
write.table(sumstats, file=sumstat_out, col.names=T,row.names=F, quote=F)
system(paste0('gzip ',sumstat_out))
sink(file = paste0(opt$out,'.log'), append = T)
cat('Done!\n')
sink()
opt$gwas_new<-paste0(opt$out,'_noNA.sumstats.gz')

# Print the FUSION command to the log files
sink(file = paste0(opt$out,'.log'), append = T)
if(is.na(opt$coloc_P)){
  cat('
FUSION command:
	Rscript ',opt$fusion_dir,'/fusion_twas-master/FUSION.assoc_test.R \\
		--sumstats ',opt$gwas_new,' \\
		--weights ',opt$pos,' \\
		--weights_dir ',opt$fusion_dir,'/SNP-weights \\
		--ref_ld_chr ',opt$fusion_dir,'/LDREF/1000G.EUR. \\
		--out ',opt$out,'_res_chr$i.txt \\
		--chr $i

',sep='')
} else {
  cat('
FUSION command:
	Rscript ',opt$fusion_dir,'/fusion_twas-master/FUSION.assoc_test.R \\
		--sumstats ',opt$gwas_new,' \\
		--weights ',opt$pos,' \\
		--weights_dir ',opt$fusion_dir,'/SNP-weights \\
		--ref_ld_chr ',opt$fusion_dir,'/LDREF/1000G.EUR. \\
		--out ',opt$out,'_res_chr$i.txt \\
		--coloc_P ',opt$coloc_P,' \\
		--GWASN ',GWASN,' \\
		--chr $i

',sep='')
}
sink()

sink(file = paste0(opt$out,'.log'), append = T)
	cat('Running TWAS across', opt$ncores,'cores ... ')
sink()

# Run TWAS across all chromosomes in parallel
is.odd <- function(x){x %% 2 != 0}
CHROMS<-c(1:22, rep(NA,30)) 
CHROMS_mat<-matrix(CHROMS,nrow=opt$ncores, ncol=ceiling(22/opt$ncores)) 
for(i in which(is.odd(1:dim(CHROMS_mat)[2]))){CHROMS_mat[,i]<-rev(CHROMS_mat[,i])} 
CHROMS<-as.numeric(CHROMS_mat) 
CHROMS<-CHROMS[!is.na(CHROMS)] 
print(CHROMS)

if(is.na(opt$coloc_P)){
  TWAS_log<-foreach(i=CHROMS) %dopar% {
    log<-system(paste0(
      'Rscript ',opt$fusion_dir,'/fusion_twas-master/FUSION.assoc_test.R ',
      '--sumstats ',opt$gwas_new,' ',
      '--weights ',opt$pos,' ',
      '--weights_dir ',opt$fusion_dir,'/SNP-weights ',
      '--ref_ld_chr ',opt$fusion_dir,'/LDREF/1000G.EUR. ',
      '--out ',opt$out,'_res_chr',i,'.txt ',
      '--chr ',i
    ),intern = TRUE)
  	write.table(log, paste0(opt$out,'_logs/chr',i,'.log'), col.names=F, row.names=F, quote=F)
  	log
  }
} else {
  TWAS_log<-foreach(i=CHROMS) %dopar% {
    log<-system(paste0(
      'Rscript ',opt$fusion_dir,'/fusion_twas-master/FUSION.assoc_test.R ',
      '--sumstats ',opt$gwas_new,' ',
      '--weights ',opt$pos,' ',
      '--weights_dir ',opt$fusion_dir,'/SNP-weights ',
      '--ref_ld_chr ',opt$fusion_dir,'/LDREF/1000G.EUR. ',
      '--out ',opt$out,'_res_chr',i,'.txt ',
      '--coloc_P ',opt$coloc_P,' ',
      '--GWASN ',GWASN,' ',
      '--chr ',i
    ),intern = TRUE)
    write.table(log, paste0(opt$out,'_logs/chr',i,'.log'), col.names=F, row.names=F, quote=F)
    log
  }
}
sink(file = paste0(opt$out,'.log'), append = T)
	cat('Done!\n')
sink()

# Check there are 23 results files (There is one extra for the MHC region).
temp<-list.files(pattern = paste0(opt$out,'_res_*'))
if(length(temp) != 23){
	sink(file = paste0(opt$out,'.log'), append = T)
		cat('There are only',length(temp),'results files!\n')
	sink()
}

# Combine the TWAS results across all chromomsomes
system(paste0('head -1 ',opt$out,'_res_chr',22,'.txt > ',opt$out,'_res_GW.txt'))
system(paste0('tail -n +2 -q ',opt$out,'_res_chr* >> ',opt$out,'_res_GW.txt'))

# Delete the per chromosome files
system(paste0('rm ',opt$out,'_res_chr*'))

# Read in the twas results to count how many genes were skipped
res<-read.table(paste0(opt$out,'_res_GW.txt'), header=T, stringsAsFactors=F)
sink(file = paste0(opt$out,'.log'), append = T)
	cat('In total,',sum(is.na(res$TWAS.P)),'/', dim(res)[1],'features were skipped.\n')
sink()

end.time <- Sys.time()
time.taken <- end.time - start.time

sink(file = paste0(opt$out,'.log'), append = T)
cat('Analysis finished at',as.character(end.time),'\n')
cat('Analysis duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),sep=,'\n')
sink()
