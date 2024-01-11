#############
# Use OpenSNP data in GenoPred for testing
#############

# Download OpenSNP data
mkdir -p /users/k1806347/oliverpainfel/Data/OpenSNP/raw
cd /users/k1806347/oliverpainfel/Data/OpenSNP/raw
wget https://zenodo.org/records/1442755/files/CrowdAI_v3.tar.gz?download=1 
mv 'CrowdAI_v3.tar.gz?download=1' CrowdAI_v3.tar.gz
tar -xvzf CrowdAI_v3.tar.gz
rm CrowdAI_v3.tar.gz

# Use fullset (training) data as the target for GenoPred, as the sample size is larger

##
# Prepare phenotype data
##
R

library(data.table)
dat<-fread('/users/k1806347/oliverpainfel/Data/OpenSNP/raw/CrowdAI_v3/training_set_details.txt')
dat$FID<-0
dat$IID<-dat$id
dat<-dat[,c('FID','IID','height'), with=F]

dir.create('/users/k1806347/oliverpainfel/Data/OpenSNP/processed/pheno', recursive=T)
write.table(dat, '/users/k1806347/oliverpainfel/Data/OpenSNP/processed/pheno/height.txt', col.names=T, row.names=F, quote=F)

q()
n

##
# Split the vcf files into per chromosome files
##
# Create index
module add bcftools/1.14-gcc-10.3.0-python3+-chk-version
bcftools index /users/k1806347/oliverpainfel/Data/OpenSNP/raw/CrowdAI_v3/fullset/genotyping_data_fullset_train.vcf.gz

# Now, split by chromosome using plink2
# Run on the command line within pipeline conda environment
mkdir /users/k1806347/oliverpainfel/Data/OpenSNP/processed/geno
for chr in $(seq 1 22);do
    plink2 \
        --vcf /users/k1806347/oliverpainfel/Data/OpenSNP/raw/CrowdAI_v3/fullset/genotyping_data_fullset_train.vcf.gz \
        --chr ${chr} \
        --out /users/k1806347/oliverpainfel/Data/OpenSNP/processed/geno/opensnp_train.chr${chr} \
        --export vcf bgz
done

##
# Create a config file for opensnp test analyses
##

R
setwd('/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/GenoPredPipe')
conf<-readLines('config.yaml')

conf<-gsub('outdir: test_data/output/test1', 'outdir: /users/k1806347/oliverpainfel/Data/OpenSNP/GenoPred/test1', conf)
conf<-gsub('config_file: config.yaml', 'config_file: misc/opensnp/config.yaml', conf)
conf<-gsub('gwas_list: example_input/gwas_list_example.txt', 'gwas_list: misc/opensnp/gwas_list.txt', conf)
conf<-gsub('target_list: example_input/target_list_example.txt', 'target_list: misc/opensnp/target_list.txt', conf)
conf<-conf[conf != '']
conf<-conf[!grepl('\\#', conf)]
conf<-conf[!grepl('score_list', conf)]

dir.create('misc/opensnp/', recursive = T)
dir.create('/users/k1806347/oliverpainfel/Data/OpenSNP/GenoPred', recursive = T)
write.table(conf, 'misc/opensnp/config.yaml', col.names=F, row.names=F, quote=F)

##
# Create target_list file
##

library(data.table)
target_list<-fread('example_input/target_list_example.txt')
target_list<-rbind(target_list, data.table(
    name='opensnp',
    path='/users/k1806347/oliverpainfel/Data/OpenSNP/processed/geno/opensnp_train',
    type='samp_imp_vcf',
    indiv_report=F))

target_list<-target_list[target_list$name == 'opensnp',]

write.table(target_list, 'misc/opensnp/target_list.txt', col.names=T, row.names=F, quote=F, sep=' ')

q()
n

##
# Download height GWAS sumstats
##
# These are from the Yengo 2022 paper
mkdir -p /users/k1806347/oliverpainfel/Data/GWAS_sumstats/opensnp_test
wget --no-check-certificate -O /users/k1806347/oliverpainfel/Data/GWAS_sumstats/opensnp_test/yengo_2022_height_eur.txt https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90245001-GCST90246000/GCST90245992/GCST90245992_buildGRCh37.tsv
wget --no-check-certificate -O /users/k1806347/oliverpainfel/Data/GWAS_sumstats/opensnp_test/yengo_2022_height_eas.txt https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90245001-GCST90246000/GCST90245991/GCST90245991_buildGRCh37.tsv

##
# Create gwas_list file
##
R

setwd('/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/GenoPredPipe')

gwas_list<-fread('example_input/gwas_list_example.txt')
gwas_list<-rbind(gwas_list, data.table(
    name='yengo_eur',
    path='/users/k1806347/oliverpainfel/Data/GWAS_sumstats/opensnp_test/yengo_2022_height_eur.txt',
    population='EUR',
    n=NA,
    sampling=NA,
    prevalence=NA,
    mean=NA,
    sd=NA,
    label="\"Yengo 2022 Height EUR\""))

gwas_list<-rbind(gwas_list, data.table(
    name='yengo_eas',
    path='/users/k1806347/oliverpainfel/Data/GWAS_sumstats/opensnp_test/yengo_2022_height_eas.txt',
    population='EAS',
    n=NA,
    sampling=NA,
    prevalence=NA,
    mean=NA,
    sd=NA,
    label="\"Yengo 2022 Height EAS\""))

gwas_list<-gwas_list[gwas_list$name %in% c('yengo_eur','yengo_eas'),]

write.table(gwas_list, 'misc/opensnp/gwas_list.txt', col.names=T, row.names=F, quote=F, sep=' ')

q()
n

####
# Download a precomputed score file
####

wget --no-check-certificate -O /users/k1806347/oliverpainfel/Data/GWAS_sumstats/opensnp_test/yengo_2022_height_eur_pgs.gz https://cnsgenomics.com/data/giant_2022/GIANT_HEIGHT_YENGO_2022_PGS_WEIGHTS_EUR.gz
# Going to pause on this for now until I have figured out a way to harmonise score files.

#######
# Run the snakmake pipeline
#######

# Run on the command line
# Do a dry run to avoid recreating reference data
snakemake -n --use-conda --configfile misc/opensnp/config.yaml output_all

#######
# Test correlation between PGS and observed height
#######

R
setwd('/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/GenoPredPipe/')
library(data.table)
library(Hmisc)

# Read in pheno data
pheno <- fread('/users/k1806347/oliverpainfel/Data/OpenSNP/processed/pheno/height.txt')
pheno$FID<-0

# Read in ancestry data
keep_list <- fread('/users/k1806347/oliverpainfel/Data/OpenSNP/GenoPred/test1/opensnp/ancestry/keep_list.txt')

# Read in pgs
gwas_list <- fread('misc/opensnp/gwas_list.txt')
pgs_methods <- c('ptclump','dbslmm','prscs','sbayesr','lassosum','ldpred2','megaprs')
pgs_methods_eur <- c('ptclump','lassosum','megaprs')

pgs <- list()
for(pop_i in keep_list$POP){
  pgs[[pop_i]] <- list()
  for(gwas_i in gwas_list$name){
    pgs[[pop_i]][[gwas_i]] <- list()
    for(pgs_method_i in pgs_methods){
      if(gwas_list$population[gwas_list$name == gwas_i] == 'EUR' | (gwas_list$population[gwas_list$name == gwas_i] != 'EUR' & (pgs_method_i %in% pgs_methods_eur))){
        pgs[[pop_i]][[gwas_i]][[pgs_method_i]] <- fread(paste0('/users/k1806347/oliverpainfel/Data/OpenSNP/GenoPred/test1/opensnp/pgs/',pop_i,'/',pgs_method_i,'/',gwas_i,'/opensnp-',gwas_i,'-',pop_i,'.profiles'))
      }
    }
  }
}

# Estimate correlation between pheno and pgs
assoc_all<-NULL
assoc <- list()
for(pop_i in keep_list$POP){
  assoc[[pop_i]] <- list()
  for(gwas_i in gwas_list$name){
    assoc[[pop_i]][[gwas_i]] <- list()
    for(pgs_method_i in pgs_methods){
      if(gwas_list$population[gwas_list$name == gwas_i] == 'EUR' | (gwas_list$population[gwas_list$name == gwas_i] != 'EUR' & (pgs_method_i %in% pgs_methods_eur))){
        pgs_i <- pgs[[pop_i]][[gwas_i]][[pgs_method_i]]
        pheno_pgs<-merge(pheno, pgs_i, by = c('FID','IID'))
        pheno_pgs$FID <- NULL
        pheno_pgs$IID <- NULL
        res <- rcorr(as.matrix(pheno_pgs))
        tmp <- 
          data.table(pop = pop_i,
                     gwas = gwas_i,
                     pgs_method = pgs_method_i,
                     name = colnames(res$r)[-1],
                     cor = as.numeric(res$r[1,-1]),
                     p = as.numeric(res$P[1,-1]),
                     n = as.numeric(res$n[1,-1])
          )
        assoc_all <- rbind(assoc_all, tmp[which(tmp$cor == max(tmp$cor, na.rm = T)),])
        assoc[[pop_i]][[gwas_i]][[pgs_method_i]] <- tmp[which(tmp$cor == max(tmp$cor, na.rm = T)),]
      }
    }
  }
}

# Its looking good. 
# Make a plot showing the R2 of the best hyperparameter and pseudovalidation/infintesimal models.
# Next steps.
# - Check these values match the expected
# - Put opensnp data through previous version of GenoPred for comparison
# - Could think about combining the opensnp datasets for increased precision, or use as external validation of results.
# 

q()
n

#####################
# Run opensnp data through the previous version of GenoPred for comparison
#####################

# Go to another version of the repo on CREATE
cd /users/k1806347/oliverpainfel/test/GenoPred

# Checkout to the v1 version of the repo
git checkout v1

####
# Prepare input files
####

R
library(data.table)
setwd('/users/k1806347/oliverpainfel/test/GenoPred/GenoPredPipe')

# target_list
target_list <- fread('target_list_example.txt')
target_list <- data.frame(
  name = 'opensnp',
  path='/users/k1806347/oliverpainfel/Data/OpenSNP/processed/geno/opensnp_train',
  type = 'samp_imp_vcf',
  output = '/users/k1806347/oliverpainfel/Data/OpenSNP/GenoPred_v1/test1',
  indiv_report = F)

dir.create('misc/opensnp', recursive = T)
write.table(target_list, 'misc/opensnp/target_list.txt', col.names=T, row.names=F, quote=F, sep=' ')

# gwas_list
yengo_eur <- fread('/users/k1806347/oliverpainfel/Data/GWAS_sumstats/opensnp_test/yengo_2022_height_eur.txt')
yengo_eur <- yengo_eur[, c('variant_id','effect_allele','other_allele','beta','standard_error','effect_allele_frequency','p_value','n'), with=F]
names(yengo_eur) <- c('SNP','A1','A2','BETA','SE','FREQ','P','N')
fwrite(yengo_eur, '/users/k1806347/oliverpainfel/Data/GWAS_sumstats/opensnp_test/yengo_2022_height_eur.format.txt', sep = ' ', quote = F, na = 'NA')

yengo_eas <- fread('/users/k1806347/oliverpainfel/Data/GWAS_sumstats/opensnp_test/yengo_2022_height_eas.txt')
yengo_eas <- yengo_eas[, c('variant_id','effect_allele','other_allele','beta','standard_error','effect_allele_frequency','p_value','n'), with=F]
names(yengo_eas) <- c('SNP','A1','A2','BETA','SE','FREQ','P','N')
fwrite(yengo_eas, '/users/k1806347/oliverpainfel/Data/GWAS_sumstats/opensnp_test/yengo_2022_height_eas.format.txt', sep = ' ', quote = F, na = 'NA')

gwas_list <- fread('gwas_list_example.txt')
gwas_list<-rbind(gwas_list, data.table(
    name='yengoeur',
    path='/users/k1806347/oliverpainfel/Data/GWAS_sumstats/opensnp_test/yengo_2022_height_eur.format.txt',
    population='EUR',
    sampling=NA,
    prevalence=NA,
    mean=NA,
    sd=NA,
    label="\"Yengo 2022 Height EUR\""))

gwas_list<-rbind(gwas_list, data.table(
    name='yengoeas',
    path='/users/k1806347/oliverpainfel/Data/GWAS_sumstats/opensnp_test/yengo_2022_height_eas.format.txt',
    population='EAS',
    sampling=NA,
    prevalence=NA,
    mean=NA,
    sd=NA,
    label="\"Yengo 2022 Height EAS\""))

gwas_list<-gwas_list[gwas_list$name %in% c('yengoeur','yengoeas'),]

write.table(gwas_list, 'misc/opensnp/gwas_list.txt', col.names=T, row.names=F, quote=F, sep=' ')

# score_list
score_list <- fread('score_list_example.txt')
score_list <- score_list[-1,]
write.table(score_list, 'misc/opensnp/score_list.txt', col.names=T, row.names=F, quote=F, sep=' ')

# config
config <- c(
  "gwas_list: misc/opensnp/gwas_list.txt",
  "target_list: misc/opensnp/target_list.txt",
  "score_list: misc/opensnp/score_list.txt"   
)
write.table(config, 'misc/opensnp/config.yaml', col.names=F, row.names=F, quote=F, sep=' ')

q()
n

# Run GenoPredPipe v1
snakemake -n --profile slurm --configfile=misc/opensnp/config.yaml --use-conda run_create_reports

# Calculate score using sbayesr, lassosum, ldpred2, prscs, and megaprs
snakemake -n --profile slurm --configfile=misc/opensnp/config.yaml --use-conda run_target_prs_all


#########
# Test correlation between PGS and phenotype
#########

setwd('/scratch//prj/oliverpainfel/test/GenoPred/GenoPredPipe/')
library(data.table)
library(Hmisc)

# Read in pheno data
pheno <- fread('/users/k1806347/oliverpainfel/Data/OpenSNP/processed/pheno/height.txt')
pheno$FID<-0

# Read in ancestry data
keep_list <- fread('/users/k1806347/oliverpainfel/Data/OpenSNP/GenoPred/test1/opensnp/ancestry/keep_list.txt')

# Read in pgs
gwas_list <- fread('misc/opensnp/gwas_list.txt')
# pgs_methods <- c('pt_clump','dbslmm','prscs','sbayesr','lassosum','ldpred2','megaprs')
pgs_methods <- c('pt_clump','dbslmm')
pgs_methods_eur <- c('pt_clump','dbslmm','lassosum','megaprs')

pgs <- list()
for(pop_i in keep_list$POP){
  pgs[[pop_i]] <- list()
  for(gwas_i in gwas_list$name){
    pgs[[pop_i]][[gwas_i]] <- list()
    for(pgs_method_i in pgs_methods){
      if(gwas_list$population[gwas_list$name == gwas_i] == 'EUR' | (gwas_list$population[gwas_list$name == gwas_i] != 'EUR' & (pgs_method_i %in% pgs_methods_eur))){
        pgs[[pop_i]][[gwas_i]][[pgs_method_i]] <- fread(paste0('/users/k1806347/oliverpainfel/Data/OpenSNP/GenoPred_v1/test1/opensnp/prs/',pop_i,'/',pgs_method_i,'/',gwas_i,'/opensnp.',gwas_i,'.',pop_i,'.profiles'))
      }
    }
  }
}

# Estimate correlation between pheno and pgs
assoc_all<-NULL
assoc <- list()
for(pop_i in keep_list$POP){
  assoc[[pop_i]] <- list()
  for(gwas_i in gwas_list$name){
    assoc[[pop_i]][[gwas_i]] <- list()
    for(pgs_method_i in pgs_methods){
      if(gwas_list$population[gwas_list$name == gwas_i] == 'EUR' | (gwas_list$population[gwas_list$name == gwas_i] != 'EUR' & (pgs_method_i %in% pgs_methods_eur))){
        pgs_i <- pgs[[pop_i]][[gwas_i]][[pgs_method_i]]
        pheno_pgs<-merge(pheno, pgs_i, by = c('FID','IID'))
        pheno_pgs$FID <- NULL
        pheno_pgs$IID <- NULL
        res <- rcorr(as.matrix(pheno_pgs))
        tmp <- 
          data.table(pop = pop_i,
                     gwas = gwas_i,
                     pgs_method = pgs_method_i,
                     name = colnames(res$r)[-1],
                     cor = as.numeric(res$r[1,-1]),
                     p = as.numeric(res$P[1,-1]),
                     n = as.numeric(res$n[1,-1])
          )
        assoc_all <- rbind(assoc_all, tmp[which(tmp$cor == max(tmp$cor, na.rm = T)),])
        assoc[[pop_i]][[gwas_i]][[pgs_method_i]] <- tmp[which(tmp$cor == max(tmp$cor, na.rm = T)),]
      }
    }
  }
}

# Results are consistent. Hard to compare exactly due to different ancestry_cut off.
# Check correlation when using more stringent ancestry cut off
# Read in ancestry predictions
model_pred <- fread('/users/k1806347/oliverpainfel/Data/OpenSNP/GenoPred_v1/test1/opensnp/ancestry/ancestry_all/opensnp.Ancestry.model_pred')

pgs_strict <- list()
for(pop_i in keep_list$POP){
  pop_i_keep <- model_pred[model_pred[[pop_i]] > 0.95, ]
  pgs_strict[[pop_i]] <- list()
  for(gwas_i in gwas_list$name){
    pgs_strict[[pop_i]][[gwas_i]] <- list()
    for(pgs_method_i in pgs_methods){
      if(gwas_list$population[gwas_list$name == gwas_i] == 'EUR' | (gwas_list$population[gwas_list$name == gwas_i] != 'EUR' & (pgs_method_i %in% pgs_methods_eur))){
        pgs_strict[[pop_i]][[gwas_i]][[pgs_method_i]] <- merge(pgs[[pop_i]][[gwas_i]][[pgs_method_i]], pop_i_keep[, c('FID','IID'), with = F], by = c('FID','IID'))
      }
    }
  }
}

assoc_all_strict<-NULL
assoc_strict <- list()
for(pop_i in keep_list$POP){
  assoc_strict[[pop_i]] <- list()
  for(gwas_i in gwas_list$name){
    assoc_strict[[pop_i]][[gwas_i]] <- list()
    for(pgs_method_i in pgs_methods){
      if(gwas_list$population[gwas_list$name == gwas_i] == 'EUR' | (gwas_list$population[gwas_list$name == gwas_i] != 'EUR' & (pgs_method_i %in% pgs_methods_eur))){
        pgs_i <- pgs_strict[[pop_i]][[gwas_i]][[pgs_method_i]]
        pheno_pgs<-merge(pheno, pgs_i, by = c('FID','IID'))
        pheno_pgs$FID <- NULL
        pheno_pgs$IID <- NULL
        res <- rcorr(as.matrix(pheno_pgs))
        tmp <- 
          data.table(pop = pop_i,
                     gwas = gwas_i,
                     pgs_method = pgs_method_i,
                     name = colnames(res$r)[-1],
                     cor = as.numeric(res$r[1,-1]),
                     p = as.numeric(res$P[1,-1]),
                     n = as.numeric(res$n[1,-1])
          )
        assoc_all_strict <- rbind(assoc_all_strict, tmp[which(tmp$cor == max(tmp$cor, na.rm = T)),])
        assoc_strict[[pop_i]][[gwas_i]][[pgs_method_i]] <- tmp[which(tmp$cor == max(tmp$cor, na.rm = T)),]
      }
    }
  }
}

# Results are highly concordant

q()
n




