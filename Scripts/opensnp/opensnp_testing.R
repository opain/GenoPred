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


