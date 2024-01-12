##################
# Prepare input GenoPred when using 23andMe test data
##################

library(data.table)

####
# gwas_list
####
# We will use one of the GWAS provided in the test_data package
gwas_list <- fread('example_input/gwas_list_example.txt')
gwas_list <- gwas_list[1, ]

gwas_list$label <- paste0("\"", gwas_list$label, "\"")

write.table(gwas_list, 'misc/23andMe/gwas_list.txt', col.names=T, row.names=F, quote=F, sep=' ')

####
# target_list
####

target_list <- data.frame(
    name = 'Joe_Bloggs',
    path = 'test_data/target/23andMe_individual/Joe_Bloggs_genome_0123456789.zip',
    type = '23andMe',
    indiv_report = T
)

write.table(target_list, 'misc/23andMe/target_list.txt', col.names=T, row.names=F, quote=F, sep=' ')

####
# config
####

config <- c(
    'outdir: test_data/output/23andMe',
    'config_file: misc/23andMe/config.yaml',
    'gwas_list: misc/23andMe/gwas_list.txt',
    'target_list: misc/23andMe/target_list.txt',
    'score_list: NA',
    "pgs_methods: ['dbslmm']",
    'testing: chr22'
)

write.table(config, 'misc/23andMe/config.yaml', col.names=F, row.names=F, quote=F)
