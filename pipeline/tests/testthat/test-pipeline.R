library(data.table)
library(testthat)

tempdir <- function(prefix = "tmpdir") {
  tmpdir <- tempfile(pattern = prefix)
  dir.create(tmpdir)
  return(tmpdir)
}

# Create temp_dir
temp_dir<-tempdir()

################
# Prepare configuration files for testing
################

#######
# Test 1: Run pipeline with standard set up
#######

# Update config to use temp_dir as outdir and resdir
# Note. Using read_yaml and write_yaml converts single line lists to multiline lists, which causes an error in snakemake
config <- readLines('misc/dev/test_data/config/config.yaml')
config[grepl('^resdir', config)]<-paste0('resdir: ', temp_dir)
config[grepl('^outdir', config)]<-paste0('outdir: ', temp_dir)
config[grepl('^config_file', config)]<-paste0('config_file: ', temp_dir, '/config.yaml')
write.table(config, paste0(temp_dir, '/config1.yaml'), col.names = F, row.names = F, quote = F)

#######
# Test 2: Update gwas_list
#######

# Add new gwas to gwas_list
gwas_list <- fread(gsub('gwas_list: ','', config[grepl('^gwas_list', config)]))
gwas_list <- rbind(gwas_list, gwas_list)
gwas_list$name[2] <- 'new_BODY04'
gwas_list$label<-paste0("\"", gwas_list$label, "\"")
write.table(gwas_list, paste0(temp_dir, '/gwas_list.txt'), col.names=T, row.names=F, quote=F)

# Update config file
config_tmp<-config
config_tmp[grepl('^gwas_list', config_tmp)]<-paste0('gwas_list: ', temp_dir, '/gwas_list.txt')
write.table(config_tmp, paste0(temp_dir, '/config2.yaml'), col.names = F, row.names = F, quote = F)

#######
# Test 3: Update target_list
#######

# Add new target to target_list
target_list <- fread(gsub('target_list: ','', config[grepl('^target_list', config)]))
target_list <- rbind(target_list, target_list)
target_list$name[2] <- 'new_example_plink2'
write.table(target_list, paste0(temp_dir, '/target_list.txt'), col.names=T, row.names=F, quote=F)

# Update config file
config_tmp<-config
config_tmp[grepl('^target_list', config_tmp)]<-paste0('target_list: ', temp_dir, '/target_list.txt')
write.table(config_tmp, paste0(temp_dir, '/config3.yaml'), col.names = F, row.names = F, quote = F)

#######
# Test 4: Update score_list
#######

# Add new target to target_list
score_list <- fread(gsub('score_list: ','', config[grepl('^score_list', config)]))
score_list <- rbind(score_list, score_list)
score_list$name[2] <- 'new_PGS002804'
score_list$label<-paste0("\"", score_list$label, "\"")
write.table(score_list, paste0(temp_dir, '/score_list.txt'), col.names=T, row.names=F, quote=F)

# Update config file
config_tmp<-config
config_tmp[grepl('^score_list', config_tmp)]<-paste0('score_list: ', temp_dir, '/score_list.txt')
write.table(config_tmp, paste0(temp_dir, '/config4.yaml'), col.names = F, row.names = F, quote = F)

#######
# Test 5: Update pgs_methods
#######

# Add lassosum to pgs_methods list
config_tmp<-config
config_tmp[grepl('^pgs_methods', config_tmp)]<-paste0("pgs_methods: ['ptclump','lassosum']")
write.table(config_tmp, paste0(temp_dir, '/config5.yaml'), col.names = F, row.names = F, quote = F)

#################
# Run pipeline commands inside container
#################

system(paste0(
  "singularity exec --writable-tmpfs /users/k1806347/oliverpainfel/Software/singularity/genopred_pipeline_latest.sif bash -c \"
      source /opt/mambaforge/etc/profile.d/conda.sh &&
      conda activate genopred &&
      cd /tools/GenoPred/pipeline &&
      # Clear output directory
      rm -r ", temp_dir, "/*
      # Test 1
      cp ", temp_dir, "/config1.yaml ", temp_dir, "/config.yaml &&
      snakemake -j1 --use-conda output_all --configfile=", temp_dir, "/config.yaml > ", temp_dir, "/snakemake1.log 2>&1 &&
      # Test 2
      cp ", temp_dir, "/config2.yaml ", temp_dir, "/config.yaml &&
      snakemake -j1 -n --use-conda output_all --configfile=", temp_dir, "/config.yaml > ", temp_dir, "/snakemake2.log 2>&1 &&
      # Test 3
      cp ", temp_dir, "/config3.yaml ", temp_dir, "/config.yaml &&
      snakemake -j1 -n --use-conda output_all --configfile=", temp_dir, "/config.yaml > ", temp_dir, "/snakemake3.log 2>&1 &&
      # Test 4
      cp ", temp_dir, "/config4.yaml ", temp_dir, "/config.yaml &&
      snakemake -j1 -n --use-conda output_all --configfile=", temp_dir, "/config.yaml > ", temp_dir, "/snakemake4.log 2>&1
      # Test 5
      cp ", temp_dir, "/config5.yaml ", temp_dir, "/config.yaml &&
      snakemake -j1 -n --use-conda output_all --configfile=", temp_dir, "/config.yaml > ", temp_dir, "/snakemake5.log 2>&1
    \""
))

#################
# Perform tests
#################

######
# Test 1
######

test_that("Run pipeline with standard config", {

  # Read in log file
  log<-readLines(paste0(temp_dir, '/snakemake1.log'))

  # Check it finished without errors
  expect_true(any(grepl("steps \\(100%\\) done", log)))

  ###############
  # Check the output of each step
  ###############

  #######
  # target_qc
  #######

  ###
  # ancestry_inference_i
  ###
  results<-fread(paste0(temp_dir,'/example_plink2/ancestry/example_plink2.Ancestry.model_pred'))
  expected<-fread('tests/expected/output/example_plink2/ancestry/example_plink2.Ancestry.model_pred')
  expect_equal(expected, results)

  ###
  # ancestry_reporter
  ###
  results<-fread(paste0(temp_dir,'/example_plink2/ancestry/ancestry_report.txt'))
  expected<-fread('tests/expected/output/example_plink2/ancestry/ancestry_report.txt')
  expect_equal(expected, results)

  ###
  # format_target_i
  ###
  for(i in c('pgen','psam','pvar')){
    results<-paste0(temp_dir,'/example_plink2/geno/example_plink2.ref.chr22.', i)
    expected<-paste0('tests/expected/output/example_plink2/geno/example_plink2.ref.chr22.', i)

    diff_result <- system2("diff", args = c("-q", expected, results), stdout = TRUE, stderr = TRUE)
    expect_true(identical(diff_result, character(0)))
  }

  #######
  # pgs_methods
  #######

  ###
  # sumstat_prep_i
  ###
  results<-fread(paste0(temp_dir,'/reference/gwas_sumstat/BODY04/BODY04-cleaned.gz'))
  expected<-fread('tests/expected/output/reference/gwas_sumstat/BODY04/BODY04-cleaned.gz')
  expect_equal(expected, results)

  ###
  # prep_pgs_ptclump_i
  ###
  # score.gz
  results <- fread(paste0(temp_dir,'/reference/pgs_score_files/ptclump/BODY04/ref-BODY04.score.gz'))
  expected<-fread('tests/expected/output/reference/pgs_score_files/ptclump/BODY04/ref-BODY04.score.gz')
  expect_equal(expected, results)

  # AFR.scale
  results<-fread(paste0(temp_dir,'/reference/pgs_score_files/ptclump/BODY04/ref-BODY04-AFR.scale'))
  expected<-fread('tests/expected/output/reference/pgs_score_files/ptclump/BODY04/ref-BODY04-AFR.scale')
  expect_equal(expected, results)

  ###
  # prep_pgs_external_i
  ###
  # score.gz
  results<-fread(paste0(temp_dir,'/reference/pgs_score_files/external/PGS002804/ref-PGS002804.score.gz'))
  expected<-fread('tests/expected/output/reference/pgs_score_files/external/PGS002804/ref-PGS002804.score.gz')
  expect_equal(expected, results)

  # AFR.scale
  results<-fread(paste0(temp_dir,'/reference/pgs_score_files/external/PGS002804/ref-PGS002804-AFR.scale'))
  expected<-fread('tests/expected/output/reference/pgs_score_files/external/PGS002804/ref-PGS002804-AFR.scale')
  expect_equal(expected, results)

  #######
  # target_scoring
  #######

  ###
  # target_pgs_i
  ###

  # external
  results <- fread(paste0(temp_dir,'/example_plink2/pgs/AFR/external/PGS002804/example_plink2-PGS002804-AFR.profiles'))
  expected <- fread('tests/expected/output/example_plink2/pgs/AFR/external/PGS002804/example_plink2-PGS002804-AFR.profiles')
  expect_equal(expected, results)

  # ptclump
  results <- fread(paste0(temp_dir,'/example_plink2/pgs/AFR/ptclump/BODY04/example_plink2-BODY04-AFR.profiles'))
  expected <- fread('tests/expected/output/example_plink2/pgs/AFR/ptclump/BODY04/example_plink2-BODY04-AFR.profiles')
  expect_equal(expected, results)

  #######
  # report
  #######

  ###
  # sample_report_i
  ###
  results<-paste0(temp_dir,'/example_plink2/reports/example_plink2-report.html')
  expect_true(file.exists(results))

})

######
# Test 2
######

# Modify gwas_list and check jobs to be rerun
test_that("Modify gwas_list and check jobs to be rerun", {

  # Read in log file
  log<-readLines(paste0(temp_dir, '/snakemake2.log'))

  # List rules to be run
  to_do<-gsub(' .*','', log[(grep('^Job stats:', log)[2]+3):(grep('^Reasons:', log)-3)])

  # Check additional output from sumstat_prep_i, prep_pgs_ptclump_i etc. are required
  expect_true(all(c("output_all", "prep_pgs_ptclump_i", "sample_report_i", "sumstat_prep_i", "target_pgs_all", "target_pgs_i") %in% to_do), log)
})

######
# Test 3
######

# Modify gwas_list and check jobs to be rerun
test_that("Modify target_list and check jobs to be rerun", {

  # Read in log file
  log<-readLines(paste0(temp_dir, '/snakemake3.log'))

  # List rules to be run
  to_do<-gsub(' .*','', log[(grep('^Job stats:', log)[2]+3):(grep('^Reasons:', log)-3)])

  # Check additional output from format_target_i etc. are required
  expect_true(all(c("ancestry_inference_i", "ancestry_reporter", "format_target_all", "format_target_i", "output_all", "sample_report_i", "target_pgs_all") %in% to_do), log)
})

######
# Test 4
######

# Modify gwas_list and check jobs to be rerun
test_that("Modify score_list and check jobs to be rerun", {

  # Read in log file
  log<-readLines(paste0(temp_dir, '/snakemake4.log'))

  # List rules to be run
  to_do<-gsub(' .*','', log[(grep('^Job stats:', log)[2]+3):(grep('^Reasons:', log)-3)])

  # Check additional output from format_target_i etc. are required
  expect_true(all(c("output_all", "prep_pgs_external_i", "sample_report_i", "score_reporter", "target_pgs_all", "target_pgs_i") %in% to_do), log)
})
