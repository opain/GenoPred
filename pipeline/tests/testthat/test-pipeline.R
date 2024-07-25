library(data.table)
library(testthat)

# Change directory
setwd('../../')

tempdir <- function(prefix = "tmpdir") {
  tmpdir <- tempfile(pattern = prefix)
  dir.create(tmpdir)
  return(tmpdir)
}

# Create temp_dir
temp_dir<-tempdir()

# Retrieve .SIF file location
sif_file <- Sys.getenv("SIF_FILE", '/users/k1806347/oliverpainfel/Software/singularity/genopred_pipeline_latest.sif')

if(sif_file == ''){
  stop("SIF_FILE is not specified")
}
if(!file.exists(sif_file)){
  stop("SIF_FILE does not exist")
}

# Retrieve branch name
repo_branch <- Sys.getenv("BRANCH_NAME", 'dev')

if(repo_branch == ''){
  stop("BRANCH_NAME is not specified")
}

################
# Prepare configuration files for testing
################

#######
# Test 1: Run pipeline with standard set up
#######

# Update config to use temp_dir as outdir and resdir
# Note. Using read_yaml and write_yaml converts single line lists to multiline lists, which causes an error in snakemake
config <- readLines('misc/dev/test_data/config/config.yaml')
config[grepl('^resdir', config)]<-paste0('resdir: ', temp_dir, '/resources')
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
config_tmp[grepl('^pgs_methods', config_tmp)]<-paste0("pgs_methods: ['ptclump','lassosum','dbslmm']")
write.table(config_tmp, paste0(temp_dir, '/config5.yaml'), col.names = F, row.names = F, quote = F)

#######
# Test 6: No target_list specified
#######

# Create new temp_dir so the pipeline starts a fresh
temp_dir2<-tempdir()

config <- readLines('misc/dev/test_data/config/config.yaml')
config[grepl('^resdir', config)]<-paste0('resdir: ', temp_dir2, '/resources')
config[grepl('^outdir', config)]<-paste0('outdir: ', temp_dir2)
config[grepl('^config_file', config)]<-paste0('config_file: ', temp_dir2, '/config6.yaml')

# Remove target_list
config<-config[!grepl('^target_list', config)]

write.table(config, paste0(temp_dir2, '/config6.yaml'), col.names = F, row.names = F, quote = F)

#######
# Test 7: No score_list specified
#######

config <- readLines('misc/dev/test_data/config/config.yaml')
config[grepl('^resdir', config)]<-paste0('resdir: ', temp_dir2, '/resources')
config[grepl('^outdir', config)]<-paste0('outdir: ', temp_dir2)
config[grepl('^config_file', config)]<-paste0('config_file: ', temp_dir2, '/config6.yaml')

# Remove target_list
config<-config[!grepl('^score_list', config)]

write.table(config, paste0(temp_dir2, '/config7.yaml'), col.names = F, row.names = F, quote = F)

#######
# Test 8: No gwas_list specified
#######

config <- readLines('misc/dev/test_data/config/config.yaml')
config[grepl('^resdir', config)]<-paste0('resdir: ', temp_dir2, '/resources')
config[grepl('^outdir', config)]<-paste0('outdir: ', temp_dir2)
config[grepl('^config_file', config)]<-paste0('config_file: ', temp_dir2, '/config8.yaml')

# Remove target_list
config<-config[!grepl('^gwas_list', config)]

write.table(config, paste0(temp_dir2, '/config8.yaml'), col.names = F, row.names = F, quote = F)

#######
# Test 9: No gwas_list or score_list specified
#######

config <- readLines('misc/dev/test_data/config/config.yaml')
config[grepl('^resdir', config)]<-paste0('resdir: ', temp_dir2, '/resources')
config[grepl('^outdir', config)]<-paste0('outdir: ', temp_dir2)
config[grepl('^config_file', config)]<-paste0('config_file: ', temp_dir2, '/config9.yaml')

# Remove target_list
config<-config[!grepl('^gwas_list', config)]
config<-config[!grepl('^score_list', config)]

write.table(config, paste0(temp_dir2, '/config9.yaml'), col.names = F, row.names = F, quote = F)

#######
# Future tests to be implemented
#######
# other pgs methods - It will take a lot more time to run due to download of reference data and slower PGS methods.

#################
# Run pipeline commands inside container
#################

requested_output <- paste0("output_all pc_projection ", temp_dir, "/reference/target_checks/example_plink2/indiv_report-4_EAS.4_EAS-report.done")

exit_status <- system(paste0(
  "singularity exec --writable-tmpfs ", sif_file, " bash -c \"
      # Set to exit if any errors incurred
      set -e &&
      # Initiate conda
      source /opt/mambaforge/etc/profile.d/conda.sh &&
      # Activate genopred environment
      conda activate genopred &&
      # Go to repo
      cd /tools/GenoPred/pipeline &&
      # Checkout to dev branch
      git checkout ", repo_branch, " &&
      # Fetch latest commits
      git pull &&
      # Clear output directory
      rm -r -f ", temp_dir, "/reference &&
      rm -r -f ", temp_dir, "/example_plink2 &&
      rm -r -f ", temp_dir, "/resources &&
      # Test 1
      cp ", temp_dir, "/config1.yaml ", temp_dir, "/config.yaml &&
      snakemake -j1 --use-conda ", requested_output, " --configfile=", temp_dir, "/config.yaml > ", temp_dir, "/snakemake1.log 2>&1 &&
      # Test 2
      cp ", temp_dir, "/config2.yaml ", temp_dir, "/config.yaml &&
      snakemake -j1 -n --use-conda ", requested_output, " --configfile=", temp_dir, "/config.yaml > ", temp_dir, "/snakemake2.log 2>&1 &&
      # Test 3
      cp ", temp_dir, "/config3.yaml ", temp_dir, "/config.yaml &&
      snakemake -j1 -n --use-conda ", requested_output, " --configfile=", temp_dir, "/config.yaml > ", temp_dir, "/snakemake3.log 2>&1 &&
      # Test 4
      cp ", temp_dir, "/config4.yaml ", temp_dir, "/config.yaml &&
      snakemake -j1 -n --use-conda ", requested_output, " --configfile=", temp_dir, "/config.yaml > ", temp_dir, "/snakemake4.log 2>&1 &&
      # Test 5
      cp ", temp_dir, "/config5.yaml ", temp_dir, "/config.yaml &&
      snakemake -j1 -n --use-conda ", requested_output, " --configfile=", temp_dir, "/config.yaml > ", temp_dir, "/snakemake5.log 2>&1 &&
      # Test 6
      snakemake -j1 -n --use-conda prep_pgs --configfile=", temp_dir2, "/config6.yaml > ", temp_dir2, "/snakemake6.log 2>&1 &&
      # Test 7
      snakemake -j1 -n --use-conda output_all --configfile=", temp_dir2, "/config7.yaml > ", temp_dir2, "/snakemake7.log 2>&1 &&
      # Test 8
      snakemake -j1 -n --use-conda output_all --configfile=", temp_dir2, "/config8.yaml > ", temp_dir2, "/snakemake8.log 2>&1 &&
      # Test 9
      snakemake -j1 -n --use-conda output_all --configfile=", temp_dir2, "/config9.yaml > ", temp_dir2, "/snakemake9.log 2>&1
    \""
))

# Check system command completed successfully
if (exit_status != 0) {
  stop("The script failed. Check the logs for details.")
}

#################
# Perform tests
#################

######
# Test 1
######

test_that("Check pipeline runs without error", {

  # Read in log file
  log<-readLines(paste0(temp_dir, '/snakemake1.log'))

  # Check it finished without errors
  expect_true(any(grepl("steps \\(100%\\) done", log)))

})

###############
# Check the output of each step
###############

#######
# target_qc
#######

###
# ancestry_inference_i
###

test_that("Check ancestry_inference_i output", {
  results<-fread(paste0(temp_dir,'/example_plink2/ancestry/example_plink2.Ancestry.model_pred'))
  expected<-fread('misc/dev/test_data/output/example_plink2/ancestry/example_plink2.Ancestry.model_pred')
  expect_equal(expected, results)
})

###
# ancestry_reporter
###

test_that("Check ancestry_reporter output", {
  results<-fread(paste0(temp_dir,'/example_plink2/ancestry/ancestry_report.txt'))
  expected<-fread('misc/dev/test_data/output/example_plink2/ancestry/ancestry_report.txt')
  expect_equal(expected, results)
})

###
# format_target_i
###

for(i in c('pgen','psam','pvar')){
  test_that(paste0("Check format_target_i output: ", i), {
    results<-paste0(temp_dir,'/example_plink2/geno/example_plink2.ref.chr22.', i)
    expected<-paste0('misc/dev/test_data/output/example_plink2/geno/example_plink2.ref.chr22.', i)

    diff_result <- system2("diff", args = c("-q", expected, results), stdout = TRUE, stderr = TRUE)
    expect_true(identical(diff_result, character(0)))
  })
}

#######
# pgs_methods
#######

###
# ref_pca_i
###

for(i in c('eigenvec.var.gz','AFR.scale')){
  test_that(paste0("Check ref_pca_i output: ", i), {
    results<-fread(paste0(temp_dir,'/resources/data/ref/pc_score_files/AFR/ref-AFR-pcs.', i))
    expected<-fread(paste0('misc/dev/test_data/output/resources/data/ref/pc_score_files/AFR/ref-AFR-pcs.', i))
    expect_equal(expected, results)
  })
}

###
# sumstat_prep_i
###

test_that("Check sumstat_prep_i output", {
  results<-fread(paste0(temp_dir,'/reference/gwas_sumstat/BODY04/BODY04-cleaned.gz'))
  expected<-fread('misc/dev/test_data/output/reference/gwas_sumstat/BODY04/BODY04-cleaned.gz')
  expect_equal(expected, results)
})

###
# internal
###
for(i in c('ptclump','lassosum')){
  for(j in c('.score.gz','-AFR.scale')){
    test_that(paste0("Check prep_pgs_", i,"_i output: ", j), {
      results <- fread(paste0(temp_dir,'/reference/pgs_score_files/', i,'/BODY04/ref-BODY04', j))
      expected<-fread(paste0('misc/dev/test_data/output/reference/pgs_score_files/', i, '/BODY04/ref-BODY04', j))
      expect_equal(expected, results)
    })
  }
}

###
# prep_pgs_lassosum_i
###
for(i in c('.score.gz','-AFR.scale')){
  test_that(paste0("Check prep_pgs_ptclump_i output: ", i), {
    results <- fread(paste0(temp_dir,'/reference/pgs_score_files/ptclump/BODY04/ref-BODY04', i))
    expected<-fread(paste0('misc/dev/test_data/output/reference/pgs_score_files/ptclump/BODY04/ref-BODY04', i))
    expect_equal(expected, results)
  })
}

###
# prep_pgs_external_i
###
for(i in c('.score.gz','-AFR.scale')){
  test_that(paste0("Check prep_pgs_external_i output: ", i), {
    results<-fread(paste0(temp_dir,'/reference/pgs_score_files/external/PGS002804/ref-PGS002804', i))
    expected<-fread(paste0('misc/dev/test_data/output/reference/pgs_score_files/external/PGS002804/ref-PGS002804', i))
    expect_equal(expected, results)
  })
}

#######
# target_scoring
#######

###
# target_pgs_i
###

# external
test_that("Check target_pgs_i output: external", {
  results <- fread(paste0(temp_dir,'/example_plink2/pgs/AFR/external/PGS002804/example_plink2-PGS002804-AFR.profiles'))
  expected <- fread('misc/dev/test_data/output/example_plink2/pgs/AFR/external/PGS002804/example_plink2-PGS002804-AFR.profiles')
  expect_equal(expected, results)
})

# internal
for(i in c('ptclump','lassosum')){
  test_that(paste0("Check target_pgs_i output: ", i), {
    results <- fread(paste0(temp_dir,'/example_plink2/pgs/AFR/', i, '/BODY04/example_plink2-BODY04-AFR.profiles'))
    expected <- fread(paste0('misc/dev/test_data/output/example_plink2/pgs/AFR/', i, '/BODY04/example_plink2-BODY04-AFR.profiles'))
    expect_equal(expected, results)
  })
}

#######
# report
#######

###
# sample_report_i
###
test_that("Check sample_report_i output", {
  results<-paste0(temp_dir,'/example_plink2/reports/example_plink2-report.html')
  expect_true(file.exists(results))
})

###
# indiv_report_i
###
test_that("Check indiv_report_i output", {
  results<-paste0(temp_dir,'/example_plink2/reports/individual/example_plink2-4_EAS.4_EAS-report.html')
  expect_true(file.exists(results))
})

######
# Test 2
######

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

test_that("Modify score_list and check jobs to be rerun", {

  # Read in log file
  log<-readLines(paste0(temp_dir, '/snakemake4.log'))

  # List rules to be run
  to_do<-gsub(' .*','', log[(grep('^Job stats:', log)[2]+3):(grep('^Reasons:', log)-3)])

  # Check additional output from format_target_i etc. are required
  expect_true(all(c("output_all", "prep_pgs_external_i", "sample_report_i", "score_reporter", "target_pgs_all", "target_pgs_i") %in% to_do), log)
})

######
# Test 5
######

test_that("Modify pgs_methods and check jobs to be rerun", {

  # Read in log file
  log<-readLines(paste0(temp_dir, '/snakemake5.log'))

  # List rules to be run
  to_do<-gsub(' .*','', log[(grep('^Job stats:', log)[2]+3):(grep('^Reasons:', log)-3)])

  # Check additional output from format_target_i etc. are required
  expect_true(all(c("download_dbslmm", "download_hm3_snplist", "download_ld_blocks", "download_ldsc", "download_ldscores_panukb", "download_plink", "indiv_report_i", "output_all", "prep_pgs_dbslmm_i", "sample_report_i", "target_pgs_all", "target_pgs_i") %in% to_do), log)
})

######
# Test 6
######

test_that("Check jobs to be run when no target_list", {

  # Read in log file
  log<-readLines(paste0(temp_dir2, '/snakemake6.log'))

  # List rules to be run
  to_do<-gsub(' .*','', log[(grep('^Job stats:', log)[2]+3):(grep('^Reasons:', log)-3)])

  # Check additional output from format_target_i etc. are required
  expect_true(all(c("prep_pgs", "prep_pgs_external_i", "prep_pgs_lassosum_i", "prep_pgs_ptclump_i", "score_reporter", "sumstat_prep_i") %in% to_do), log)
})

######
# Test 7
######

test_that("Check jobs to be run when no score_list", {

  # Read in log file
  log<-readLines(paste0(temp_dir2, '/snakemake7.log'))

  # List rules to be run
  to_do<-gsub(' .*','', log[(grep('^Job stats:', log)[2]+3):(grep('^Reasons:', log)-3)])

  # Check additional output from format_target_i etc. are required
  expect_true(all(c("ancestry_inference_i", "ancestry_reporter", "format_target_all", "format_target_i", "output_all", "sample_report_i", "sumstat_prep_i", "target_pgs_all") %in% to_do), log)
})

######
# Test 8
######

test_that("Check jobs to be run when no gwas_list", {

  # Read in log file
  log<-readLines(paste0(temp_dir2, '/snakemake8.log'))

  # List rules to be run
  to_do<-gsub(' .*','', log[(grep('^Job stats:', log)[2]+3):(grep('^Reasons:', log)-3)])

  # Check additional output from format_target_i etc. are required
  expect_true(all(c("ancestry_inference_i", "ancestry_reporter", "format_target_all", "format_target_i", "output_all", "prep_pgs_external_i", "sample_report_i", "score_reporter", "target_pgs_all") %in% to_do), log)
})

######
# Test 9
######

test_that("Check jobs to be run when no score_list or gwas_list", {

  # Read in log file
  log<-readLines(paste0(temp_dir2, '/snakemake9.log'))

  # List rules to be run
  to_do<-gsub(' .*','', log[(grep('^Job stats:', log)[2]+3):(grep('^Reasons:', log)-3)])

  # Check additional output from format_target_i etc. are required
  expect_true(all(c("ancestry_inference_i", "ancestry_reporter", "format_target_all", "format_target_i", "output_all", "sample_report_i", "target_pgs_all") %in% to_do), log)
})

