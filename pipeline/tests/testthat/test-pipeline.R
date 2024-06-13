library(data.table)

tempdir <- function(prefix = "tmpdir") {
  tmpdir <- tempfile(pattern = prefix)
  dir.create(tmpdir)
  return(tmpdir)
}

# Build conda environments
test_that("Building conda environments", {
  # Create resdir
  resdir<-tempdir()

  # Update config to use resdir
  config<-readLines('misc/dev/test_data/config/config.yaml')
  config[grepl('^resdir', config)]<-paste0('resdir: ', resdir)
  write.table(config, paste0(resdir, '/config.yaml'), col.names = F, row.names = F, quote = F)

  # Run pipeline
  system(paste0(
    'tests/run_snakemake.sh --configfile=',resdir,'/config.yaml --restart-times 3 -j1 --use-conda --conda-frontend mamba install_r_packages resources/software/pgscatalog_utils/download_pgscatalog_utils.done > ',resdir,'/snakemake.log 2>&1'
  ))

  # Check it finished
  log<-readLines(paste0(resdir, '/snakemake.log'))
  expect_true(any(grepl("steps \\(100%\\) done", log)))
})

# Run pipeline using mini test data and compare outputs to the expected
test_that("Run pipeline with standard config", {
  # Create temp_dir
  temp_dir<-tempdir()

  # Update config to use temp_dir as outdir
  config<-readLines(paste0(resdir,'/config.yaml'))
  config[grepl('^outdir', config)]<-paste0('outdir: ', temp_dir)
  config[grepl('^config_file', config)]<-paste0('config_file: ', temp_dir, '/config.yaml')
  write.table(config, paste0(temp_dir, '/config.yaml'), col.names = F, row.names = F, quote = F)

  # Run pipeline
  system(paste0(
    'tests/run_snakemake.sh -j1 --use-conda output_all --configfile=',temp_dir,'/config.yaml > ',temp_dir,'/snakemake.log 2>&1'
  ))

  log<-readLines(paste0(temp_dir, '/snakemake.log'))

  # Check it finished without errors
  expect_true(any(grepl("steps \\(100%\\) done", log)))

  ###############
  # Check the output of each step
  ###############

  ###
  # ancestry_inference_i
  ###
  results<-read.table(paste0(temp_dir,'/example_plink2/ancestry/example_plink2.Ancestry.model_pred'), header = T, stringsAsFactors = F)
  expected<-read.table('tests/expected/output/example_plink2/ancestry/example_plink2.Ancestry.model_pred', header = T, stringsAsFactors = F)
  expect_equal(expected, results)

  ###
  # ancestry_reporter
  ###
  results<-read.table(paste0(temp_dir,'/example_plink2/ancestry/ancestry_report.txt'), header = T, stringsAsFactors = F)
  expected<-read.table('tests/expected/output/example_plink2/ancestry/ancestry_report.txt', header = T, stringsAsFactors = F)
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

  ###
  # prep_pgs_external_i
  ###
  # score.gz
  results<-read.table(paste0(temp_dir,'/reference/pgs_score_files/external/PGS002804/ref-PGS002804.score.gz'), header = T, stringsAsFactors = F)
  expected<-read.table('tests/expected/output/reference/pgs_score_files/external/PGS002804/ref-PGS002804.score.gz', header = T, stringsAsFactors = F)
  expect_equal(expected, results)

  # AFR.scale
  results<-read.table(paste0(temp_dir,'/reference/pgs_score_files/external/PGS002804/ref-PGS002804-AFR.scale'), header = T, stringsAsFactors = F)
  expected<-read.table('tests/expected/output/reference/pgs_score_files/external/PGS002804/ref-PGS002804-AFR.scale', header = T, stringsAsFactors = F)
  expect_equal(expected, results)

  ###
  # sample_report_i
  ###
  results<-paste0(temp_dir,'/example_plink2/reports/example_plink2-report.html')
  expect_true(file.exists(results))

  ###
  # sumstat_prep_i
  ###
  results<-read.table(paste0(temp_dir,'/reference/gwas_sumstat/BODY04/BODY04-cleaned.gz'), header = T, stringsAsFactors = F)
  expected<-read.table('tests/expected/output/reference/gwas_sumstat/BODY04/BODY04-cleaned.gz', header = T, stringsAsFactors = F)
  expect_equal(expected, results)
})
