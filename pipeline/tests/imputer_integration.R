#!/usr/bin/env Rscript
# Run in the pipeline analysis R environment from the repository root.
root <- normalizePath('.')
work <- tempfile('imputer integration '); dir.create(work)
on.exit(unlink(work, recursive = TRUE), add = TRUE)
dir.create(file.path(work, 'run'))
stopifnot(file.symlink(file.path(root, 'functions'), file.path(work, 'functions')))
ref <- file.path(work, 'ref folder'); dir.create(ref)
legend <- gzfile(file.path(ref, 'panel_chr22.legend.gz'), 'wt')
writeLines(c('id position a0 a1', 'rs123 100 A C'), legend); close(legend)
genotype <- file.path(work, 'genotype $(unsafe).txt'); writeLines('example', genotype)
output <- file.path(work, 'output folder', 'result')
mock <- file.path(root, 'pipeline/tests/mock_imputation_tool.py')
# The test does not require the checkout to have executable test files.
executable <- file.path(work, 'mock tool')
stopifnot(file.copy(mock, executable)); Sys.chmod(executable, '0755')
script <- file.path(root, 'Scripts/23andMe_imputer/23andMe_imputer.R')
old <- setwd(file.path(work, 'run')); on.exit(setwd(old), add = TRUE)
args <- c(script, '--geno', genotype, '--ref', ref, '--ref_name', 'panel', '--ref5', ref,
          '--map5', ref, '--output', output, '--chr', '22', '--plink', executable,
          '--plink2', executable, '--shapeit', executable, '--bcftools', executable, '--impute5', executable)
status <- system2('Rscript', shQuote(args))
stopifnot(status == 0, all(file.exists(paste0(output, c('.bed', '.bim', '.fam', '.gen', '.sample')))),
          !file.exists(file.path(work, 'run', 'unsafe')))
cat('23andMe imputer integration passed\n')
