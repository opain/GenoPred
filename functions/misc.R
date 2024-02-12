#!/usr/bin/Rscript

# Make a function to source all files in a directory
source_all <- function(directory) {
  # List all .R files in the specified directory
  r_files <- list.files(directory, pattern = "\\.R$", full.names = TRUE)

  # Source each file
  for (file in r_files) {
    source(file)
  }
}

# Make a version of cat with sep=''
cat0 <- function(...) {
  cat(..., sep = '')
}

# Detect wether object refers to a data.frame or a file
obj_or_file<-function(x, header = F, return_file = T){
  if(is.null(x)){
    output <- NULL
  } else {
    if(is.vector(x)){
      if(length(x) == 1 && file.exists(x)){
        if(return_file){
          output <- x
        } else {
          output <- fread(x, header = header)
          if(ncol(output) == 1){
            output <- output$V1
          }
        }
      } else {
        if(return_file){
          output<-tempfile()
          write.table(x, output, col.names = header, row.names=F, quote=F)
        } else {
          output <- x
        }
      }
    } else if(is.data.frame(x)){
      if(return_file){
        output<-tempfile()
        write.table(x, output, col.names = header, row.names=F, quote=F)
      } else {
        output <- x
      }
    } else {
      stop("Input must be a vector or data frame object, or a valid file path.")
    }
  }
  return(output)
}

# Record start time
test_start<-function(log_file){
  test_start.time <- Sys.time()
  log_add(log_file = log_file, message = paste0('Test started at ',as.character(test_start.time)))
  return(test_start.time)
}

# Record end time
test_finish<-function(log_file, test_start.time){
  end.time <- Sys.time()
  time.taken <- end.time - test_start.time
  sink(file = log_file, append = T)
  cat('Test run finished at',as.character(end.time),'\n')
  cat('Test duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
  sink()
}

# Create log file with standard header including script name, command line options and git repo version and commit
log_header <- function(log_file, opt, script, start.time) {
  options(width = 1000)
  sink(file = log_file, append = FALSE)

  # Fetch git repository name and latest tag
  repo_path <- system("git rev-parse --show-toplevel", intern = TRUE)
  repo_name <- basename(repo_path)
  git_tag <- system("git describe --tags", intern = TRUE)

  sink(file = log_file, append = FALSE)
  cat0(
    '#################################################################\n',
    '# ', script, '\n',
    '# For questions contact Oliver Pain (oliver.pain@kcl.ac.uk)\n',
    '#################################################################\n',
    '# Repository: ', repo_name, '\n',
    '# Version (tag): ', git_tag, '\n'
  )
  cat0('---------------\n')
  print.data.frame(opt_to_df(opt), row.names = FALSE, quote = FALSE, right = FALSE)
  cat0('---------------\n')
  cat0('Analysis started at ', as.character(start.time), '\n')
  sink()
}
