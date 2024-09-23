#!/usr/bin/Rscript

# Create the function for converting R2 into liability R2
h2l_R2 <- function(k, r2, p) {
  # K baseline disease risk
  # r2 from a linear regression model attributable to genomic profile risk score
  # P proportion of sample that are cases
  # calculates proportion of variance explained on the liability scale
  #from ABC at http://www.complextraitgenomics.com/software/
  #Lee SH, Goddard ME, Wray NR, Visscher PM. (2012) A better coefficient of determination for genetic profile analysis. Genet Epidemiol. 2012 Apr;36(3):214-24.
  x= qnorm(1-k)
  z= dnorm(x)
  i=z/k
  C= k*(1-k)*k*(1-k)/(z^2*p*(1-p))
  theta= i*((p-k)/(1-k))*(i*((p-k)/(1-k))-x)
  h2l_R2 = C*r2 / (1 + C*theta*r2)
}

# Function for merging by IID
predictor_merger <- function(x, y) {
  return(merge(x, y, by = 'IID'))
}

# Function to combine FID and IID, into a variable called IID
combine_fid_iid <- function(data, fid_col = "FID", iid_col = "IID") {
  # Check if the specified columns exist in the data
  if (!all(c(fid_col, iid_col) %in% names(data))) {
    stop("Both FID and IID columns must be present in the data.")
  }

  # Combine FID and IID into a new IID column
  data[[iid_col]] <- paste0(data[[fid_col]], ':', data[[iid_col]])

  # Remove the FID column
  data[[fid_col]] <- NULL

  return(data)
}

# Remove variables with > opt$pred_miss missing values
filter_columns_by_missing <- function(data, threshold, first_col_keep = TRUE) {
  # Identify columns to keep based on missing value threshold
  col_keep <- c(first_col_keep, sapply(data[, -1, with = FALSE], function(col) {
    mean(!is.finite(col) | is.na(col)) < threshold
  }))

  # Filter the data by keeping only the columns that meet the criteria
  data_filtered <- data[, col_keep, with = FALSE]

  return(data_filtered)
}

# Set the grid search for the elastic net
# This is expanded from the default to allow more aggressive shrinkage
enet_grid <- expand.grid(
  alpha = seq(0, 1, length = 5),          # Explore alpha values: 0 (Ridge) to 1 (Lasso)
  lambda = 10^seq(-4, 1, length = 10)     # Explore lambda values: 0.0001 to 10
)
