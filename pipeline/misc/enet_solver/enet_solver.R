# Load required package
library(glmnet)

# Function to compute elastic net weights using correlation information with soft thresholding
compute_elastic_net_weights_with_thresholding <- function(correlations, predictor_corr_matrix, alpha, lambda, max_iter = 100, tol = 1e-6) {
  # Check that dimensions match
  n <- length(correlations)
  if (ncol(predictor_corr_matrix) != n || nrow(predictor_corr_matrix) != n) {
    stop("Dimensions of predictor correlation matrix must match the length of the correlation vector.")
  }

  # Define the Ridge penalty matrix
  ridge_penalty <- (1 - alpha) * lambda * diag(n)

  # Initialize weights
  beta <- rep(0, n)

  # Iterate for coordinate descent-like updates with soft-thresholding
  for (iter in 1:max_iter) {
    beta_old <- beta

    # Update each coefficient using soft-thresholding
    for (j in 1:n) {
      # Calculate the adjusted correlation with residuals for variable j
      residual_correlation <- correlations[j] - sum(predictor_corr_matrix[j, -j] * beta[-j])

      # Apply the soft-thresholding step
      beta[j] <- soft_threshold(residual_correlation / (1 + ridge_penalty[j, j]), alpha * lambda)
    }

    # Check for convergence
    if (max(abs(beta - beta_old)) < tol) {
      break
    }
  }

  return(as.vector(beta))
}

# Soft thresholding function for Lasso-like behavior
soft_threshold <- function(x, threshold) {
  sign(x) * pmax(abs(x) - threshold, 0)
}

# Example usage:
# Simulate data for comparison
set.seed(123)
n_samples <- 10000
n_predictors <- 100

# Generate random predictors
X <- matrix(rnorm(n_samples * n_predictors), ncol = n_predictors)
# Generate a random outcome with some noise
y <- X %*% rnorm(n_predictors) + rnorm(n_samples)

# Compute correlations between predictors and outcome
correlations <- as.vector(cor(X, y))

# Compute the correlation matrix of predictors
predictor_corr_matrix <- cor(X)

# Define a range of alphas and lambdas to test
alphas <- c(0.5)  # Use a fixed alpha value for simplicity in this comparison
lambdas <- c(0.4)  # Different regularization strengths

# Compute custom elastic net weights
custom_weights <- compute_elastic_net_weights_with_thresholding(correlations, predictor_corr_matrix, alphas, lambdas)

# Compare with glmnet
# Use a single alpha and a range of lambda values
glmnet_fit <- glmnet(X, y, alpha = alphas, lambda = lambdas, standardize = TRUE)

# Extract glmnet coefficients for the selected lambda values
glmnet_weights <- as.matrix(coef(glmnet_fit))[-1, ]  # Drop the intercept

plot(custom_weights, glmnet_weights)
cor(cbind(custom_weights, glmnet_weights))

# Convert custom weights data frame to match glmnet output structure for comparison
custom_weights_glmnet_format <- custom_weights[, paste0("alpha_0.5_lambda_", lambdas)]

# Combine into a single data frame for comparison
comparison <- data.frame(
  Predictor = paste0("X", 1:n_predictors),
  Custom_Lambda_0.1 = custom_weights_glmnet_format[, 1],
  Custom_Lambda_0.2 = custom_weights_glmnet_format[, 2],
  Custom_Lambda_0.5 = custom_weights_glmnet_format[, 3],
  GLMNet_Lambda_0.1 = glmnet_weights[, 1],
  GLMNet_Lambda_0.2 = glmnet_weights[, 2],
  GLMNet_Lambda_0.5 = glmnet_weights[, 3]
)

# Display the comparison
print(comparison)
