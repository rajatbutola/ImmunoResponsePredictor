################################################################################
# preprocess.R
# Script to apply "ComBat + Quantile" normalization on test data (Count matrix)
# Input:
# - train.csv (zipped)
# - test.csv (zipped)
# Output:
# - test.csv (normalized and ready for prediction)
################################################################################

# Load required libraries
suppressPackageStartupMessages({
  library(sva)               # For ComBat
  library(preprocessCore)    # For quantile normalization
  library(utils)
})

# Set paths (adjust as needed)
input_dir <- "."  # Working directory with zipped expression data
train_file <- unz(file.path(input_dir, "log2TPMp1_train.zip"), "log2TPMp1_train.csv")
test_file  <- unz(file.path(input_dir, "log2TPMp1_test.zip"),  "log2TPMp1_test.csv")

# Load expression data
expr_train <- read.csv(train_file, row.names = 1, check.names = FALSE)
expr_test  <- read.csv(test_file, row.names = 1, check.names = FALSE)

# Transpose to samples x genes (as required by ComBat)
expr_train_t <- t(expr_train)
expr_test_t  <- t(expr_test)

# Create batch vector
batch_vector <- c(rep("train", nrow(expr_train_t)), rep("test", nrow(expr_test_t)))

# Combine data
combined_expr <- rbind(expr_train_t, expr_test_t)

# Apply ComBat
cat("Applying ComBat normalization...\n")
combat_expr <- ComBat(
  dat = t(t(combined_expr)), 
  batch = batch_vector, 
  mod = NULL, 
  par.prior = TRUE, 
  prior.plots = FALSE
)

# Apply Quantile Normalization
cat("Applying Quantile normalization...\n")
qn_expr <- normalize.quantiles(combat_expr, copy = TRUE)

# Restore row and column names
colnames(qn_expr) <- colnames(combat_expr)
rownames(qn_expr) <- gsub("^X", "", rownames(combat_expr))

# Transpose back to gene x sample
qn_expr_t <- t(qn_expr)

# Extract only test samples
n_train <- nrow(expr_train_t)
qn_test <- qn_expr_t[(n_train + 1):nrow(qn_expr_t), ]

# Save normalized test set
write.csv(qn_test, file = "test.csv", row.names = TRUE)  # Now this file is ready for uploading into the GUI on server. 
cat("Preprocessed test data saved to test.csv\n")
