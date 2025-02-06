INT <- function(data, pheno, k = 0.375) {
  # Check for missing values in the phenotype column
  missing_index <- which(is.na(data[[pheno]]))
  
  if (length(missing_index) > 0) {
    # Split the data into subsets: complete (non-missing) and missing
    data.complete <- data[-missing_index, ]
    data.missing <- data[missing_index, ]
  } else {
    # If no missing values, treat the whole dataset as complete
    data.complete <- data
    data.missing <- NULL
  }
  
  # Compute the rank and apply the inverse normal transformation
  n <- nrow(data.complete) # Number of complete cases
  r <- rank(data.complete[[pheno]]) # Ranks of the phenotype
  data.complete[[paste0(pheno, "_int")]] <- qnorm((r - k) / (n - 2 * k + 1))
  
  # If there are missing values, add a column with NA for the missing subset
  if (!is.null(data.missing)) {
    data.missing[[paste0(pheno, "_int")]] <- NA
    # Combine the complete and missing subsets
    return(rbind(data.complete, data.missing))
  }
  
  # If no missing values, return the complete dataset
  return(data.complete)
}