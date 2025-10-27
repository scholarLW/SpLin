process_regions <- function(nregion, range_min = 1, range_max = 100) {
  if (length(nregion) == 0) {
    stop("Input 'nregion' is empty. Please provide a non-empty vector.")
  }
  if (any(is.na(nregion))) {
    stop("Input 'nregion' contains missing values (NA). Please remove or handle them.")
  }
  if (!is.numeric(nregion)) {
    stop("Input 'nregion' must be a numeric vector.")
  }
  nregion <- sort(nregion)
  
  merged_regions <- list()
  
  if (length(nregion) == 1) {
    merged_regions <- append(merged_regions, list(nregion))
  } else {
    current_region <- nregion[1]
    current_subset <- c(current_region)
    
    for (i in 2:length(nregion)) {
      if (!is.na(nregion[i]) && nregion[i] == nregion[i - 1] + 1) {
        current_subset <- c(current_subset, nregion[i])
      } else {
        merged_regions <- append(merged_regions, list(current_subset))
        current_subset <- c(nregion[i])
      }
    }
    merged_regions <- append(merged_regions, list(current_subset))
  }
  
  expanded_regions <- list()
  
  for (subset in merged_regions) {
    start <- subset[1] - 1
    end <- subset[length(subset)] + 1
    expanded_subset <- c(start, subset, end)
    expanded_regions <- append(expanded_regions, list(expanded_subset))
  }
  
  final_regions <- list()
  
  for (subset in expanded_regions) {
    valid_subset <- subset[subset >= range_min & subset <= range_max]
    final_regions <- append(final_regions, list(valid_subset))
  }
  
  return(final_regions)
}
