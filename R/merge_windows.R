merge_windows <- function(data, pv = 0.000001)
{
  windows <- unique(data$window)
  regions <- list()
  current_region <- list(windows[1])  
  region_map <- list()  
  
  for (i in 1:(length(windows) - 1)) {
    current_window <- windows[i]
    next_window <- windows[i + 1]
    
    current_scores <- data %>% filter(window == current_window) %>% pull(Score)
    next_scores <- data %>% filter(window == next_window) %>% pull(Score)
    
    if (length(current_scores) > 1 && length(next_scores) > 1) {
      test_result <- wilcox.test(current_scores, next_scores)
      p_value <- test_result$p.value
    } else {
      p_value <- 0
    }
    
    if (p_value > pv) {
      current_region <- c(current_region, next_window)
    } else {
      regions <- append(regions, list(current_region))
      current_region <- list(next_window)
    }
  }
  
  if (length(current_region) > 0) {
    regions <- append(regions, list(current_region))
  }
  
  for (i in seq_along(regions)) {
    for (window in regions[[i]]) {
      region_map[[as.character(window)]] <- i
    }
  }
  
  return(region_map)
}
