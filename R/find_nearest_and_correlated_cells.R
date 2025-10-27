find_nearest_and_correlated_cells <- function(points, data, expGeneList, k = 10, threshold = 0.5)
{
  # Get all valid cell IDs from expGeneList
  valid_cells <- colnames(expGeneList)
  
  # Filter points and data to only include valid cells
  points <- points[points$cell %in% valid_cells, ]
  data <- data[data$cell %in% valid_cells, ]
  
  # Check if there are any valid points or data
  if (nrow(points) == 0 || nrow(data) == 0) {
    warning("No valid cells found after filtering. Returning empty list.")
    return(list())
  }else{
    results <- list()
    
    for (i in 1:nrow(points)) 
    {
      point <- points[i, ]
      point_x <- point$x
      point_y <- point$y
      point_cell <- point$cell
      
      data$distance <- sqrt((data$x - point_x)^2 + (data$y - point_y)^2)
      
      nearest_points <- data %>% 
        arrange(distance) %>% 
        head(k + 1)
      
      nearest_cells <- nearest_points$cell
      
      point_expression <- expGeneList[, point_cell]
      nearest_expressions <- expGeneList[, nearest_cells, drop = FALSE]  # Keep as matrix
      
      correlations <- apply(nearest_expressions, 2, function(x) cor(x, point_expression, use = "pairwise.complete.obs"))
      correlated_cells <- nearest_cells[correlations >= threshold]
      correlated_cells <- na.omit(correlated_cells)
      
      results[[point_cell]] <- correlated_cells
    }
    
    return(results)
  }
}
