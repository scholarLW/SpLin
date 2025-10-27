find_points_within_distance <- function(points, hull_points, median_distance)
{
  results <- data.frame(cell = character(), x = numeric(), y = numeric())
  for (i in 1:nrow(points)) {
    point <- points[i, ]
    point_cell <- point$cell
    point_x <- point$x
    point_y <- point$y
    distances <- sqrt((hull_points$x - point_x)^2 + (hull_points$y - point_y)^2)

    if (any(distances <= median_distance)) 
    {
      results <- rbind(results, point)
    }
  }
  
  return(results)
}
