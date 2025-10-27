find_convex_regions <- function(data) 
{
  data <- data[order(data$region), ]
  convex_regions <- c()
  
  if(nrow(data)>2)
  {
    for (i in 2:(nrow(data) - 1)) 
    {
      current_mean <- data$mean[i]
      left_mean <- data$mean[i - 1]
      right_mean <- data$mean[i + 1]
     
      if (current_mean > left_mean && current_mean > right_mean) 
      {
        convex_regions <- c(convex_regions, data$region[i])
      }
    }
  }
  
  return(convex_regions)
}
