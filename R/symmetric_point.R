symmetric_point = function(D, O)
{
  x <- D[1]
  y <- D[2]
  x1 <- O[1]
  y1 <- O[2]
  
  x_s <- 2 * x1 - x
  y_s <- 2 * y1 - y
  
  return(c(x_s, y_s))
}
