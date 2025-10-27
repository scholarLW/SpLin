# Sobel
calculate_gradient <- function(raster_data)
{
  #sobel_x <- matrix(c(-1, 0, 1, -2, 0, 2, -1, 0, 1), nrow = 3)
  sobel_x <- matrix(c(-1, 0, 1, -2, 0, 2, -1, 0, 1)/4, nrow = 3)
  sobel_y <- t(sobel_x)
  grad_x <- focal(raster_data, w = sobel_x, pad = TRUE)
  grad_y <- focal(raster_data, w = sobel_y, pad = TRUE)
  
  return(list(grad_x = grad_x, grad_y = grad_y))
}