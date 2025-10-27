kriging_interpolation <- function(data, grid_density) 
{
  df <- data.frame(
    x = data$x,
    y = data$y,
    value = data$Score
  )
  
  coordinates(df) <- ~x + y
  semivariogram <- variogram(value ~ 1, df)
  vgm_model <- fit.variogram(semivariogram, vgm(model = "Sph"))
  print(vgm_model)
  x_min <- min(df@coords[, 1])
  x_max <- max(df@coords[, 1])
  y_min <- min(df@coords[, 2])
  y_max <- max(df@coords[, 2])
  
  #grid_density×grid_density
  x_step <- (x_max - x_min) / grid_density
  y_step <- (y_max - y_min) / grid_density
  
  grid <- expand.grid(x = seq(from = x_min, to = x_max, by = x_step),
                      y = seq(from = y_min, to = y_max, by = y_step))
  
  coordinates(grid) <- ~x + y
  gridded(grid) <- TRUE
  kriging_res <- krige(value ~ 1, df, newdata = grid, model = vgm_model)
  return(kriging_res)
}
