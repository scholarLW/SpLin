polygon_centroid <- function(vertices)
{
  if (!is.matrix(vertices) || nrow(vertices) < 3 || ncol(vertices) != 2) {
    stop("The input must be a matrix containing coordinates of at least 3 vertices, with each vertex having two coordinate values, x and y.")
  }
  
  sum_x <- sum(vertices[, 1])
  sum_y <- sum(vertices[, 2])
  num_vertices <- nrow(vertices)
  centroid_x <- sum_x / num_vertices
  centroid_y <- sum_y / num_vertices
  return(c(centroid_x, centroid_y))
}
