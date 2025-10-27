get_points_between <- function(P, D, O)
{
  index1 = which(P[, 1] == D[1] & P[, 2] == D[2])
  index2 = which(P[, 1] == O[1] & P[, 2] == O[2])
  return(P[index1:index2, ])  
}
