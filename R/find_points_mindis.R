find_points_mindis <- function(P, B)
{
  index = NULL
  if(is.null(nrow(P)))
  {
    P = t(P)
  }
  n = nrow(P)
  for(i in 1:n)
  {
    P0 = P[i, ]
    index = c(index, find_points_dis(P0, B))
  }
  return(P[which.min(index), ])  
}
