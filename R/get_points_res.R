get_points_res = function(A, B)
{
  index = NULL
  if(is.null(nrow(A)))
  {
    A = t(A)
  }
  if(is.null(nrow(B)))
  {
    B = t(B)
  }  
  n = nrow(A)
  for(i in 1:n)
  {
    point_A_i = t(A[i, ])
    diss = apply(B, 1, find_points_dis, point_A_i)
    if(sum(which(diss == 0)) == 0)
    {
      index = c(index, i)
    }
  }
  return(A[index, ])
}
