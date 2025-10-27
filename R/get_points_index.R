get_points_index = function(B, A)
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
    diss = find_points_dis(B, point_A_i)
    if(diss == 0)
    {
      index = i
      break
    }
  }
  return(index)
}
