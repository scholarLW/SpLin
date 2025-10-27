calculate_intersection_point = function(C, A, D)
{
  if(C[1] == A[1])
  {
    x_O = D[1]
    k = (C[2] - D[2]) / (C[1] - D[1])
    b = A[2] - k * A[1]
    y_O = k * x_O + b
  }else
    if(C[1] == D[1])
    {
      x_O = A[1]
      k = (C[2] - A[2]) / (C[1] - A[1])
      b = D[2] - k * D[1]
      y_O = k * x_O + b
    }else{
      k = (C[2] - A[2]) / (C[1] - A[1])
      b = D[2] - k * D[1]
      k1 = (C[2] - D[2]) / (C[1] - D[1])
      b1 = A[2] - k1 * A[1]
      x_O = (b1-b)/(k-k1)
      y_O = k * x_O + b
    }
  return(c(x_O, y_O)) 
}
