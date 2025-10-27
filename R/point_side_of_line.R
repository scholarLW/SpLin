point_side_of_line <- function(A, C, O) 
{
  vector_AC <- C - A
  vector_AO <- O - A
  
  cross_product <- vector_AC[1] * vector_AO[2] - vector_AC[2] * vector_AO[1]
  
  if (cross_product > 0) {
    return("left")
  } else if (cross_product < 0) {
    return("right")
  } else {
    return("on_line")
  }
}
