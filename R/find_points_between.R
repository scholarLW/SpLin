find_points_between <- function(P, D, O)
{   
  if (!is.matrix(P) || ncol(P) != 2) {  
    stop("P must be a matrix with two columns")  
  }  
  
  D <- as.numeric(D)  
  O <- as.numeric(O)  
  if (length(D) != 2 || length(O) != 2) {  
    stop("D and O must be two-dimensional numeric vectors")  
  }  
  
  D_mat <- matrix(D, nrow = 1, byrow = TRUE)  
  O_mat <- matrix(O, nrow = 1, byrow = TRUE)  
  DO <- O_mat - D_mat  
  P_vec <- P - rep(D_mat, each = nrow(P), byrow = TRUE)  
  DO_norm <- sqrt(sum(DO^2))    
  proj_lengths <- P_vec %*% t(DO) / DO_norm   
  between_indices <- which(proj_lengths >= 0 & proj_lengths <= DO_norm)  
  
  return(P[between_indices, ])  
}  
