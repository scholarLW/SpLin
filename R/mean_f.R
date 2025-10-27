mean_f = function(x, Name, metadata)
{
  m = max(metadata$window, na.rm = TRUE)
  num = rep(0 , m)
  for(i in 1:m)
  {
    ix = metadata$window == i
    iy = which(Name %in% rownames(metadata)[ix]) 
    num[i] = mean(x[iy])
  }
  return(num)
}