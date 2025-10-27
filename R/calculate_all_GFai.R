calculate_all_GFai = function(x, Name, metadata, window_size, maxWindow = 1000)
{
  m = max(metadata$window, na.rm = TRUE)
  num = rep(0 , m)
  for(i in 1:m)
  {
    ix = metadata$window == i
    iy = which(Name %in% rownames(metadata)[ix]) 
    num[i] = mean(x[iy]) * log2(maxWindow) / log2(window_size)
  }
  return(num)
}