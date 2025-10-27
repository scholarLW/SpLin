calculate_all_GFai_2 = function(x, Name, metadata, window_size, maxWindow = 1000)
{
  Grid = unique(sort(metadata$window))
  m = length(Grid)
  num = rep(0 , m)
  for(i in 1:m)
  {
    ix = metadata$window == Grid[i]
    iy = which(Name %in% rownames(metadata)[ix]) 
    num[i] = mean(x[iy]) * log2(maxWindow) / log2(window_size)
  }
  return(num)
}