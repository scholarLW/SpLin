jsonFilter = function(rds, json, nlim)
{
  suppressMessages(library(Seurat))
  suppressMessages(library(jsonlite))
  suppressMessages(library(mgcv))
  obj = readRDS(rds)
  jsonData = fromJSON(json)
  metaDta = obj@meta.data[, c('x', 'y')]
  i = which(jsonData$shapes$label == 'Wai')
  j = which(jsonData$shapes$label == 'Nei')	
  n = c(i, j)
  newlist = list()
  for(ij in 1:2)
  {
    i = n[ij] 
    tmp = as.matrix(jsonData$shapes$points[[i]])
    if(ij == 1)
    {
      max_length = nrow(tmp)
      j = 1
      while(j <= (max_length - 2))
      {
        if((j+2)<=max_length)
        {
          bnd = tmp[j:(j+2), ]
          inside <- in.out(bnd, as.matrix(metaDta))
          if(sum(inside)<=nlim)
          {
            tmp = tmp[-(j+1), ]
            max_length = nrow(tmp)
          }else{
            j = j + 1
          } 
        }
      } 
    }
    newlist[[i]] = tmp
  }
  jsonData$shapes$points = newlist
  
  return(jsonData)
}
