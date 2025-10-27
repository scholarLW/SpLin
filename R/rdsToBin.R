rdsToBin = function(rds, json, nlim = 0, dropPoints = TRUE, triangleMerge = TRUE, triangle.probs = 0.8)
{
  suppressMessages(library(Seurat))
  suppressMessages(library(mgcv))
  BinData = xyToBin(rds, json, nlim = nlim, dropPoints = dropPoints, triangleMerge = triangleMerge, triangle.probs = triangle.probs)
  
  obj = readRDS(rds)	
  metaDta = obj@meta.data
  absolute_path <- dirname(normalizePath(json))
  
  metaDta$bin = NA
  n = length(BinData$binI$bin)
  tmpBinI = NULL
  for(i in 1:n)
  {
    if(!is.null(BinData$binI$bin[[i]]))
    {
      atmp = t(BinData$binI$bin[[i]])
      tmpBinI = as.data.frame(rbind(tmpBinI, t(atmp)))
      ID = atmp[1]
      neitmp = BinData$binI$nei[[i]]
      if(is.null(nrow(neitmp)))
      {
        neitmp = t(as.matrix(neitmp))
      }
      waitmp = BinData$binI$wai[[i]]
      bnd = rbind(waitmp, neitmp[nrow(neitmp):1,])
      inside <- in.out(bnd, as.matrix(metaDta[, c('x', 'y')])) 
      metaDta$bin[inside] = ID 
    }
  }
  
  if(!is.null(BinData$binF))
  {
    tmp = BinData$binF
    inside <- in.out(tmp, as.matrix(metaDta[, c('x', 'y')])) 
    metaDta$bin[inside] = 'binF'
  }
  
  tmp = tmpBinI
  nbin = nrow(tmp)
  binlenALL = rep(0, nbin)
  sumCA = 0
  
  for(i in 1:nbin)
  {
    binlenALL[i] = sqrt((tmp$Cx[i]-tmp$Ax[i])^2+(tmp$Cy[i]-tmp$Ay[i])^2)
    tmp$Cx1[i] = sumCA
    tmp$Cy1[i] = 0
    sumCA = sumCA + binlenALL[i]
    tmp$Ax1[i] = sumCA
    tmp$Ay1[i] = 0
    tmp$Bx1[i] = NA
    tmp$By1[i] = NA
    tmp$Dx1[i] = NA
    tmp$Dy1[i] = NA
  }
  
  metaDta = cbind(rownames(metaDta), metaDta)
  colnames(metaDta)[1] = 'Cell'
  metaDta$xnew = NA
  metaDta$ynew = NA
  pdf(paste0(absolute_path, '/BinRegionData.pdf'), 8, 8)
  plot(metaDta$x, metaDta$y, type = "n", xlab = "", ylab = "")
  binName = sort(unique(metaDta$bin))
  ni = length(binName)
  for(i in 1:ni)
  {
    ix = metaDta$bin == binName[i]
    tmpdaa = as.data.frame(metaDta[ix,])
    tmpdaa = tmpdaa[!is.na(tmpdaa$Cell), ]
    
    x_coords <- tmpdaa$x
    y_coords <- tmpdaa$y	
    chull_indices <- chull(x_coords, y_coords)
    chull_x <- x_coords[chull_indices]
    chull_y <- y_coords[chull_indices]
    
    if(binName[i] == 'binF')
    {
      polygon(chull_x, chull_y, col="lightblue", border="red") 
    }else{
      polygon(chull_x, chull_y, col = "grey90", border = "black")
    }
    
  }
  dev.off()
  res = list('binI'=tmp, 'binF'=BinData$binF, 'metaDta'=metaDta)	
  return(res)
}
