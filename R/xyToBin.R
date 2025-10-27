xyToBin = function(rds, json, nlim = 0, dropPoints = TRUE, triangleMerge = TRUE, triangle.probs = 0.8)
{
  suppressMessages(library(jsonlite))
  suppressMessages(library(mgcv))
  jsonData = jsonFilter(rds, json, nlim)
  i = which(jsonData$shapes$label == 'Wai')
  j = which(jsonData$shapes$label == 'Nei')
  waiP = jsonData$shapes$points[[i]]
  neiP = jsonData$shapes$points[[j]]	
  n1 = nrow(waiP) 
  n2 = nrow(neiP) 
  m = min(n1, n2)
  
  # bini
  binI = list()
  ID = 1
  # CA0
  C = waiP[1, ]
  D = neiP[1, ]
  for(i in 2:m)
  {
    # A
    A = waiP[i, ]
    O = calculate_intersection_point(C, A, D)
    D1 = symmetric_point(D, O)
    houDots = find_points_between(neiP, D, D1)
    if(is.null(nrow(houDots)))
    {
      houDots = t(houDots)
    }  
    if(nrow(houDots)==0)
    {
      B = D
    }else{  
      indexALLnum = apply(houDots, 1, get_points_index, neiP)
      Dindex = get_points_index(D, neiP)
      
      if(length(indexALLnum)==1)
      {
        B = D
      }else{
        indexALLnum = indexALLnum[indexALLnum>=Dindex]
        indexALLnum_ix = which(diff(indexALLnum) != 1)
        if(length(indexALLnum_ix) == 0)
        {
          indexALLnum_ix = length(indexALLnum)
        }else{
          indexALLnum_ix = indexALLnum_ix[1]
        }
        houDots = neiP[indexALLnum[1:indexALLnum_ix ], ]
        B = find_points_mindis(houDots, A)
        Bindex = get_points_index(B, neiP)
        if(Bindex < Dindex)
        {
          B = D
        }  
      }
    }
    setwai = rbind(C, A)
    setnei = get_points_between(neiP, D, B)
    OC = polygon_centroid(rbind(setwai, setnei))
    binItmp = as.data.frame(cbind(ID,t(C),t(A),t(B),t(D), t(OC)))
    colnames(binItmp) = c('ID','Cx','Cy','Ax','Ay','Bx','By','Dx','Dy','Ox','Oy')
    C = A
    D = B
    binI$nei[[i - 1]] = setnei
    binI$wai[[i - 1]] = setwai
    binI$bin[[i - 1]] = binItmp 
    ID = ID + 1  
  }
  
  waiPuse = NULL
  for(i in 1:length(binI$wai))
  {
    waiPuse = unique(rbind(waiPuse, binI$wai[[i]]))
  }
  neiPuse = NULL
  for(i in 1:length(binI$nei))
  {
    neiPuse = unique(rbind(neiPuse, binI$nei[[i]]))
  }
  
  # binF
  binF = NULL
  B = t(neiPuse[nrow(neiPuse),])
  A = t(waiPuse[nrow(waiPuse),])
  neiPres = get_points_res(neiP, neiPuse)
  waiPres = get_points_res(waiP, waiPuse)
  if(is.vector(neiPres))
  {
    neiPres = t(neiPres)
  }
  if(is.vector(waiPres))
  {
    waiPres = t(waiPres)
  }	
  if(nrow(neiPres)>0 || nrow(waiPres)>0)
  {
    if(nrow(waiPres)>0)
    {
      if(nrow(neiPres)>0)
      {
        binF = rbind(A, waiPres, neiPres[nrow(neiPres):1,], B) 
      }else{
        binF = rbind(A, waiPres, B) 
      }
    }else{
      binF = rbind(A, neiPres[nrow(neiPres):1,], B)
    }
    colnames(binF) = c('x', 'y')
  }
  
  if(triangleMerge)
  {
    nbin = length(binI$bin)
    rowS = NULL
    for(i in 1:length(binI$nei))
    {
      add = binI$nei[[i]]
      if(is.null(nrow(add)))
      {
        rowS = c(rowS, 0)
      }else{
        rowS = c(rowS, find_points_dis(add[1,], add[nrow(add),]))
      }
    }
    medlen = quantile(rowS, probs = 1-triangle.probs, na.rm = TRUE)[1]
    Index0 = which(rowS <= medlen)
    
    duiNum = list()
    j = 1
    while(length(Index0)) 
    {
      index = Index0
      index1 = NULL
      for(i in 1:length(index))
      {
        if(! index[i] %in% index1)
        {
          index1 = c(index1, index[i])
          Jmerge = index[i]
          if((index[i]+1) %in% index)
          {
            for(ij in index[i]+1:nbin)
            {
              if(! ij %in% index)
              {
                break
              }
              if(sum(rowS[Jmerge]) >= medlen)
              {
                break
              }          
              Jmerge = c(Jmerge, ij)
              index1 = c(index1, Jmerge)
            }
          }else{
            if(index[i] == 1)
            {
              Jmerge = c(Jmerge, index[i] + 1)
            }else{
              if(index[i] < nbin)
              {
                if((index[i] + 1) %in% sort(unlist(duiNum)))
                {
                  suu = index[i] - 1
                }else{
                  if((index[i] - 1) %in% sort(unlist(duiNum)))
                  {
                    suu = index[i] + 1
                  }else{
                    if(rowS[index[i] + 1] < rowS[index[i] - 1])
                    {
                      suu = index[i] + 1
                    }else{
                      suu = index[i] - 1
                    }
                  }
                }
                Jmerge = c(Jmerge, suu)
              }
            }
          }
          if(length(Jmerge)>1)
          {
            duiNum[[j]] = sort(Jmerge)
            j = j + 1 
          }
        }
      }
      Index0 = setdiff(index, index1)
    }
    Index = unique(sort(unlist(duiNum)))
    
    binInew = list()
    keepindex = setdiff(1:length(binI$bin), Index)
    for(i in keepindex)
    {
      binInew$nei[[i]] = binI$nei[[i]]
      binInew$wai[[i]] = binI$wai[[i]]
      binInew$bin[[i]] = binI$bin[[i]]
    }
    nm = length(duiNum)
    
    for(mi in 1:nm)
    {
      tmddd = duiNum[[mi]]
      neimitmp = NULL
      waimitmp = NULL
      for(pij in 1:length(tmddd))
      {
        atde = binI$nei[[tmddd[pij]]]
        if(is.null(nrow(atde)))
        {
          atde = t(atde)
        }
        neimitmp = unique(rbind(neimitmp, atde))
        waimitmp = unique(rbind(waimitmp, binI$wai[[tmddd[pij]]]))
      }
      binmitmp = binI$bin[[tmddd[1]]]
      OC = polygon_centroid(rbind(neimitmp, waimitmp))
      binmitmp$Ax = binI$bin[[tmddd[length(tmddd)]]]$Ax
      binmitmp$Ay = binI$bin[[tmddd[length(tmddd)]]]$Ay
      binmitmp$Bx = binI$bin[[tmddd[length(tmddd)]]]$Bx
      binmitmp$By = binI$bin[[tmddd[length(tmddd)]]]$By
      binmitmp$Ox = OC[1]
      binmitmp$Oy = OC[2]
      binInew$nei[[tmddd[1]]] = neimitmp
      binInew$wai[[tmddd[1]]] = waimitmp
      binInew$bin[[tmddd[1]]] = binmitmp
    }
    binI = binInew   
  }
  
  binInew = list()
  j = 1
  for(i in 1:length(binI$bin))
  {
    tmpik = binI$bin[[i]]
    if(!is.null(tmpik))
    {
      binInew$nei[[j]] = binI$nei[[i]]
      binInew$wai[[j]] = binI$wai[[i]]
      binInew$bin[[j]] = binI$bin[[i]]
      j = j + 1
    }
  }
  C = binInew$bin[[1]][,c('Cx','Cy')]
  A = binInew$bin[[1]][,c('Ax','Ay')]
  D = binInew$bin[[1]][,c('Dx','Dy')]
  pointsside = point_side_of_line(C, A, D)
  duiNum = list()
  j = 1
  index = NULL
  for(i in 1:(length(binInew$bin)-1))
  {
    index = c(index, i)
    ij = i + 1
    O = binInew$bin[[ij]][,c('Ax','Ay')]
    pointssidetmp = point_side_of_line(C, A, O)
    A = binInew$bin[[ij]][,c('Ax','Ay')]
    if(pointssidetmp == pointsside)
    {
      duiNum[[j]] = index
      C = binInew$bin[[ij]][,c('Cx','Cy')]
      j = j + 1
      index = NULL
    }
  }
  duiNum[[j]] = c(index, length(binInew$bin))
  
  binInew1 = list()
  for(i in 1:length(duiNum))
  {
    tmddd = duiNum[[i]]
    if(length(tmddd)==1)
    {
      binInew1$nei[[i]] = binInew$nei[[tmddd]]
      binInew1$wai[[i]] = binInew$wai[[tmddd]]
      binInew1$bin[[i]] = binInew$bin[[tmddd]]
    }else{
      neimitmp = NULL
      waimitmp = NULL
      for(pij in 1:length(tmddd))
      {
        atde = binInew$nei[[tmddd[pij]]]
        if(is.null(nrow(atde)))
        {
          atde = t(atde)
        }
        neimitmp = unique(rbind(neimitmp, atde))
        waimitmp = unique(rbind(waimitmp, binInew$wai[[tmddd[pij]]]))
      }
      binmitmp = binInew$bin[[tmddd[1]]]
      OC = polygon_centroid(rbind(neimitmp, waimitmp))
      binmitmp$Ax = binInew$bin[[tmddd[length(tmddd)]]]$Ax
      binmitmp$Ay = binInew$bin[[tmddd[length(tmddd)]]]$Ay
      binmitmp$Bx = binInew$bin[[tmddd[length(tmddd)]]]$Bx
      binmitmp$By = binInew$bin[[tmddd[length(tmddd)]]]$By
      binmitmp$Ox = OC[1]
      binmitmp$Oy = OC[2]
      binInew1$nei[[i]] = neimitmp
      binInew1$wai[[i]] = waimitmp
      binInew1$bin[[i]] = binmitmp
    }
  }
  binI = binInew1 
  
  absolute_path <- dirname(normalizePath(json))
  
  if(dropPoints)
  {	  
    pdf(paste0(absolute_path, '/BinRegion.pdf'), 8, 8)
    plot(waiPuse, xlab = '', ylab = '', col = 'black', pch = 18)
    points(neiPuse,col='grey', pch = 20)
    ni = length(binI$bin)
    for(i in 1:ni)
    {
      tmpdaa = binI$bin[[i]]
      x_coords <- c(tmpdaa$Cx, tmpdaa$Ax, tmpdaa$Bx, tmpdaa$Dx)
      y_coords <- c(tmpdaa$Cy, tmpdaa$Ay, tmpdaa$By, tmpdaa$Dy)		
      polygon(x_coords, y_coords, col="grey90", border="black") 
      # Sys.sleep(1) 
    }
    if(!is.null(binF))
    {
      polygon(binF[,1], binF[, 2], col="lightblue", border="red") 
    }
    dev.off()  
  }
  
  BinData = list('binI' = binI, 'binF' = binF)
  save(BinData, file = paste0(absolute_path,'/xyToBin.RData'))
  
  return(BinData)
}
