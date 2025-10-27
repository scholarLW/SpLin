SpotBinPixel = function(pixelPointsFile, binRD, BinFPointsFile = NULL, BinIPointsFile = NULL)
{
  suppressMessages(library(data.table))
  suppressMessages(library(jsonlite))
  suppressMessages(library(mgcv))
  
  absolute_path <- dirname(normalizePath(pixelPointsFile))
  metaDta = as.data.frame(fread(pixelPointsFile))
  metaDtaPointsFile = paste0(absolute_path, '/pixelPoints.txt')
  load(binRD)
  colnames(metaDta) = c('x', 'y', 'Cell')
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
  metaDta$xnew = NA
  metaDta$ynew = NA
  write.table(metaDta[, c('Cell','x','y','bin','xnew','ynew')], file = metaDtaPointsFile, quote = F, col.names = T, row.names = F, sep= '\t')
  
  script <- system.file("python/coordLinear.py", package = "SpLin")
  outputFile = paste0(absolute_path, '/pixelPointsUpdate.txt')
  cmd <- paste0("python ", script, " -mdp ", metaDtaPointsFile, " -o ", outputFile, " -bip ", BinIPointsFile)
  
  if (!is.null(BinFPointsFile)) 
  {
    cmd <- paste0(cmd, " -bfp ", shQuote(BinFPointsFile))
  } else {
    cmd <- cmd
  }
  system(cmd)	
  
}
