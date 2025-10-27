RapidRegisterOUT = function(rds, json, idjson, pixelssDNAFile = NULL, pixelHEFile = NULL, index = 1)
{
  labelfont = 1
  labelfontsize = 2
  axistitlefontsize = 3
  spatialpointsize = 0.1
  pointalpha = 1
  pointshape = 16
  legendtextsize = 10
  legendtitlesize = 25
  guidelegendsize = 6
  
  suppressMessages(library(Seurat))
  suppressMessages(library(ggplot2))
  suppressMessages(library(RColorBrewer))
  suppressMessages(library(jsonlite))
  suppressMessages(library(mgcv))
  suppressMessages(library(data.table))
  
  idjsonData = fromJSON(idjson)
  ow = idjsonData$original_dimensions$width
  oh = idjsonData$original_dimensions$height
  absolute_path <- dirname(normalizePath(json))
  JsonFile = paste0(absolute_path, '/Adjusted_Output_Manual.json')
  script <- system.file("python/reRapidRegister.py", package = "SpLin") 
  
  cmd <- paste0("python ", script, " -j ", json, " -s ", ow, " ", oh, " -i ", index)
  
  if (!is.null(pixelssDNAFile)) 
  {
    cmd <- paste0(cmd, " -ssDNA ", shQuote(pixelssDNAFile))
  }
  if (!is.null(pixelHEFile)) 
  {
    cmd <- paste0(cmd, " -HE ", shQuote(pixelHEFile))
  }	
  system(cmd)
  
  
  obj = readRDS(rds)
  jsonData = fromJSON(JsonFile)
  i = which(jsonData$shapes$label == 'Wai')
  j = which(jsonData$shapes$label == 'Nei')
  n = nrow(jsonData$shapes$points[[j]])
  bnd = rbind(jsonData$shapes$points[[i]], jsonData$shapes$points[[j]][n:1,])	
  
  inside <- in.out(bnd, as.matrix(obj@meta.data[, c('x', 'y')]))
  cellNames = rownames(obj@meta.data)[inside]
  tmpobj = subset(obj, cells = cellNames)
  outrds = paste0(absolute_path, '/Update_Manual_',basename(normalizePath(rds)))
  saveRDS(tmpobj, file = outrds)
  if(!is.null(pixelssDNAFile))
  {
    absolute_path <- dirname(normalizePath(pixelssDNAFile))	
    ssDNA = paste0(absolute_path, '/pixel_markingPoints_XY.txt')
    ssDNAda = as.data.frame(fread(ssDNA))	
    png(paste0(absolute_path, '/pixel_markingPointsMatch.png'), width = 500, height = 500)
    plot(ssDNAda[,c(1,2)], col = rgb(1, 0, 0, 0.3), xlim = range(c(ssDNAda[,1], bnd[,1])), ylim = range(c(bnd[,2], ssDNAda[,2])), xaxt = "n", yaxt = "n", bty = "n", xlab = '', ylab = '')
    polygon(bnd[,c(1,2)], col = rgb(0.85, 0.85, 0.85, 0.3))
    dev.off()	
  }
  if(!is.null(pixelHEFile))
  {
    absolute_path <- dirname(normalizePath(pixelHEFile))		
    HE = paste0(absolute_path, '/pixel_markingPoints_XY.txt')
    HEda = as.data.frame(fread(HE))	
    png(paste0(absolute_path, '/pixel_markingPointsMatch.png'), width = 500, height = 500)
    plot(HEda[,c(1,2)], col = rgb(1, 0, 0, 0.3), xlim = range(c(HEda[,1], bnd[,1])), ylim = range(c(bnd[,2], HEda[,2])), xaxt = "n", yaxt = "n", bty = "n", xlab = '', ylab = '')
    polygon(bnd[,c(1,2)], col = rgb(0.85, 0.85, 0.85, 0.3))
    dev.off()	
  }	
}

