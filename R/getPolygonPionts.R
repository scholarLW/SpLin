getPolygonPionts = function(rds, imageFile, json, idjson, Signal = FALSE)
{
  suppressMessages(library(Seurat))
  suppressMessages(library(jsonlite))
  idjsonData = fromJSON(idjson)
  ow = idjsonData$original_dimensions$width
  oh = idjsonData$original_dimensions$height
  rw = idjsonData$scaled_dimensions$width
  rh = idjsonData$scaled_dimensions$height	
  
  obj = readRDS(rds)
  absolute_path <- dirname(normalizePath(json))
  JsonFile = paste0(absolute_path, '/Adjusted_',basename(normalizePath(json)))
  
  metaData = obj@meta.data[, c('x', 'y')]
  pdf(paste0(absolute_path, '/cellPoints.pdf'))
  plot(metaData, axes = F, xlab = '', ylab = '', pch = 20)
  dev.off()	
  write.table(metaData, file = paste0(absolute_path, '/cellPoints.txt'), quote = F, col.names = F, row.names = F, sep= '\t')	
  
  if(is.null(oh) || is.null(ow))
  {
    if(Signal)
    {
      oh = ceiling(rh*1.2)
      ow = ceiling(rw*1.2)
    }else{
      tmpn = max(max(metaData$x), max(metaData$y))
      if(rw > rh)
      {
        oh = ceiling(1.2 * tmpn)
        ow = ceiling(1.5 * tmpn)
      }else{
        oh = ceiling(1.5 * tmpn)
        ow = ceiling(1.2 * tmpn)
      }
    }
    idjsonDataN <- list(
      scaled_dimensions = list(
        width = as.numeric(rw),
        height = as.numeric(rh)
      ),
      original_dimensions = list(
        width = as.numeric(ow),
        height = as.numeric(oh)
      )
    )
    write_json(idjsonDataN, idjson, pretty = TRUE)
  }
  
  script <- system.file("python/coordMap.py", package = "SpLin")
  system(paste0("python ", script, " -j ", json, " -o ", JsonFile, " -nw ", rw, " -nh ", rh, " -ow ", ow, " -oh ", oh))
  
  jsonData = fromJSON(JsonFile)	
  if(length(jsonData$shapes$points)>1)
  {
    if(isTRUE(Signal))
    {
      i = which(jsonData$shapes$label == 'Sign')
      bnd = jsonData$shapes$points[[i]]
    }else{
      i = which(jsonData$shapes$label == 'Wai')
      j = which(jsonData$shapes$label == 'Nei')
      n = nrow(jsonData$shapes$points[[j]])
      bnd = rbind(jsonData$shapes$points[[i]], jsonData$shapes$points[[j]][n:1,])			
    }
    pdf(paste0(absolute_path, '/markingPoints.pdf'))
    plot(bnd, axes = F, xlab = '', ylab = '', pch = 20)
    dev.off()
    write.table(bnd, file = paste0(absolute_path, '/markingPoints.txt'), quote = F, col.names = F, row.names = F, sep= '\t')
  }else{
    bnd = jsonData$shapes$points[[1]]
    pdf(paste0(absolute_path, '/markingPoints.pdf'))
    plot(bnd, axes = F, xlab = '', ylab = '', pch = 20, cex = 1)
    dev.off()
    write.table(bnd, file = paste0(absolute_path, '/markingPoints.txt'), quote = F, col.names = F, row.names = F, sep= '\t')
  }
  
  if(!isTRUE(Signal))
  {	
    script <- system.file("python/coordXY.py", package = "SpLin")
    markingPointsFile = paste0(absolute_path, '/Adjusted_output_coordinates_and_colors.txt')
    system(paste0("python ", script, " -j ", json, " -i ", imageFile, " -o ", markingPointsFile, " -ww ", rw, " -wh ", rh, " -ow ", ow, " -oh ", oh))
  }
  cat("Job Done!\n")
}