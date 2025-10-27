coordLinearEXP = function(rds, json, nlim = 0, drop = TRUE, triangleMerge = TRUE, triangle.probs = 0.2)
{
  suppressMessages(library(Seurat))
  suppressMessages(library(data.table))
  res = rdsToBin(rds, json, nlim = nlim, dropPoints = drop, triangleMerge = triangleMerge, triangle.probs = 1 - triangle.probs)
  tmp = res$binI
  binF = res$binF
  metaDta = res$metaDta
  
  absolute_path <- dirname(normalizePath(json))
  metaDtaPointsFile = paste0(absolute_path, '/metaDtaPoints.txt')
  outputFile = paste0(absolute_path, '/metaDtaPointsUpdate.txt')
  BinFPointsFile = paste0(absolute_path, '/BinFPoints.txt')
  BinIPointsFile = paste0(absolute_path, '/BinIPoints.txt')
  
  write.table(tmp, file = BinIPointsFile, quote = F, col.names = T, row.names = F, sep= '\t')
  if(file.exists(BinFPointsFile))
  {
    file.remove(BinFPointsFile)
  }
  if(!is.null(binF))
  {
    write.table(binF, file = BinFPointsFile, quote = F, col.names = T, row.names = F, sep= '\t')
  }else{
    BinFPointsFile = NULL
  }
  write.table(metaDta[, c('Cell','x','y','bin','xnew','ynew')], file = metaDtaPointsFile, quote = F, col.names = T, row.names = F, sep= '\t')
  
  script <- system.file("python/coordLinear.py", package = "SpLin")
  cmd <- paste0("python ", script, " -mdp ", metaDtaPointsFile, " -o ", outputFile, " -bip ", BinIPointsFile)
  
  if (!is.null(BinFPointsFile)) 
  {
    cmd <- paste0(cmd, " -bfp ", shQuote(BinFPointsFile))
  } else {
    cmd <- cmd
  }
  system(cmd)	
  
  obj = readRDS(rds)
  BinData = as.data.frame(fread(outputFile))
  rownames(BinData) = BinData$Cell
  BinData = BinData[, -c(1,2,3)]
  merged_data <- merge(obj@meta.data, BinData, by = "row.names", all = TRUE, sort = FALSE)
  rownames(merged_data) <- merged_data$Row.names
  merged_data = merged_data[, -1]
  obj@meta.data = merged_data
  absolute_path <- dirname(normalizePath(rds))
  saveRDS(obj, file = paste0(absolute_path, '/SpLin_EXP_output.rds'))
}
