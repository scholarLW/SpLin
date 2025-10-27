preAutoRegisterHE = function(markingPointsFile, cellPointsFile, idjson)
{
  suppressMessages(library(jsonlite))
  idjsonData = fromJSON(idjson)
  ow = idjsonData$original_dimensions$width
  oh = idjsonData$original_dimensions$height
  script <- system.file("python/preAutoRegisterHE.py", package = "SpLin")
  cmd <- paste0("python ", script, " -mp ", markingPointsFile, " -cp ", cellPointsFile, " -s ", ow, " ", oh)
  
  system(cmd)
  
}

