preAutoRegister = function(json, cellPointsFile, idjson, Signal = FALSE)
{
  suppressMessages(library(jsonlite))
  idjsonData = fromJSON(idjson)
  ow = idjsonData$original_dimensions$width
  oh = idjsonData$original_dimensions$height
  script <- system.file("python/preAutoRegister.py", package = "SpLin")
  
  cmd <- paste0("python ", script, " -j ", json, " -cp ", cellPointsFile, " -s ", ow, " ", oh)
  if (Signal) 
  {
    cmd <- paste0(cmd, " -S")
  } else {
    cmd <- cmd
  }
  system(cmd)
}

