reSizeIMG = function(ori, resi = NULL, dpi = 300, it = 'Signal')
{
  script <- system.file("python/reSIMG.py", package = "SpLin")
  
  cmd <- paste0("python ", script, " -I ", ori, " -d ", dpi, " -t ", it)
  
  if (!is.null(resi)) 
  {
    cmd <- paste0(cmd, " -O ", shQuote(resi))
  } else {
    cmd <- cmd
  }
  system(cmd)  
  
}
















