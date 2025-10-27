rdsMerge = function(lassoRDS, celltypeRDS)
{
  suppressMessages(library(Seurat))
  seuratobj <- readRDS(lassoRDS)
  celltypeobj <- readRDS(celltypeRDS)
  metac = celltypeobj@meta.data
  metal = seuratobj@meta.data
  matches = NULL
  aname = NULL
  bname = NULL  
  istep = 1
  while(!sum(matches) && istep <= 100000)
  {
    a <- sample(rownames(metal), size = 10, replace = FALSE)
    b <- sample(rownames(metac), size = 10, replace = FALSE)
    matches <- sapply(a, function(x) any(sapply(b, function(y) grepl(x, y))))
    if(sum(matches))
    {
      aname <- a[matches][1]
      bname <- b[grepl(aname, b)]
    }else{
      istep = istep + 1
    }
  }
  if(is.null(aname))
  {
    stop("The input RDS file does not match, please confirm.")
  }
  dif = sub(aname,'',bname)
  rownames(metal) = paste0(rownames(metal), dif)
  
  samnames = intersect(rownames(metal),rownames(metac))
  obj <- subset(celltypeobj, cells = samnames)
  absolute_path <- dirname(normalizePath(celltypeRDS))
  
  saveRDS(obj, file = paste0(absolute_path,"/SpLin_input.rds"))
}
















