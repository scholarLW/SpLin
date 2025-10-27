wcellDue = function(rds, Assay = 'Spatial', method = NULL, windowSV = 250, windowSH = 10)
{
  suppressMessages(library(Seurat))
  options(stringsAsFactors = FALSE)
  obj = readRDS(rds)
  absolute_path <- paste0(dirname(normalizePath(rds)),"/IGSI/winCell")
  dir.create(absolute_path, showWarnings = FALSE, recursive = TRUE, mode = "0777")
  window_size <- windowSV	
  if(is.null(method))
  {
    method <- 'Vertical' 
  }
  metadata = obj@meta.data
  metadata$window = NA
  metadata$windowVertical = NA
  metadata$windowHorizontal = NA
  
  if(windowSV < 100 || windowSV > 1000)
  {
    stop("The vertical window size has exceeded the specified range. The standard range is from 100 to 1000. Please select a value within this range, and it must be an integer.")
  }
  if(windowSH < 10 || windowSH > 100)
  {
    stop("The horizontal window size has exceeded the specified range. The standard range is from 10 to 100. Please select a value within this range, and it must be an integer.")
  }	
  
  maxWindow <- 1000
  metadata = metadata[order(metadata$xnew), ]
  max_x = max(max(metadata$xnew), max(metadata$xnew) - min(metadata$xnew))
  num_windowsVertical = max(ceiling(max_x / windowSV), 50) 
  
  for (i in 1:num_windowsVertical)
  {
    start_index = min(metadata$xnew) + (i - 1) * windowSV - 1
    end_index = min(metadata$xnew) + i * windowSV + 1
    metadata$windowVertical[(metadata$xnew >= start_index) & (metadata$xnew <= end_index)] = i
  }
  if(end_index < max_x + 1)
  {
    metadata$windowVertical[metadata$xnew >= end_index] = i
  }
  if(method %in% c('Horizontal', 'horizontal', 'h', 'H'))
  {
    window_size <- windowSH	
    maxWindow <- 100
    metadata = metadata[order(metadata$ynew), ]
    max_y = max(max(metadata$ynew), max(metadata$ynew) - min(metadata$ynew))
    num_windowsHorizontal = max(ceiling(max_y / windowSH), 20) 
    windowVertical = unique(metadata$windowVertical)
    dataTMP = NULL
    for(wd in windowVertical)
    {
      ix = metadata$windowVertical == wd
      tmpDa = metadata[ix,]
      ymin = min(tmpDa$ynew)
      ymax = max(tmpDa$ynew)
      window_sizey = (ymax - ymin)/num_windowsHorizontal
      for (i in 1:num_windowsHorizontal)
      {
        start_index = min(tmpDa$ynew) + (i - 1) * window_sizey - 1
        end_index = min(tmpDa$ynew) + i * window_sizey + 1
        tmpDa$windowHorizontal[(tmpDa$ynew >= start_index) & (tmpDa$ynew <= end_index)] = i
      }
      if(end_index < max_y + 1)
      {
        tmpDa$windowHorizontal[tmpDa$ynew >= end_index] = i
      }	
      dataTMP = as.data.frame(rbind(dataTMP, tmpDa))
    }
    dataTMP = dataTMP[rownames(metadata), ]
    metadata$window = dataTMP$windowHorizontal
  }else{
    metadata$window = metadata$windowVertical
  }
  obj@meta.data = metadata[rownames(obj@meta.data), ]
  count_matrix <- GetAssayData(obj, assay = Assay, slot = "counts")
  # Mean
  #data = apply(count_matrix, 1, mean_f, colnames(count_matrix), metadata)
  # GFai
  data = apply(count_matrix, 1, calculate_all_GFai, colnames(count_matrix), metadata, window_size, maxWindow = maxWindow)	
  rownames(data) = paste("Window", 1:nrow(data), sep = "")
  data = t(data)
  objwd <- CreateSeuratObject(counts = data, assay = "Spatial", min.cells = 1, min.features = 1)
  saveRDS(objwd, file = paste0(absolute_path, "/windowCell.RDS"))
  saveRDS(obj, file = paste0(absolute_path, "/SpLin_EXP_output_withWindow.RDS"))
}
