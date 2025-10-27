signalToFig = function(rds, colN = NULL, pointcellsize = NULL, dpi = 300)
{
  suppressMessages(library(Seurat))
  suppressMessages(library(RColorBrewer))
  obj = readRDS(rds)
  
  absolute_path <- dirname(normalizePath(rds))
  colors <- c(brewer.pal(9, "Set1"), 
              "#FF34B3","#BC8F8F","#20B2AA","#00F5FF","#FFA500","#ADFF2F",
              "#FF6A6A","#7FFFD4", "#AB82FF","#90EE90","#00CD00","#008B8B",
              "#6495ED","#FFC1C1","#CD5C5C","#8B008B","#FF3030", "#7CFC00",
              '#1C386A','#5691B9','#D6E1FD','#087149','#51A186','#A6B5D6',
              '#90D4BB','#6B4470','#896F8C','#D7AAB1','#295286','#C0BDD2',
              '#BDC4DE','#3A5E92','#EAC49D','#566B32','#6E8557','#A7AA73',
              '#CED0A9','#D2C4A9','#7B3B1F','#A1725E','#C58F77','#F9C0A2',
              '#F4CCA8','#EEACB6','#D76C7C','#E45E52','#C84058','#99152C') 
  
  if(is.null(colN))
  {
    color = adjustcolor("black", alpha.f = 0.8) 
  }else{
    if(!colN %in% colnames(obj@meta.data))
    {
      colN = NULL
      cat("Warning: The provided column name does not exist in the metadata of the Seurat object file. \n")
    }else{	
      n_groups <- length(unique(obj@meta.data[, colN]))
      group_colors <- rainbow(n_groups)
      color <- group_colors[match(obj@meta.data[, colN], unique(obj@meta.data[, colN]))]
    }
  }
  tmpn = max(max(obj@meta.data$x), max(obj@meta.data$y))
  png(
    filename = paste0(absolute_path, '/signalMatrix.png'),
    width = tmpn,    
    height = tmpn,
    units = "px",     
    res = dpi,             
    type = "cairo"         
  )
  if(is.null(pointcellsize))
  {
    pointcellsize <- 0.0005*tmpn
  }
  par(mar = c(0, 0, 0, 0)) 
  plot(
    obj@meta.data$x, 
    obj@meta.data$y,
    axes = F, 
    xlab = '', 
    ylab = '', 
    pch = 20,              
    cex = pointcellsize,             
    col = color
  )
  dev.off()
}
















