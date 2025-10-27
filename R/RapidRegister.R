RapidRegister = function(json, cellPointsFile, markingPointsImage, markingPoints_MASKJSON, cellPointsImage, cellPoints_MASKJSON, idjson, transform = NULL, kpixel = 10, epsilon = 0.001, interval = 0.05, intervalAngle = 30, up = 0, down = 0, left = 0, right = 0, theta = 0, spatialpointsize = 0.1, Signal = FALSE)
{
  thetaStep = 0.5
  miuStep = 0.03
  ngrid = 1
  labelfont = 1
  labelfontsize = 2
  axistitlefontsize = 3
  #spatialpointsize = 0.1
  pointalpha = 1
  polygonalpha = 0.8
  pointshape = 16
  legendtextsize = 10
  legendtitlesize = 25
  guidelegendsize = 6
  
  suppressMessages(library(data.table))
  suppressMessages(library(ggplot2))
  suppressMessages(library(RColorBrewer))
  suppressMessages(library(jsonlite))
  suppressMessages(library(mgcv))
  suppressMessages(library(dplyr))	
  idjsonData = fromJSON(idjson)
  ow = idjsonData$original_dimensions$width
  oh = idjsonData$original_dimensions$height
  absolute_path <- dirname(normalizePath(json))
  JsonFile = paste0(absolute_path, '/AllRegistrationSchemes.json')
  script <- system.file("python/rapidRegEXP.py", package = "SpLin")
  cmd <- paste0("python ", script, " -j ", json, " -o ", JsonFile, " -cp ", cellPointsFile, " -mpi ", markingPointsImage, " -mpj ", markingPoints_MASKJSON, " -cpi ", cellPointsImage, " -cpj ", cellPoints_MASKJSON, " -s ", ow, " ", oh, " -th ", thetaStep, " -m ", miuStep, " -n ", ngrid, " -k ", kpixel, " -e ", epsilon, " --up ", up, " --down ", down, " --left ", left, " --right ", right, " --theta ", theta)
  
  if (!is.null(transform)) 
  {
    cmd <- paste0(cmd, " -t ", shQuote(transform))
  } else {
    cmd <- cmd
  }
  if (Signal) 
  {
    cmd <- paste0(cmd, " -S")
  } else {
    cmd <- cmd
  }	
  system(cmd)
  
  rmFile = paste0(absolute_path, '/markingPoints_MASK_update.png')
  if (file.exists(rmFile))
  {
    file.remove(rmFile)
  }	
  
  obj = as.data.frame(fread(cellPointsFile,header=F,sep="\t"))
  colnames(obj) = c('x', 'y')
  jsonData = fromJSON(JsonFile)
  if(length(jsonData) == 0)
  {
    stop("No registration scheme meets the requirements, please reset the threshold parameter minRC!")
  }	
  pdf_file_path <- paste0(absolute_path, '/Multiple_AutoRegister_ST_plot.pdf')
  pdf(file = pdf_file_path, width = 8, height = 6)
  
  m = length(jsonData$RC)
  for(index in 1:m)
  {
    tmp = obj			
    i = which(jsonData$JSON$shapes[[index]]$label == 'Wai')
    j = which(jsonData$JSON$shapes[[index]]$label == 'Nei')
    n = nrow(jsonData$JSON$shapes[[index]]$points[[j]])
    bnd = rbind(jsonData$JSON$shapes[[index]]$points[[i]], jsonData$JSON$shapes[[index]]$points[[j]][n:1,])			
    
    inside <- in.out(bnd, as.matrix(tmp))
    ratio <- round(jsonData$RC[index], 4)
    pards <- round(jsonData$PARDS[index], 4)
    freq = round(sum(inside)/length(inside),4) * 100
    
    tmp$seurat_clusters = 'Excluded'
    tmp$seurat_clusters[inside] = 'Included'
    color = colorRampPalette(brewer.pal(9, "Set1"))(2)
    names(color) = c('Included', 'Excluded')
    
    umapdata <- tmp
    umapdata$seurat_clusters = factor(umapdata$seurat_clusters, levels = c('Included','Excluded'))
    umapdata = umapdata[order(umapdata$seurat_clusters), ]
    
    
    p <- autoReplot(umapdata, 
                    clusters='seurat_clusters', 
                    color,
                    labelfontsize, 
                    axistitlefontsize, 
                    spatialpointsize, 
                    pointalpha, 
                    pointshape, 
                    legendtextsize, 
                    legendtitlesize, 
                    guidelegendsize,
                    ratio,
                    freq,
                    index = index,
                    pards,
                    jsonData$CENT[[index]][1],
                    jsonData$CENT[[index]][2],
                    interval = interval,
                    intervalAngle = intervalAngle)				
    print(p)
  }
  dev.off()
  
  absolute_path <- paste0(absolute_path, "/Scheme")
  for(index in 1:m)
  {		
    dir.create(paste0(absolute_path, '/', index), showWarnings = FALSE, recursive = TRUE)
    png_file_path <- paste0(absolute_path, '/', index, "/Multiple_AutoRegister_ST_plot.png")
    
    if(isTRUE(Signal))
    {
      i = which(jsonData$JSON$shapes[[index]]$label == 'Sign')
      bnd = jsonData$JSON$shapes[[index]]$points[[i]]
    }else{
      i = which(jsonData$JSON$shapes[[index]]$label == 'Wai')
      j = which(jsonData$JSON$shapes[[index]]$label == 'Nei')
      n = nrow(jsonData$JSON$shapes[[index]]$points[[j]])
      bnd = rbind(jsonData$JSON$shapes[[index]]$points[[i]], jsonData$JSON$shapes[[index]]$points[[j]][n:1,])			
    }			
    
    colnames(bnd) = c('x', 'y')			
    ratio <- round(jsonData$RC[index], 4)
    pards <- round(jsonData$PARDS[index], 4)
    inside <- in.out(bnd, as.matrix(obj))
    freq = round(sum(inside)/length(inside),4) * 100
    
    if(isTRUE(Signal))
    {
      i = which(jsonData$shapes[[index]]$label == 'Sign')
      all_data <- rbind(
        transform(bnd, dataset = "SignalPanel"),
        transform(obj, dataset = "Signal")
      )
      all_data$dataset = factor(all_data$dataset, levels = c('Signal', 'SignalPanel'))
      df <- all_data[order(all_data$dataset), ]
      color = colorRampPalette(brewer.pal(9, "Set1"))(2)
      names(color) = c('Signal', 'SignalPanel')		
      
    }else{
      all_data <- rbind(
        transform(bnd, dataset = "HE/ssDNA"),
        transform(obj, dataset = "EXP")
      )
      all_data$dataset = factor(all_data$dataset, levels = c('EXP', 'HE/ssDNA'))
      df <- all_data[order(all_data$dataset), ]
      color = colorRampPalette(brewer.pal(9, "Set1"))(2)
      names(color) = c('EXP', 'HE/ssDNA')		
    }
    
    
    p <- autoReplotHE(df, 
                      clusters=levels(df$dataset), 
                      color,
                      labelfontsize, 
                      axistitlefontsize, 
                      spatialpointsize, 
                      pointalpha,
                      polygonalpha,						
                      pointshape, 
                      legendtextsize, 
                      legendtitlesize, 
                      guidelegendsize,
                      ratio,
                      freq,
                      index = index,
                      pards,
                      jsonData$CENT[[index]][1],
                      jsonData$CENT[[index]][2],
                      interval = interval,
                      intervalAngle = intervalAngle,
                      Signal = Signal)			
    ggsave(filename = png_file_path, plot = p, device = "png", width = 10, height = 8, units = "in", dpi = 300)
  }		
}

