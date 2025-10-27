RapidRegisterHE = function(pixelssDNAFile, pixelHEFile, markingPointsImage, markingPoints_MASKJSON, cellPointsImage, cellPoints_MASKJSON, idjson, transform = NULL, kpixel = 10, epsilon = 0.001, interval = 0.05, intervalAngle = 30, up = 0, down = 0, left = 0, right = 0, theta = 0)
{
  thetaStep = 0.5  
  miuStep = 0.01 
  ngrid = 10 
  labelfont = 1
  labelfontsize = 2
  axistitlefontsize = 3
  spatialpointsize = 0.01
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
  suppressMessages(library(sp))	
  suppressMessages(library(dplyr))
  
  idjsonData = fromJSON(idjson)
  ow = idjsonData$original_dimensions$width
  oh = idjsonData$original_dimensions$height
  absolute_path <- dirname(normalizePath(pixelHEFile))
  JsonFile = paste0(absolute_path, '/AllRegistrationSchemes.json')
  script <- system.file("python/rapidRegPixel.py", package = "SpLin")
  
  cmd <- paste0("python ", script, " -o ", JsonFile, " -p ", pixelHEFile, " -mpi ", markingPointsImage, " -mpj ", markingPoints_MASKJSON, " -cpi ", cellPointsImage, " -cpj ", cellPoints_MASKJSON, " -s ", ow, " ", oh, " -th ", thetaStep, " -m ", miuStep, " -n ", ngrid, " -k ", kpixel, " -e ", epsilon, " --up ", up, " --down ", down, " --left ", left, " --right ", right, " --theta ", theta)
  
  if (!is.null(transform)) 
  {
    cmd <- paste0(cmd, " -t ", shQuote(transform))
  } else {
    cmd <- cmd
  }
  system(cmd)	
  
  rmFile = paste0(absolute_path, '/HE_pixel_MASK_update.png')
  if (file.exists(rmFile))
  {
    file.remove(rmFile)
  }	
  
  ssDNAda = as.data.frame(fread(pixelssDNAFile))
  colnames(ssDNAda)[1:2] = c('x', 'y')
  jsonData = fromJSON(JsonFile)
  if(length(jsonData) == 0)
  {
    stop("No registration scheme meets the requirements, please reset the threshold parameter minRC!")
  }
  
  m = length(jsonData$RC)
  for(index in 1:m)
  {
    ssDNAHEFile = paste0(absolute_path, '/Scheme/', index, "/ssDNAHE_Adjusted_output_coordinates_and_colors.txt")
    png_file_path <- paste0(absolute_path, '/Scheme/', index, "/Multiple_AutoRegister_ST_plot.png")
    ithemp <- jsonData$MASKPOINT[[index]]		
    colnames(ithemp) = c('x', 'y')
    bnd <- as.data.frame(rbind(ithemp, ithemp[1, ]))
    ratio <- round(jsonData$RC[index], 4)
    pards <- round(jsonData$PARDS[index], 4)
    inside <- point.in.polygon(ssDNAda$x, ssDNAda$y, bnd$x, bnd$y)
    freq = round(sum(inside)/length(inside),4) * 100			
    
    HEda <- as.data.frame(fread(ssDNAHEFile))
    colnames(HEda)[1:2] = c('x', 'y')
    all_data <- rbind(
      transform(HEda, dataset = "HE"),
      transform(ssDNAda, dataset = "ssDNA")
    )
    all_data$dataset = factor(all_data$dataset, levels = c('ssDNA', 'HE'))
    df <- all_data[order(all_data$dataset), ]
    color = colorRampPalette(brewer.pal(9, "Set1"))(2)
    names(color) = c('ssDNA', 'HE')
    
    
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
                      intervalAngle = intervalAngle)			
    
    ggsave(filename = png_file_path, plot = p, device = "png", width = 10, height = 8, units = "in", dpi = 300)
  }	
}

