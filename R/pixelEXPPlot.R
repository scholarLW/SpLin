SpatialFeaturePlot <- function(vln.df, vln.df1, SpatialFeaturePlotpointsize, SpatialFeaturePlotpointalpha, SpatialFeaturePlotmycolor, y_min = NA, y_max = NA) 
{  
  if (is.na(y_min)) {  
    y_min <- min(vln.df$y)
  }  
  
  if (is.na(y_max)) {  
    y_max <- max(0.05 * max(vln.df$x), max(vln.df$y, vln.df1$ynew))
  }  
  color_mapping <- setNames(vln.df1$Cell,vln.df1$Cell)
  
  p1 <- ggplot(data = vln.df1, aes(x = xnew, y = ynew, color = Cell)) +
    geom_point(size = 0.01, alpha = 1, shape = 20) +
    scale_y_continuous(limits = c(y_min, y_max)) +
    scale_color_manual(name = "Cell", values = color_mapping) +
    theme_bw() +
    xlab('') +
    ylab('Pixel') + 
    theme(
      legend.position = "none",
      strip.background = element_rect(fill = NA, color = NA),
      strip.text = element_text(),
      panel.grid = element_blank(),  
      plot.background = element_blank(),
      axis.line = element_line(color = "black"), 
      axis.text = element_text(),
      axis.text.y = element_blank(),
      axis.ticks = element_line(color = "black"), 
      axis.ticks.y = element_blank(), 
      panel.grid.major.x = element_blank(),  
      panel.grid.minor.x = element_blank()   
    )
  
  p2 <- ggplot(data = vln.df, aes(x = x, y = y, colour = Expression)) +  
    geom_point(size = SpatialFeaturePlotpointsize, alpha = SpatialFeaturePlotpointalpha, shape = 16) +   
    scale_color_gradientn(colours = rev(SpatialFeaturePlotmycolor),  
                          breaks = c(0, floor(max(vln.df$Expression))),  
                          labels = c("Low", "High")) +  
    scale_y_continuous(limits = c(y_min, y_max)) +  
    theme_bw() +  
    xlab('') +
    ylab('Exp.') + 
    theme(  
      legend.position = "top",
      strip.background = element_rect(fill = NA, color = NA),
      strip.text = element_text(),
      panel.grid = element_blank(),  
      plot.background = element_blank(),
      axis.line = element_line(color = "black"), 
      axis.text = element_text(), 
      axis.text.y = element_blank(),
      axis.ticks = element_line(color = "black"),
      axis.ticks.y = element_blank(), 
      panel.grid.major.x = element_blank(), 
      panel.grid.minor.x = element_blank()   
    )  
  
  return(p2/p1)  
}

hex_to_rgb <- function(hex_color) {
  rgb_color <- col2rgb(hex_color)
  return(as.vector(rgb_color))
}

is_near_black <- function(x, tolerance = 10) 
{
  hex_color = x[1]
  rgb_color <- hex_to_rgb(hex_color)
  all(rgb_color <= tolerance)
}

is_near_white <- function(x, tolerance = 10)
{
  hex_color = x[1]
  rgb_color <- hex_to_rgb(hex_color)
  return(all(rgb_color >= 255 - tolerance))
}

pixelEXPPlot = function(rds, pixelPointsFile, markerGenes = c('Reg3b','Dcn','Hmgcs2'), imageType = 'HE', tolerance = 10, xmin = NULL, xmax = NULL, 
                        SpatialFeaturePlotstriptextsize = 20, SpatialFeaturePlotpointsize = 0.5, SpatialFeaturePlotpointalpha = 1, sv = 'V3')
{
  suppressMessages(library(Seurat))
  suppressMessages(library(ggplot2))
  suppressMessages(library(dplyr))
  suppressMessages(library(reshape2))
  suppressMessages(library(patchwork))
  suppressMessages(library(data.table))
 	if (!(sv %in% c("V3", "V5", "v3", "v5"))) 
	{
		stop("sv must be either V3 or V5")
	} 
  SpatialFeaturePlotmycolor =  rev(c('grey90', 'grey100', "lavenderblush", "lightcoral", "red","red4"))
  seuratobj <- readRDS(rds)
  absolute_path <- paste0(dirname(normalizePath(pixelPointsFile)),"/markerPlot")
  dir.create(absolute_path, showWarnings = FALSE, recursive = TRUE, mode = "0777")
	if(sv %in%  c("V5", "v5"))
	{
		counts_matrix <- LayerData(object = seuratobj, assay = "Spatial", layer = "data")
	}else{
		counts_matrix <- seuratobj@assays$Spatial@data
	}  
  for(gene in markerGenes)
  {
    tmp = data.frame(counts_matrix[gene, ])
    da = merge(seuratobj@meta.data, tmp, by="row.names", all=F)
    colnames(da)[ncol(da)] = gene
    ix = is.na(da$bin)
    da = da[!ix, c('xnew','ynew',gene)]
    colnames(da)[1:2] = c('x','y')
    da = da[order(da$x), ]
    vln.df = da %>% reshape2::melt(, gene)
    vln.df = data.frame(vln.df[, c('x','y','variable','value')])
    colnames(vln.df) = c("x","y","gene","Expression")
    if(!is.null(xmin))
    {
      ix = vln.df$x >= xmin
      vln.df = vln.df[ix,]
    }
    if(!is.null(xmax))
    {
      ix = vln.df$x <= xmax
      vln.df = vln.df[ix,]
    }
    
    data = as.data.frame(fread(pixelPointsFile))
    if(!is.null(xmin))
    {
      ix = data$xnew >= xmin
      data = data[ix,]
    }
    if(!is.null(xmax))
    {
      ix = data$xnew <= xmax
      data = data[ix,]
    }
    if(imageType %in% 'HE')
    {
      isX = apply(data, 1, is_near_white, tolerance)
    }else{
      isX = apply(data, 1, is_near_black, tolerance)
    }
    data = data[!isX, ]
    
    p <- SpatialFeaturePlot(vln.df, data, SpatialFeaturePlotpointsize, SpatialFeaturePlotpointalpha, SpatialFeaturePlotmycolor)
    
    if(!is.null(xmax) || !is.null(xmin))
    {
      ggsave(paste0(absolute_path,"/",gene,".ST_Pixel_CL.pdf"), p, limitsize = FALSE, width = 8, height = 4)
      ggsave(filename = paste0(absolute_path,"/",gene,".ST_Pixel_CL.png"), plot = p, limitsize = FALSE, device = "png", width = 8, height = 4, units = "in", dpi = 300)
    }else{
      ggsave(paste0(absolute_path,"/",gene,".ST_Pixel_CL_allregion.pdf"), p, limitsize = FALSE, width = 20, height = 8)	
      ggsave(filename = paste0(absolute_path,"/",gene,".ST_Pixel_CL_allregion.png"), plot = p, limitsize = FALSE, device = "png", width = 20, height = 8, units = "in", dpi = 300)
    }
  }
}