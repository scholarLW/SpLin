SpatialFeaturePlot_lw <- function(vln.df, SpatialFeaturePlotpointsize, SpatialFeaturePlotpointalpha, SpatialFeaturePlotstriptextsize, SpatialFeaturePlotmycolor)
{
  p <- ggplot(data = vln.df, aes(x = x, y = y, colour=Expression)) +
    geom_point(size=SpatialFeaturePlotpointsize, alpha = SpatialFeaturePlotpointalpha, shape=19) + 
    scale_color_gradientn(colours=rev(SpatialFeaturePlotmycolor),
                          breaks=c(0,floor(max(vln.df$Expression))),
                          labels=c("Low","High")) +
    theme_bw() +
    theme(
      strip.background = element_rect(fill=NA, color=NA),
      strip.text = element_text(size=SpatialFeaturePlotstriptextsize),		
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      panel.border = element_blank(),		
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  return(p)
}

SpatialFeaturePlot_lwCL <- function(vln.df, y_min = NA, y_max = NA, SpatialFeaturePlotpointsize, SpatialFeaturePlotpointalpha, SpatialFeaturePlotstriptextsize, SpatialFeaturePlotmycolor)
{  
  if (is.na(y_min)) {  
    y_min <- min(vln.df$y)
  }  
  
  if (is.na(y_max)) {  
    y_max <- 0.2 * max(vln.df$x)
  }  
  
  p <- ggplot(data = vln.df, aes(x = x, y = y, colour = Expression)) +  
    geom_point(size = SpatialFeaturePlotpointsize, alpha = SpatialFeaturePlotpointalpha, shape = 19) +   
    scale_color_gradientn(colours = rev(SpatialFeaturePlotmycolor),  
                          breaks = c(0, floor(max(vln.df$Expression))),  
                          labels = c("Low", "High")) +  
    scale_y_continuous(limits = c(y_min, y_max)) +  
    theme_bw() +  
    theme(  
      legend.position = "none", 
      strip.background = element_rect(fill = NA, color = NA),  
      strip.text = element_text(size = SpatialFeaturePlotstriptextsize),		  
      panel.grid.major = element_blank(),   
      panel.grid.minor = element_blank(),  
      panel.background = element_blank(),  
      panel.border = element_blank(),		  
      axis.title.x = element_blank(),  
      axis.title.y = element_blank(),  
      axis.text.x = element_blank(),  
      axis.ticks.x = element_blank(),  
      axis.text.y = element_blank(),  
      axis.ticks.y = element_blank()  
    )  
  return(p)  
}


markPlot = function(rds, markerGenes = c('Reg3b','Dcn','Hmgcs2'), SpatialFeaturePlotstriptextsize = 20, SpatialFeaturePlotpointsize = 0.1, SpatialFeaturePlotpointalpha = 0.5, pdf.width = 16, pdf.height = 4, png.width = 50, png.height = 8, png.dpi = 300, sv = 'V3')
{
  suppressMessages(library(Seurat))
  suppressMessages(library(ggplot2))
  suppressMessages(library(dplyr))
  suppressMessages(library(reshape2))
 	if (!(sv %in% c("V3", "V5", "v3", "v5"))) 
	{
		stop("sv must be either V3 or V5")
	} 
  SpatialFeaturePlotmycolor =  rev(c('grey90', 'grey100', "lavenderblush", "lightcoral", "red","red4"))
  seuratobj <- readRDS(rds)
  absolute_path <- paste0(dirname(normalizePath(rds)),"/markerPlot")
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
    da = da[!ix, c('x','y',gene)]
    vln.df = da %>% reshape2::melt(, gene)
    vln.df = data.frame(vln.df[, c('x','y','variable','value')])
    colnames(vln.df) = c("x","y","gene","Expression")
    p <- SpatialFeaturePlot_lw(vln.df, SpatialFeaturePlotpointsize, SpatialFeaturePlotpointalpha, SpatialFeaturePlotstriptextsize, SpatialFeaturePlotmycolor)
    ggsave(paste0(absolute_path,"/",gene,".ST.pdf"),p)
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
    p <- SpatialFeaturePlot_lwCL(vln.df, y_min = NA, y_max = NA, SpatialFeaturePlotpointsize, SpatialFeaturePlotpointalpha, SpatialFeaturePlotstriptextsize, SpatialFeaturePlotmycolor)
    ggsave(paste0(absolute_path,"/",gene,".ST_CL.pdf"),p,limitsize = FALSE, width = pdf.width, height = pdf.height)	
    ggsave(filename = paste0(absolute_path,"/",gene,".ST_CL.png"), plot = p, limitsize = FALSE, device = "png", width = png.width, height = png.height, units = "in", dpi = png.dpi)
  }
}