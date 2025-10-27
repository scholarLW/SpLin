moduleGene = function(rds, rdsList = FALSE, dims = 10, resolution = 0.8, geneVector = NULL, ColumnCluster = FALSE, strict = TRUE, show_rownames = TRUE, show_colnames = TRUE, angle_col = "45", width = 18, height = 18, sv = 'V3')
{
  nst = 5  
  type = "unsigned"
  corType = "pearson"
  corFnc = ifelse(corType=="pearson", cor, bicor)
  maxPOutliers = ifelse(corType=="pearson",1,0.05)
  robustY = ifelse(corType=="pearson",T,F)
  MADmin = 0.01 # 0.01
  #probs = 0.75  # 0.75
  minModuleSize = 20
 	if (!(sv %in% c("V3", "V5", "v3", "v5"))) 
	{
		stop("sv must be either V3 or V5")
	} 
  suppressMessages(library(Seurat))
  suppressMessages(library(WGCNA))
  suppressMessages(library(dplyr))
  suppressMessages(library(ggplot2))	
  suppressMessages(library(stringr))	
  suppressMessages(library(RColorBrewer))
  suppressMessages(library(pheatmap))
  options(stringsAsFactors = FALSE)
  obj = readRDS(rds)
  colors <- c(brewer.pal(9, "Set1"), 
              "#FF34B3","#BC8F8F","#20B2AA","#00F5FF","#FFA500","#ADFF2F",
              "#FF6A6A","#7FFFD4", "#AB82FF","#90EE90","#00CD00","#008B8B",
              "#6495ED","#FFC1C1","#CD5C5C","#8B008B","#FF3030", "#7CFC00",
              '#1C386A','#5691B9','#D6E1FD','#087149','#51A186','#A6B5D6',
              '#90D4BB','#6B4470','#896F8C','#D7AAB1','#295286','#C0BDD2',
              '#BDC4DE','#3A5E92','#EAC49D','#566B32','#6E8557','#A7AA73',
              '#CED0A9','#D2C4A9','#7B3B1F','#A1725E','#C58F77','#F9C0A2',
              '#F4CCA8','#EEACB6','#D76C7C','#E45E52','#C84058','#99152C') 
  
  if(rdsList)
  {
    for (i in 1:length(obj)) {
      obj[[i]] <- SCTransform(obj[[i]], assay = "Spatial", return.only.var.genes = FALSE, verbose = FALSE)} 			
    features <- SelectIntegrationFeatures(object.list = obj, nfeatures = 3000)
    obj <- PrepSCTIntegration(object.list = obj, anchor.features = features)
    Em.anchors <- FindIntegrationAnchors(object.list = obj, normalization.method = "SCT", anchor.features = features)
    pbmc <- IntegrateData(anchorset = Em.anchors, normalization.method = "SCT")
  }else{
    pbmc <- SCTransform(obj, assay = "Spatial", new.assay.name = "SCT", verbose = FALSE)
  }
  pbmc <- RunPCA(object = pbmc, npcs = dims, verbose = FALSE)
  #pbmc <- RunUMAP(object = pbmc, reduction = "pca", dims = 1:dims)
  #pbmc = FindNeighbors(object = pbmc, k.param=10, dims = 1:dims)
  pbmc = FindNeighbors(object = pbmc, k.param=10, dims = 1:dims, reduction = "pca")
  subobj <- FindClusters(object = pbmc, resolution = resolution)	
  
  f <- FindAllMarkers(subobj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  ALLUSE = TRUE
  if(nrow(f)>0)
  {
    if(strict)
    {
      ix = f$p_val_adj < 0.05 & abs(f$avg_log2FC) > 0.5
      if(sum(ix)>0)
      {
        f = data.frame(f[ix,])
      }
    }
    genes = unique(f$gene) 
    if(length(genes)>=10)
    {
      ALLUSE = FALSE
    }
  }
  if(!is.null(geneVector))
  {
    if(!is.vector(geneVector))
    {
      stop("Error: The geneVector parameter should be a vector composed of gene names. Please re-enter.")
    }
    genes = unique(c(genes, geneVector))
    genes = genes[genes %in% rownames(subobj)]
    if(length(genes)>=10)
    {
      ALLUSE = FALSE
    }		
  }	
  if(ALLUSE)
  {	
    genes = rownames(subobj)
  }
  
  absolute_path <- paste0(dirname(normalizePath(rds)), "/Modules")
  dir.create(absolute_path, showWarnings = FALSE, recursive = TRUE)
  enableWGCNAThreads(nst)
  if(sv %in%  c("V5", "v5"))
	{
		dataM <- LayerData(object = subobj, assay = "SCT", layer = "data")
	}else{
		  dataM = subobj@assays$SCT@data  # counts	
	}
  dataf = data.frame(dataM[genes, ])
  col_sums <- colSums(dataf)  
  nonzero_cols <- which(col_sums != 0)  
  dataf <- dataf[, nonzero_cols]
  binGroup = data.frame(cbind(rownames(subobj@meta.data), subobj@meta.data$seurat_clusters))
  nsg = unique(binGroup$X2)
  dataExprraw = t(apply(dataf, 1, getCfreq, binGroup, nsg, colnames(dataf)))
  if(nrow(dataExprraw) == 1)
  {
    dataExprraw = t(dataExprraw)
  }
  colnames(dataExprraw) = nsg
  if(ncol(dataExprraw)>1)
  {
    m.mad <- apply(dataExprraw, 1, mad)
    #dataExprVar <- dataExprraw[which(m.mad > max(quantile(m.mad, probs=seq(0, 1, 1-probs))[2], MADmin)),]
    ix = which(m.mad > MADmin)
  }else{
    ix = which(dataExprraw > MADmin)
  }
  if(length(ix)>0)
  {
    if(ncol(dataExprraw) == 1)
    {
      dataExprVar <- t(t(dataExprraw[ix,]))
      colnames(dataExprVar) <- colnames(dataExprraw)
    }else{
      dataExprVar <- dataExprraw[ix,]
    }
    geneMDir = paste0(absolute_path, "/Gene")
    dir.create(geneMDir, showWarnings = FALSE, recursive = TRUE)
    dataExpr <- as.data.frame(t(dataExprVar))
    gsg = try(goodSamplesGenes(dataExpr, verbose = 3), silent = TRUE)
    if(inherits(gsg, "try-error"))
    {
      stop("Too few genes with valid expression levels in the required number of samples. Please reset the window parameter. The smaller the window, the more information you can obtain, but this does not necessarily mean that smaller is better. Choose an appropriate threshold. It is recommended to select 500.")
    }
    if (!gsg$allOK){
      # Optionally, print the gene and sample names that were removed:
      if (sum(!gsg$goodGenes)>0)
        printFlush(paste("Removing genes:",
                         paste(names(dataExpr)[!gsg$goodGenes], collapse = ",")));
      if (sum(!gsg$goodSamples)>0)
        printFlush(paste("Removing samples:",
                         paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ",")));
      # Remove the offending genes and samples from the data:
      dataExpr = dataExpr[gsg$goodSamples, gsg$goodGenes]
    }
    
    nGenes = ncol(dataExpr)
    nSamples = nrow(dataExpr)
    
    powers = c(c(1:10), seq(from = 12, to=30, by=2))
    sft = pickSoftThreshold(dataExpr, powerVector=powers,
                            networkType=type, verbose=5)
    power = sft$powerEstimate
    
    if(is.na(power))
    {
      power = ifelse(nSamples<20, ifelse(type == "unsigned", 9, 18),
                     ifelse(nSamples<30, ifelse(type == "unsigned", 8, 16),
                            ifelse(nSamples<40, ifelse(type == "unsigned", 7, 14),
                                   ifelse(type == "unsigned", 6, 12))      
                     )
      )
    }
    
    net = blockwiseModules(dataExpr, power = power, maxBlockSize = nGenes,
                           TOMType = type, minModuleSize = minModuleSize,
                           reassignThreshold = 0, mergeCutHeight = 0.25,
                           numericLabels = TRUE, pamRespectsDendro = FALSE,
                           saveTOMs=TRUE, corType = corType,
                           maxPOutliers=maxPOutliers, loadTOMs=TRUE,
                           saveTOMFileBase = paste0(geneMDir, "/Modules"),
                           verbose = 3)
    
    moduleLabels = net$colors
    moduleColors = labels2colors(moduleLabels)
    tmp = as.data.frame(cbind(colnames(t(moduleLabels)), moduleLabels, moduleColors))
    colnames(tmp) = c('geneName', 'moduleLabel', 'color')
    rownames(tmp) = NULL
    for(ikg in 1:nrow(tmp))
    {
      tmp$geneName[ikg] = strsplit(tmp$geneName[ikg], '\\.')[[1]][1]
    }
    if(0 %in% unique(as.numeric(tmp$moduleLabel)))
    {
      tmp$moduleLabel = as.numeric(tmp$moduleLabel) + 1
    }
    nCluster <- length(unique(tmp$moduleLabel))
    xend <- ceiling(nCluster/15)
    xeach <- ceiling(nCluster/xend)
    xlist <- sort(rep(c(1:xend), each=xeach, length.out = nCluster))
    df_xlist <- table(xlist) %>% as.data.frame()
    ylist <- rep(c(max(df_xlist$Freq):1),nrow(df_xlist), length.out=nCluster)
    
    df_legend <- as.data.frame(unique(tmp[, 2:3]))
    df_legend$x <- xlist
    df_legend$y <- ylist
    df_legend$moduleLabel = factor(df_legend$moduleLabel, levels = df_legend$moduleLabel)
    p_lgd <- ggplot(df_legend, aes(x=x, y=y))+
      geom_point(alpha=1, size=12, aes(colour=as.factor(moduleLabel)), pch = 15)+
      labs(x='',y='',title="")+
      geom_text(aes(x=x, y=y,label=moduleLabel), fontface="bold",color = "white",size = 4)+
      #geom_text(aes(x=x+xend*0.03, y=y,label=moduleLabel), fontface="bold",color = "black",size = 5, hjust = 0)+
      theme_void()+
      scale_x_continuous(limits = c(1, xend+1))+
      # scale_y_continuous(limits = c(ymin-2, ymax))+
      theme(axis.text = element_blank())+
      theme(axis.ticks = element_blank())+
      theme(panel.border = element_blank())+
      theme(axis.title = element_blank())+
      scale_color_manual(values=df_legend$color)+
      theme(legend.position = "none")
    ggsave(paste0(geneMDir, "/legend_moduleLabel.pdf"),p_lgd,device='pdf',width = 5, height = 4)
    load(net$TOMFiles[1], verbose=T)  
    ## Loading objects:
    ##   TOM
    
    TOM <- as.matrix(TOM)	  
    dissTOM = 1-TOM
    # Transform dissTOM with a power to make moderately strong
    # connections more visible in the heatmap
    plotTOM = dissTOM^7
    # Set diagonal to NA for a nicer plot
    diag(plotTOM) = NA
    # Call the plot function
    myheatcol <- colorRampPalette(c("red", "orange", "lemonchiffon"))(250) 
    df <- data.frame(  
      x = rep(0.5, 250),  
      y = seq(0, 1 - 1/250, length.out = 250),  
      width = rep(1, 250),  
      height = rep(1/250, 250),  
      fill = myheatcol  
    )  
    
    pyanse <- ggplot(df, aes(x = x, y = y, width = width, height = height, fill = fill)) +  
      geom_tile() +  
      scale_fill_identity() +  
      theme_void() +  
      labs(title = "Score")
    ggsave(paste0(geneMDir, "/colorRamp.pdf"),pyanse,device='pdf',width = 4, height = 7)
    
    geneTree = hclust(as.dist(dissTOM), method = "average")
    geneTree$labels = names(moduleLabels)
    pdf(file = paste0(geneMDir, "/geneTree.pdf"), width = 36, height = 6)
    plot(geneTree, main = "Gene Tree", sub = "", xlab = "Genes")
    dev.off()
    tmp$geneTree = geneTree$order
    tmp = tmp[order(as.numeric(tmp$geneTree)), ]
    #if(0 %in% unique(as.numeric(tmp$moduleLabel)))
    #{
    #	tmp$moduleLabel = as.numeric(tmp$moduleLabel) + 1
    #}
    write.table(tmp, file = paste0(geneMDir, "/moduleLabel.txt"), col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
    
    
    pdf(file=paste0(geneMDir, "/TOM_geneTree.pdf"))
    TOMplot(plotTOM, geneTree, col=myheatcol, moduleColors,
            main = paste0("Network heatmap plot"))
    dev.off()	
    
    objtmp = subset(subobj, features = tmp$geneName)
    geneMDir = paste0(absolute_path, "/windowCell")
    dir.create(geneMDir, showWarnings = FALSE, recursive = TRUE)
    tmpW = as.data.frame(cbind(rownames(objtmp@meta.data), objtmp@meta.data))
    colnames(tmpW)[1] = 'WindowCell'
    tmpW = tmpW[, c('WindowCell', 'seurat_clusters')]
    nc = nlevels(tmpW$seurat_clusters)
    colorsi = colors[1:nc]
    rownames(tmpW) = NULL
    for(ikg in 1:nrow(tmpW))
    {
      tmpW$color[ikg] = colorsi[as.numeric(as.character(tmpW$seurat_clusters[ikg])) + 1]
    }
    colnames(tmpW) = c('WindowCell', 'moduleLabel', 'color')
    write.table(tmpW, file = paste0(geneMDir, "/moduleLabel.txt"), col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

    if(sv %in%  c("V5", "v5"))
    {
      dataExpr <- as.data.frame(LayerData(object = objtmp, assay = "SCT", layer = "scale.data"))
    }else{
      dataExpr <- as.data.frame(objtmp@assays$SCT@scale.data) # counts data scale.data
    }  
    if(ColumnCluster)
    {
      tmpW = tmpW[order(tmpW$moduleLabel), ]
      dataExpr = dataExpr[, tmpW$WindowCell[tmpW$WindowCell %in% colnames(dataExpr)]]
    }
    
    unique_tmpW <- unique(tmpW[, c("moduleLabel", "color")])
    unique_tmpW <- unique_tmpW[order(as.numeric(unique_tmpW$moduleLabel)), ]
    unique_tmpW$moduleLabel <- as.character(unique_tmpW$moduleLabel)
    windowCellModule_colors <- setNames(unique_tmpW$color, unique_tmpW$moduleLabel)
    tmpp = tmp[tmp$geneName %in% rownames(dataExpr), ]
    unique_tmp <- unique(tmpp[, c("moduleLabel", "color")])
    unique_tmp <- unique_tmp[order(as.numeric(unique_tmp$moduleLabel)), ]
    unique_tmp$moduleLabel <- as.character(unique_tmp$moduleLabel)
    GeneModule_colors <- setNames(unique_tmp$color, unique_tmp$moduleLabel)
    ann_colors <- list(windowCellModule = windowCellModule_colors, GeneModule = GeneModule_colors)
    annotation_col = data.frame(windowCellModule = factor(tmpW$moduleLabel, levels = sort(unique(as.numeric(as.character(tmpW$moduleLabel))))))
    rownames(annotation_col) = tmpW$WindowCell		
    annotation_row = data.frame(GeneModule = factor(tmpp$moduleLabel, levels = sort(unique(as.numeric(as.character(tmpp$moduleLabel))))))
    rownames(annotation_row) = tmpp$geneName
    
    pdf(paste0(geneMDir,"/Heatmap.pdf"), width = width, height = height)
    pheatmap(as.matrix(dataExpr), clustering_distance_rows = "correlation", clustering_method = "average", cluster_rows = TRUE, cluster_cols = FALSE, annotation_col = annotation_col, annotation_row = annotation_row, show_rownames = show_rownames, show_colnames = show_colnames, show_colnames_annotation = FALSE, show_rownames_annotation = FALSE, annotation_colors = ann_colors, angle_col = angle_col)
    dev.off()
  }
}