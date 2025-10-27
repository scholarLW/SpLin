IGSI = function(rds, GMF, topFreq = 0.05, nbin = 10, alpha = 0.001, SD = 3, topSegFreq = 0.4, TopCellFreq = 0, eps = 200, grid_density = 50, kneighbor = 10, corr = 0.5, piontsize = 0.5, sv = 'V3')
{
  suppressMessages(library(Seurat))
  suppressMessages(library(dplyr))
  suppressMessages(library(circlize))
  suppressMessages(library(reshape2))
  suppressMessages(library(stringr))
  suppressMessages(library(RColorBrewer))
  suppressMessages(library(ggplot2))
  suppressMessages(library(ggpubr))
  suppressMessages(library(ggsignif))
  suppressMessages(library(gghalves))
  suppressMessages(library(ggsci))
  suppressMessages(library(rlang))
  suppressMessages(library(data.table))
  suppressMessages(library(DNAcopy))
  suppressMessages(library(fuzzyjoin))
  suppressMessages(library(dbscan))
  suppressMessages(library(nls2))
  suppressMessages(library(ggdensity))	
  suppressMessages(library(akima))
  suppressMessages(library(spdep))
  suppressMessages(library(sf))
  suppressMessages(library(sp))
  suppressMessages(library(raster))
  suppressMessages(library(gstat))
  suppressMessages(library(ggquiver))
  suppressMessages(library(patchwork))
  
  color <- pal_npg("nrc")(2)
  options(stringsAsFactors = FALSE)
  obj = readRDS(rds)
  if (!(sv %in% c("V3", "V5", "v3", "v5"))) 
  {
    stop("sv must be either V3 or V5")
  }
  absolute_path <- paste0(dirname(normalizePath(rds)))
  absolute_path <- paste0(absolute_path,"/HGMSR_",topFreq,'_',nbin,'_',alpha,'_',SD,'_',topSegFreq,"_",TopCellFreq,"_",eps,"_",grid_density,"_",kneighbor,"_",corr)
  
  dir.create(absolute_path, showWarnings = FALSE, recursive = TRUE, mode = "0777")
  
  metadata = obj@meta.data
  metadata = metadata[order(metadata$xnew), ]
  metadata = metadata[colnames(obj), ]
  obj@meta.data = metadata
  
  tmp = as.data.frame(fread(GMF))
  tmp = tmp[order(as.numeric(tmp$moduleLabel)), ]
  moduleLabelNUM = unique(tmp$moduleLabel)
  geneSets = list()
  for(i in 1:length(moduleLabelNUM))
  {
    geneSets[[i]] = tmp$geneName[tmp$moduleLabel==moduleLabelNUM[i]]
  }
  
  objmodule =  AddModuleScore(object = obj, features = geneSets, name = 'ModuleScoreGeneSet', nbin = nbin)
  metadata = objmodule@meta.data
  metadata = metadata[order(metadata$xnew), ]
  write.table(metadata, file = paste0(absolute_path,"/metaData.txt"), quote = F, col.names=T, row.names=F, sep = "\t")
  
  module_columns <- startsWith(names(metadata), "ModuleScoreGeneSet")
  mycolor = c('#1874CD',rev(brewer.pal(11, "Spectral"))[-1]) 
  
	if(sv %in%  c("V5", "v5"))
	{
		counts_matrix <- LayerData(object = obj, assay = "Spatial", layer = "counts")
	}else{
		counts_matrix <- obj@assays$Spatial@counts
	}

  arrayMean = NULL  
  for (i in seq_along(module_columns))
  {
    if (module_columns[i])
    {	
      current_column <- names(metadata)[i]
      windowDa = metadata[, c('window', 'xnew', 'ynew', current_column, 'x', 'y')]
      colnames(windowDa)[4] = 'Score'	
      dataDD = cbind(rownames(windowDa), windowDa)
      colnames(dataDD)[1] = 'cell'		
      top_percent <- quantile(dataDD$Score, 1 - topFreq, na.rm = TRUE)		
      mu_top <- mean(dataDD$Score[dataDD$Score >= top_percent], na.rm = TRUE)
      sigma_top <- sd(dataDD$Score[dataDD$Score >= top_percent], na.rm = TRUE)	
      vln.dat <- dataDD %>%
        mutate(Group = ifelse(Score >= top_percent, "TopCell", "NonTopCell"))
      xmin = min(vln.dat$x)
      xmax = max(vln.dat$x)
      ymin = min(vln.dat$y)
      ymax = max(vln.dat$y)
      geneLIST = geneSets[i - sum(!module_columns)][[1]]
      geneLIST = geneLIST[geneLIST %in% rownames(counts_matrix)]
      if(length(geneLIST)>0)
      {		
        if(length(geneLIST) == 1)
        {
          expGeneList = t(counts_matrix[geneLIST, ])
        }else{
          expGeneList = counts_matrix[geneLIST, ]
        }
        result <- expGeneList %>%
          as.matrix() %>%
          {data.frame(
            Sum = colSums(.),
            Mean = colMeans(.)
          )} %>%
          tibble::rownames_to_column(., "Cell") %>%
          setNames(c("Cell", "Sum", "Mean"))
        rownames(result) = result$Cell
        result = result[rownames(metadata), ]
        result = as.data.frame(cbind(result, metadata[,c('x','y')]))
        p1 <- ggplot(result, aes(x, y,color=Sum))+
          geom_point(size = piontsize, alpha = 0.9, shape = 15)+
          scale_color_gradientn('Sum Expression', colours = mycolor)+
          scale_x_continuous(limits = c(xmin-2, xmax))+
          scale_y_continuous(limits = c(ymin-2, ymax))+
          theme(axis.text = element_blank())+
          theme(axis.ticks = element_blank())+
          theme(panel.border = element_blank())+
          theme(axis.title = element_blank())+
          theme(panel.background = element_blank(),
                panel.border = element_blank(), 
                legend.background = element_blank(),
                legend.position = "right",
                legend.text = element_text(size=16),
                legend.title=element_text(size = 15, face ="bold"),		
                legend.key.size = unit(0.25, "inches"))
        
        p2 <- ggplot(result, aes(x, y,color=Mean))+
          geom_point(size = piontsize, alpha = 0.9, shape = 15)+
          scale_color_gradientn('Mean Expression', colours = mycolor)+
          scale_x_continuous(limits = c(xmin-2, xmax))+
          scale_y_continuous(limits = c(ymin-2, ymax))+
          theme(axis.text = element_blank())+
          theme(axis.ticks = element_blank())+
          theme(panel.border = element_blank())+
          theme(axis.title = element_blank())+
          theme(panel.background = element_blank(),
                panel.border = element_blank(), 
                legend.background = element_blank(),
                legend.position = "right",
                legend.text = element_text(size=16),
                legend.title=element_text(size = 15, face ="bold"),		
                legend.key.size = unit(0.25, "inches"))	
        df = vln.dat  		
        p3 <- ggplot(df, aes(x, y,color=Score))+
          geom_point(size = piontsize, alpha = 0.9, shape = 15)+
          scale_color_gradientn('Score Value', colours = mycolor)+
          scale_x_continuous(limits = c(xmin-2, xmax))+
          scale_y_continuous(limits = c(ymin-2, ymax))+
          theme(axis.text = element_blank())+
          theme(axis.ticks = element_blank())+
          theme(panel.border = element_blank())+
          theme(axis.title = element_blank())+
          theme(panel.background = element_blank(),
                panel.border = element_blank(), 
                legend.background = element_blank(),
                legend.position = "right",
                legend.text = element_text(size=16),
                legend.title=element_text(size = 15, face ="bold"),		
                legend.key.size = unit(0.25, "inches"))
        
        df_high_score_cells <- df %>%
          filter(Group != "NonTopCell")			
        p4 <- ggplot(df, aes(x = x, y = y)) + 
          geom_hdr(data = df_high_score_cells, aes(fill = Group), color = "white", size = 1, show.legend = FALSE) +				
          scale_fill_manual(values = "orange") + 
          geom_point(aes(color = Group, size = Group), alpha = 0.25, shape = 20) +	
          scale_color_manual(values = c("TopCell" = "red", "NonTopCell" = "grey85"),
                             labels = c("TopCell" = "Yes", "NonTopCell" = "No")) +  
          scale_size_manual(values = c("TopCell" = 1.5*piontsize, "NonTopCell" = piontsize),
                            labels = c("TopCell" = "Yes", "NonTopCell" = "No")) + 
          labs(color = "TopCell", size = "TopCell") +  
          guides(
            color = guide_legend(override.aes = list(size = 2.5)), 
            size = guide_legend(override.aes = list(size = 2.5))   
          ) +				
          scale_x_continuous(limits = c(xmin-2, xmax))+
          scale_y_continuous(limits = c(ymin-2, ymax))+
          theme(axis.text = element_blank())+
          theme(axis.ticks = element_blank())+
          theme(panel.border = element_blank())+
          theme(axis.title = element_blank())+
          theme(panel.background = element_blank(),
                panel.border = element_blank(), 
                legend.background = element_blank(),
                legend.position = "right",
                legend.text = element_text(size=16),
                legend.title=element_text(size = 15, face ="bold"),		
                legend.key.size = unit(0.25, "inches"))
        ggsave(paste0(absolute_path, "/", current_column, ".GeneModuleDistribution.png"), (p1|p2)/(p3|p4), width = 14, height = 10, units = "in")
        
        ix = is.na(dataDD$window) 
        dataDD = dataDD[!ix, ]
        dataDD <- dataDD %>%
          group_by(window) %>%
          mutate(			
            P_c = ifelse(Score >= 0, exp(-((Score - mu_top)^2) / (2 * sigma_top^2)), 0),
            top_scores = ifelse(Score >= top_percent, Score, NA),
            N_w = sum(P_c),
            S_w = sum(P_c * Score)
          ) %>%
          ungroup() 
        
        igsi_values <- dataDD %>%
          dplyr::select(window, N_w, S_w) %>%
          unique() %>%
          mutate(
            total_sum = sum(N_w * S_w),
            IGSI = (N_w * S_w) / total_sum
          )
        dataDD <- left_join(dataDD, igsi_values %>% dplyr::select(window, IGSI), by = "window")	
        
        dataDD$xnew <- as.numeric(dataDD$xnew)
        dataDD = dataDD[order(dataDD$xnew), ]
        
        dataDD1 <- unique(dataDD[,c('xnew', 'IGSI')])		
        
        segments <- segment(
          smooth.CNA(
            CNA(
              dataDD1$IGSI,  
              as.factor('Region'), 
              dataDD1$xnew,
              data.type = 'logratio',
              sampleid = 'Region'
            ),
            smooth.region = 20 
          ),
          undo.splits = 'sdundo',
          undo.SD = SD,
          min.width = 5,       
          alpha = alpha
        )$output[, 2:6]
        
        segments <- segments %>%
          mutate(segment = row_number())
        
        quantile_value <- quantile(segments$seg.mean, 1 - topSegFreq, na.rm = TRUE)
        
        segments$HGMSR1 <- ifelse(segments$seg.mean >= quantile_value, 'HGMSR', 'Non-HGMSR')
        
        dataDD <- dataDD %>%
          interval_join(
            segments,
            by = c("xnew" = "loc.start", "xnew" = "loc.end"),
            mode = "left" 
          ) %>%
          
          dplyr::select(cell, window, xnew, ynew, Score, IGSI, segment, HGMSR1, top_scores, x, y, seg.mean) %>%
          arrange(window)
        
        dataDD <- dataDD %>%
          group_by(cell, window, xnew, ynew, Score, IGSI, x, y) %>% 
          mutate(diff = abs(Score - seg.mean)) %>%  
          filter(diff == min(diff, na.rm = TRUE)) %>% 
          dplyr::select(-diff, -seg.mean) %>%  
          ungroup() %>% 
          arrange(window)  
        
        dataDD <- dataDD %>%
          group_by(segment) %>%
          summarise(
            total_non_na = sum(!is.na(top_scores)),
            .groups = "drop"
          ) %>%
          mutate(total_non_na_all = sum(total_non_na)) %>%
          mutate(proportion = total_non_na / total_non_na_all) %>%
          left_join(dataDD, by = "segment") %>%
          mutate(
            SegType = ifelse(HGMSR1 == "HGMSR" & proportion > TopCellFreq, "HGMSR", "Non-HGMSR")
          ) %>%	  
          dplyr::select(-HGMSR1) %>%
          ungroup()
        
        region_means <- dataDD %>%
          group_by(segment) %>%
          summarise(
            mean_value = mean(IGSI, na.rm = TRUE),
            xnew_start = min(xnew),
            xnew_end = max(xnew),
            SegType = ifelse(any(SegType == "HGMSR", na.rm = TRUE), "HGMSR", "Non-HGMSR")  
          )
        
        dataDD <- dataDD %>%
          left_join(region_means %>% dplyr::select(segment, mean_value), by = "segment") %>%
          mutate(IGSI = coalesce(mean_value, IGSI)) %>%
          dplyr::select(-mean_value) %>%
          arrange(window)
        
        min_Score <- 0 
        max_Score <- max(dataDD$Score)
        min_mean_value <- 0 
        max_mean_value <- max(region_means$mean_value)
        
        scale_ratio <- (max_Score - min_Score) / (max_mean_value - min_mean_value)
        region_means$mean_value_scaled <- (region_means$mean_value - min_mean_value) * scale_ratio + min_Score
        sec_axis_inv_trans <- function(x) {
          return((x - min_Score) / scale_ratio + min_mean_value)
        }
        
        x_range <- range(dataDD$xnew)
        x_breaks <- seq(floor(x_range[1]/1000)*1000, ceiling(x_range[2]/1000)*1000, by = 1000)
        if (length(x_breaks) > 10) {
          interval <- 1000 * ceiling((max(x_breaks) - min(x_breaks)) / 10000)
          x_breaks <- seq(floor(x_range[1]/interval)*interval, ceiling(x_range[2]/interval)*interval, by = interval)
        }
        
        
        nameDD = paste0('Gene ', paste0(strsplit(current_column,'ScoreGeneSet')[[1]],collapse = " "))
        p <- ggplot(dataDD, aes(x = xnew, y = Score)) +
          geom_point(color = "chartreuse3", size = 2) +
          geom_segment(data = region_means, 
                       aes(x = xnew_start, xend = xnew_end, y = mean_value_scaled, yend = mean_value_scaled, 
                           color = SegType, size = SegType),  
                       show.legend = TRUE) + 
          scale_color_manual(values = c("HGMSR" = "magenta", "Non-HGMSR" = "#D3D3D3"), 
                             labels = c("HGMSR" = "HGMSR", "Non-HGMSR" = "Non-HGMSR")) + 
          scale_size_manual(values = c("HGMSR" = 2, "Non-HGMSR" = 0.8),  
                            guide = "none") +  
          labs(title = nameDD, 
               x = "Position", 
               y = "Integrated Genomic Score Index",
               color = "SegType") +  
          theme_void() +
          theme(
            axis.ticks = element_line(color = "black", size = 1), 
            axis.ticks.length = unit(2.5, "mm"),
            axis.text = element_text(color = "black", margin = margin(2, 2, 2, 2, "mm")),
            axis.title.x = element_text(color = "black"),
            axis.title.y = element_text(color = "black", angle = 90, vjust = 0.5),
            axis.line = element_line(color = "black", size = 1),
            legend.position = "right", 
            legend.text = element_text(color = "black"), 
            legend.title = element_text(color = "black") 
          ) +	
          scale_y_continuous(
            name = "Score Value",
            sec.axis = sec_axis(~ sec_axis_inv_trans(.), name = "Integrated Genomic Score Index")  
          )+
          scale_x_continuous(
            breaks = x_breaks,  
            minor_breaks = x_breaks  
          ) 
        ggsave(paste0(absolute_path, "/", current_column, ".png"), p, width = 20, height = 5, units = "in")
        
        
        
        region_means <- region_means %>%
          ungroup() %>%
          dplyr::select(-mean_value_scaled)		
        write.csv(region_means, paste0(absolute_path, "/", current_column, ".IGSI_Segment.csv"), row.names = FALSE)  
        
        if(is.null(arrayMean))
        {
          arrayMean = dataDD[, c('cell','SegType','xnew','IGSI')]
          colnames(arrayMean)[4] = current_column
        }else{
          arrayMean_tmp = dataDD[, c('cell', 'IGSI')]
          colnames(arrayMean_tmp)[2] = current_column
          arrayMean = merge(arrayMean, arrayMean_tmp, by = 'cell', all = F)
        }
        
        ix = dataDD$SegType == 'HGMSR'
        nregion = unique(dataDD$segment[ix])
        
        if (!"CandidateArea" %in% colnames(dataDD)) {
          dataDD$CandidateArea <- NA_character_
        }
        if (!"MCS" %in% colnames(dataDD)) {
          dataDD$MCS <- NA_character_
        }		
        if (!"MCSR" %in% colnames(dataDD)) {
          dataDD$MCSR <- NA_character_
        }	
        
        if(length(nregion)>0)
        {
          setlist <- process_regions(nregion, range_min = min(dataDD$segment), range_max = max(dataDD$segment))
          jnum <- 1
          for(j in 1:length(setlist))
          {	
            i_current_column = paste0(current_column,'.CandidateArea.',jnum)
            ix = dataDD$segment %in% setlist[[j]] 
            data = dataDD[ix, ]
            data <- data %>%
              mutate(Group = ifelse(is.na(top_scores), "NonTopCell", "TopCell"))
            data$SegType = factor(data$SegType, levels = c("Non-HGMSR", "HGMSR"))
            data <- data %>%
              mutate(scale_score = (Score - min(Score)) / (max(Score) - min(Score)))			
            
            high_score_cells <- data[data$Group == 'TopCell', ] 		
            if(nrow(high_score_cells)>0)
            {
              dbscan_res <- dbscan(high_score_cells[, c("x", "y")], 
                                   eps = eps, 
                                   minPts = 4)				 
              high_score_cells$cluster <- dbscan_res$cluster									
              clusterN_cells1 <- high_score_cells %>%
                filter(cluster != 0) %>%
                dplyr::select(cell) %>%
                pull(cell)  
              sf_data <- st_as_sf(data, coords = c("x", "y"))
              coords <- st_coordinates(sf_data)
              nb <- try(knn2nb(knearneigh(coords, k = 4)),silent = TRUE)
              if(inherits(nb, "try-error"))
              {
                clusterN_cells2 <- NULL
              }else{
                listw <- nb2listw(nb, style="B") 
                gi_star <- localG(data$Score, listw = listw)     
                data$hotspot <- ifelse(gi_star > 1.96, "Hotspot", "Non-significant")
                clusterN_cells2 <- data %>%
                  filter(hotspot == 'Hotspot') %>%
                  dplyr::select(cell) %>%
                  pull(cell) 
              }
              
              clusterN_cells = intersect(clusterN_cells1, clusterN_cells2)					  
              if(length(clusterN_cells)>0)
              {
                points <- data %>%
                  filter(cell %in% clusterN_cells) %>%
                  dplyr::select(cell, x, y) %>%
                  as.data.frame()
                correlated_cells <- find_nearest_and_correlated_cells(points, data, expGeneList, k = kneighbor, threshold = corr)
                correlated_cells_vector <- unique(as.character(unlist(correlated_cells))) 
                
                data <- data %>%
                  mutate(
                    MCS = case_when(
                      cell %in% correlated_cells_vector ~ "YES",
                      TRUE ~ "NO"
                    )
                  )
                
                data$MCS = factor(data$MCS, levels = c('YES','NO'))
                
                numMixMCS = 5	
                ix = data$MCS == 'YES'
                if(sum(ix)>=numMixMCS)
                {
                  
                  mcs_high_score_cells = data[ix, ]
                  
                  hdr_data <- ggdensity::stat_hdr(
                    data = mcs_high_score_cells,
                    aes(x = x, y = y, fill = MCS),
                    probs = c(0.95, 0.50)  
                  )
                  
                  contour_points <- layer_data(ggplot() + hdr_data)
                  contour_polygons <- contour_points %>%
                    group_by(probs, piece, group, subgroup) %>%
                    summarise(geometry = list({
                      coords <- cbind(x, y)
                      if(!all(coords[1,] == coords[nrow(coords),])) {
                        coords <- rbind(coords, coords[1,])
                      }
                      st_polygon(list(coords))
                    })) %>%
                    st_as_sf()
                  data_sf <- st_as_sf(data, coords = c("x", "y"))
                  contour_polygons <- contour_polygons %>%
                    mutate(priority = case_when(
                      probs == "50%" ~ 1,
                      probs == "95%" ~ 2,
                      TRUE ~ 3 
                    )) %>%
                    arrange(priority) %>%  
                    mutate(poly_id = row_number())
                  
                  zones <- st_join(
                    data_sf %>% mutate(point_id = row_number()),  
                    contour_polygons,
                    join = st_within,
                    largest = FALSE,
                    left = TRUE
                  ) %>%
                    group_by(point_id) %>%  
                    arrange(priority) %>%    
                    slice(1) %>%           
                    ungroup()
                  
                  zones <- zones %>%
                    mutate(zone = case_when(
                      probs == "50%" ~ "CFD",
                      probs == "95%" ~ "DAF",
                      is.na(probs) ~ "PBR",
                      TRUE ~ "PBR"  
                    )) %>% 
                    dplyr::select(-point_id, -poly_id)
                  
                  data <- data %>%
                    left_join(
                      st_drop_geometry(zones) %>% 
                        dplyr::select(cell, zone),  
                      by = "cell"  
                    ) %>%
                    mutate(MCSR = zone) %>%  
                    dplyr::select(-zone)  
                  data$MCSR = factor(data$MCSR, levels = c('CFD','DAF','PBR'))
                  #numMixMCS = 5
                  if(sum(table(data$MCSR) < numMixMCS) == 0) 						
                  {
                    ixx = data$MCS == 'YES'						
                    pvalue <- 1
                    ptest <- try(wilcox.test(data$Score[ixx], data$Score[!ixx], alternative = "greater"),silent=TRUE)
                    if(!inherits(ptest, "try-error"))
                    {
                      pvalue = ptest$p.value
                    }
                    if(pvalue < 0.05)
                    {					
                      signif_label <- ifelse(pvalue < 0.001, "***", 
                                             ifelse(pvalue < 0.01, "**", 
                                                    ifelse(pvalue < 0.05, "*", "ns")))
                      
                      p <- ggplot(data = data, aes(x = MCS, y = Score, fill = MCS)) +
                        geom_half_violin(color = NA, side = "r", trim = FALSE) +
                        geom_boxplot(aes(color = MCS), fill = "white", outlier.shape = NA, alpha = 0.8, width = 0.15) +
                        geom_half_point(side = "l", shape = 16, size = 0.5, color = 'grey50') +
                        geom_text(aes(x = 1.5, y = max(data$Score) * 1.1, label = signif_label), 
                                  size = 5, fontface = "bold", color = "black") +
                        scale_fill_manual(values = c("YES" = "darkred", "NO" = "darkgreen"), guide = "legend", drop = FALSE) +
                        scale_color_manual(values = c("YES" = "darkred", "NO" = "darkgreen"), guide = "legend", drop = FALSE) +
                        labs(x = 'Group', y = 'Score Value') +
                        theme_bw() +
                        theme(panel.grid = element_blank()) +
                        theme(axis.text.x = element_text(family = NULL, size = 8, angle = 45, hjust = 1, vjust = 1))
                      
                      ggsave(paste0(absolute_path, "/", i_current_column, '.MSC-boxPlotOneSidedTest.png'), p, width = 8, height = 6)
                      
                      gnames = levels(data$MCSR)
                      gnum = nlevels(data$MCSR)
                      groupDDDINDEX = combn(gnum, 2)
                      comparisonslist = list()
                      for(i in 1:ncol(groupDDDINDEX))
                      {
                        comparisonslist[[i]] = gnames[groupDDDINDEX[,i]]
                      }
                      
                      
                      p <- ggplot(data = data, aes(x = MCSR, y = Score, fill = MCSR)) +
                        geom_half_violin(color = NA, side = "r", trim = FALSE) +
                        geom_boxplot(aes(color = MCSR), fill = "white", outlier.shape = NA, alpha = 0.8, width = 0.15) +
                        geom_half_point(side = "l", shape = 16, size = 0.5, color = 'grey50') +
                        geom_signif(comparisons=comparisonslist,
                                    margin_top = 0.1,
                                    step_increase=0.08,
                                    test="wilcox.test",
                                    textsize=2.8,
                                    map_signif_level = TRUE) +  
                        scale_fill_manual(values=mycolor, guide="none") +
                        scale_color_manual(values=mycolor, guide="none") +	
                        scale_fill_manual(values = c("CFD" = "darkred", 'DAF' = "darkblue", "PBR" = "darkgreen"), guide = "legend", drop = FALSE) +
                        scale_color_manual(values = c("CFD" = "darkred", 'DAF' = "darkblue","PBR" = "darkgreen"), guide = "legend", drop = FALSE) +										  
                        labs(x = 'Group', y = 'Score Value') +
                        theme_bw() +
                        theme(panel.grid = element_blank()) +
                        theme(axis.text.x = element_text(family = NULL, size = 8, angle = 45, hjust = 1, vjust = 1))  
                      ggsave(paste0(absolute_path, "/", i_current_column, '.MSC-boxPlotTwoSidedTest.png'), p, width = 8, height = 6)
                      
                      center_cells <- data %>%
                        filter(MCS == "YES") %>%
                        dplyr::select(cell) %>%
                        pull(cell)
                      
                      dataDDplot <- dataDD %>%
                        mutate(MCS = ifelse(cell %in% center_cells, "YES", "NO"))
                      
                      p <- ggplot(dataDDplot, aes(x = x, y = y)) + 
                        geom_hdr(data = mcs_high_score_cells, aes(fill = MCS), color = "white", size = 0.3, show.legend = FALSE) +				
                        scale_fill_manual(values = "orange") +  
                        geom_point(aes(color = MCS, size = MCS), alpha = 0.25, shape = 20) +	
                        scale_color_manual(values = c("YES" = "red", "NO" = "grey65"),
                                           labels = c("YES" = "Yes", "NO" = "No")) +  
                        scale_size_manual(values = c("YES" = piontsize, "NO" = piontsize),
                                          labels = c("YES" = "Yes", "NO" = "No")) +  
                        labs(color = "MCS", size = "MCS") +  
                        guides(
                          color = guide_legend(override.aes = list(size = 2.5)), 
                          size = guide_legend(override.aes = list(size = 2.5))   
                        ) +				
                        scale_x_continuous(limits = c(xmin-2, xmax))+
                        scale_y_continuous(limits = c(ymin-2, ymax))+
                        theme(axis.text = element_blank())+
                        theme(axis.ticks = element_blank())+
                        theme(panel.border = element_blank())+
                        theme(axis.title = element_blank())+
                        theme(panel.background = element_blank(),
                              panel.border = element_blank(), 
                              legend.background = element_blank(),
                              legend.position = "right",
                              legend.text = element_text(size=16),
                              legend.title=element_text(size = 15, face ="bold"),		
                              legend.key.size = unit(0.25, "inches"))
                      ggsave(paste0(absolute_path, "/", i_current_column,'.MCS.pdf'), p, width = 10, height = 10)
                      ggsave(paste0(absolute_path, "/", i_current_column,'.MCS.png'), p, width = 10, height = 10)
                      
                      p1 <- ggplot(data, aes(x = x, y = y)) +  
                        geom_point(aes(color = MCS, size = 0.01 * scale_score/5), 
                                   shape = 20, 
                                   alpha = 0.2,
                                   show.legend = c(color = TRUE, size = FALSE)) +  
                        scale_size_continuous(range = c(0, 1), guide = "none") +  
                        scale_color_manual(
                          values = c("YES" = "red", "NO" = "grey85"),
                          name = "MCS",  
                          guide = guide_legend(
                            override.aes = list(
                              size = 6  
                            )
                          )	
                        ) +
                        scale_x_continuous(breaks = seq(min(data$x), max(data$x), length.out = 3)) +
                        theme(
                          plot.background = element_blank(),
                          legend.position = "right",
                          legend.key.size = unit(0.5, "cm"), 
                          legend.spacing = unit(0.2, "cm"),   
                          legend.box.spacing = unit(0.2, "cm"), 
                          legend.margin = margin(0, 0, 0, 0, "cm"), 
                          axis.line = element_line(color = "black", size = 1.2),
                          panel.border = element_blank(),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          panel.background = element_blank(),
                          axis.ticks = element_line(color = "black", size = 1.2),
                          axis.ticks.length = unit(1, "mm")
                        ) +
                        labs(x = 'x', y = 'y') +            
                        ggtitle("Module Co-Activity Spots Region")
                      
                      p1 <- p1 +   
                        geom_sf(
                          data = contour_polygons, 
                          aes(linetype = probs),      
                          fill = NA, 
                          color = "black",
                          inherit.aes = FALSE
                        ) +
                        scale_linetype_manual(
                          values = c("50%" = "solid", "95%" = "longdash"),
                          labels = c("50%" = "CFD", "95%" = "DAF"),
                          breaks = c("50%", "95%"),  
                          name = "Zone",
                          guide = guide_legend(
                            override.aes = list(
                              size = 0.8,
                              linetype = c("solid", "longdash")
                            )
                          )
                        )
                      ggsave(paste0(absolute_path, "/", i_current_column, '.MCS_Local.pdf'), p1, width = 8, height = 6)
                      
                      kriging_res <- kriging_interpolation(data, grid_density)
                      r <- rasterFromXYZ(kriging_res)
                      gradient_result <- calculate_gradient(r)
                      grad_x <- gradient_result$grad_x
                      grad_y <- gradient_result$grad_y
                      vector_df <- data.frame(
                        x = coordinates(r)[,1],
                        y = coordinates(r)[,2],
                        u = values(grad_x),
                        v = values(grad_y)
                      )
                      
                      p2 <- ggplot() +
                        geom_raster(data = raster::as.data.frame(r, xy = TRUE), 
                                    aes(x = x, y = y, fill = layer)) +	 
                        geom_quiver(data = vector_df, 
                                    aes(x = x, y = y, u = u, v = v), 
                                    color = "white", arrow.length = 0.05, size = 0.5) +			 
                        scale_fill_viridis_c(option = "inferno", name = "Energy Field Value") +
                        theme(
                          plot.background = element_blank(),  
                          axis.line = element_line(color = "black", size = 1.2), 
                          panel.border = element_blank(),   
                          panel.grid.major = element_blank(), 
                          panel.grid.minor = element_blank(),  
                          panel.background = element_blank(),  
                          axis.ticks = element_line(color = "black", size = 1.2),  
                          axis.ticks.length = unit(2, "mm")  
                        ) +    
                        labs(title = "Energy Field Gradient")
                      ggsave(paste0(absolute_path, "/", i_current_column,'.Getis-Ord_Gi.pdf'), p2, width = 6, height = 6)	
                      
                      vector_df$magnitude <- sqrt(vector_df$u^2 + vector_df$v^2)
                      p3 <- ggplot(vector_df, aes(x, y, u = u, v = v, color = magnitude)) +
                        geom_quiver(arrow.length = 0.05) +
                        scale_color_viridis_c(option = "plasma", name = "Gradient Magnitude") +
                        theme(
                          plot.background = element_blank(),  
                          axis.line = element_line(color = "black", size = 1.2),  
                          panel.border = element_blank(),   
                          panel.grid.major = element_blank(),  
                          panel.grid.minor = element_blank(),  
                          panel.background = element_blank(),  
                          axis.ticks = element_line(color = "black", size = 1.2),  
                          axis.ticks.length = unit(2, "mm")  
                        ) +			  
                        ggtitle("Gradient Magnitude Distribution")  
                      
                      combined_plot <- wrap_plots(p1, p2, p3, nrow = 1, widths = c(1, 1, 1), heights = c(1))
                      ggsave(paste0(absolute_path, "/", i_current_column,'.MCS-TrendDirection.pdf'), combined_plot, width = 20, height = 5)
                      ggsave(paste0(absolute_path, "/", i_current_column,'.MCS-TrendDirection.png'), combined_plot, width = 20, height = 5)
                      
                      
                      dataDD <- dataDD %>%
                        mutate(MCS = case_when(
                          is.na(MCS) & cell %in% data$cell ~ data$MCS[match(cell, data$cell)],  
                          TRUE ~ MCS  
                        )) %>%							
                        mutate(MCSR = case_when(
                          is.na(MCSR) & cell %in% data$cell ~ data$MCSR[match(cell, data$cell)], 
                          TRUE ~ MCSR  
                        )) %>%
                        mutate(CandidateArea = case_when(
                          is.na(CandidateArea) & cell %in% data$cell ~ as.character(jnum), 
                          TRUE ~ CandidateArea  
                        ))
                      jnum = jnum + 1							  														  
                    }
                  }
                }
              }
            }
          }
        }	 
        dataDD <- dataDD %>%
          ungroup() %>%
          dplyr::select(-top_scores)		
        write.csv(dataDD, paste0(absolute_path, "/", current_column, ".IGSI_Cell.csv"), row.names = FALSE)		
        
      }
    }
  }
  
  if(!is.null(arrayMean))
  {
    arrayMean = arrayMean[order(arrayMean$xnew),]
    arrayMean_long <- reshape2::melt(arrayMean, id.vars = c("cell", "SegType", "xnew"), 
                                     variable.name = "ModuleScore", value.name = "Score")
    arrayMean_long$ModuleScore = paste0('Gene ', gsub('ScoreGeneSet', ' ', arrayMean_long$ModuleScore))					   
    p <- ggplot(arrayMean_long, aes(x = xnew, y = Score, color = ModuleScore)) +
      geom_line(size = 1.5) + 
      labs(title = "Module Score Curves", x = "position", y = "IGSI") + 
      theme_minimal() +
      theme(
        legend.position = "none",  
        plot.title = element_text(hjust = 0.5)  
      ) +
      facet_wrap(~ ModuleScore, ncol = 1) 
    ggsave(paste0(absolute_path, "/module_score_curves.pdf"), plot = p, width = 12, height = 6, units = "in", dpi = 300)
  }
  
}