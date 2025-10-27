autoReplotHE <- function(df, clusters, color, labelfontsize, axistitlefontsize, 
                         spatialpointsize, pointalpha, polygonalpha, pointshape, legendtextsize, 
                         legendtitlesize, guidelegendsize, ratio, freq, index = 1,
                         pards = 0, point_x, point_y, interval = 0.05, intervalAngle = 30, Signal = FALSE)
{ 
  max_x <- max(df$x)
  min_x <- min(df$x)
  max_y <- max(df$y)
  min_y <- min(df$y)
  margin_x <- (max_x - min_x) * 0.1
  margin_y <- (max_y - min_y) * 0.1
  xlim_range <- c(min_x - margin_x, max_x + margin_x)
  ylim_range <- c(min_y - margin_y, max_y + margin_y)
  interval1 <- (max_x - min_x) * interval
  radii <- seq(0, max(max_x - point_x, point_x - min_x) + interval1, by = interval1)
  circle_data <- data.frame()
  for (radius in radii) {
    theta <- seq(0, 2 * pi, length.out = 100) 
    x <- point_x + radius * cos(theta)
    y <- point_y + radius * sin(theta)
    group <- rep(radius, length(theta)) 
    circle_data <- rbind(circle_data, data.frame(x = x, y = y, group = group))
  }  
  
  convex_hulls <- df %>%
    group_by(dataset) %>%
    do({
      points_matrix <- as.matrix(.[, c("x", "y")])
      hull_indices <- chull(points_matrix)
      hull_points <- points_matrix[hull_indices, ]
      data.frame(x = hull_points[, 1], y = hull_points[, 2], dataset = unique(.$dataset))
    })
  
  ix = df$dataset == clusters[1]
  df = df[ix, ]
  ix = convex_hulls$dataset == clusters[2]
  contour_df = convex_hulls[ix, ]  
  
  ix = df$dataset == clusters[1]
  df = df[ix, ]
  ix = contour_df$dataset == clusters[2]
  contour_df = contour_df[ix, ]  
  
  
  if(Signal)
  {
    p <- ggplot() +
      geom_point(data = df, aes(x = x, y = y, color = dataset), size = spatialpointsize, alpha = pointalpha) + 
      geom_polygon(data = contour_df, aes(x = x, y = y, fill = dataset), alpha = polygonalpha) +  
      scale_color_manual(values = color, labels = c("Signal" = "Signal", "SignalPanel" = "SignalPanel")) +  
      scale_fill_manual(values = color, labels = c("Signal" = "Signal", "SignalPanel" = "SignalPanel")) +  
      guides(color = guide_legend(title = ""), fill = guide_legend(title = "")) +
      theme(panel.background = element_blank(),
            panel.border = element_blank(), 
            legend.background = element_blank(),
            legend.box = "horizontal",  # Changed from element_blank() to valid value			
            legend.text = element_text(size = legendtextsize),
            legend.title = element_text(size = legendtitlesize),				
            legend.position = "right",
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            plot.title = element_text(hjust = 1)) +  
      ggtitle(paste0('Index: ',index,"\n",'Pixel Area Ratio Difference Sum: ',pards,"\n",'Registration Coefficient: ',ratio,"\n",'spot proportion: ',freq, '%',"\nConcentric Circle Radius Interval: ",interval,"\nDeflection Angle Interval: ",intervalAngle)) +  
      guides(color=guide_legend(title = '', override.aes = list(size = guidelegendsize), ncol = 2)) +
      geom_vline(xintercept = point_x, color = "yellow", size = 0.2, alpha = 0.8) + 
      geom_hline(yintercept = point_y, color = "yellow", size = 0.2, alpha = 0.8) +
      geom_text(aes(x = point_x, y = point_y + radii[length(radii)-1], label = "up"), vjust = -2, color = "black") +  
      geom_text(aes(x = point_x, y = point_y - radii[length(radii)-1], label = "down"), vjust = 2.5, color = "black") + 
      geom_text(aes(x = point_x - radii[length(radii)-1], y = point_y, label = "left"), hjust = 2, color = "black") +  
      geom_text(aes(x = point_x + radii[length(radii)-1], y = point_y, label = "right"), hjust = -0.5, color = "black") + 
      scale_x_continuous(limits = xlim_range) +  
      scale_y_continuous(limits = ylim_range) +
      geom_path(data = circle_data, aes(x = x, y = y, group = group), 
                color = "yellow", linetype = "dashed", size = 0.8, inherit.aes = FALSE)
  }else{
    p <- ggplot() +
      geom_point(data = df, aes(x = x, y = y, color = dataset), size = spatialpointsize, alpha = pointalpha) +  
      geom_polygon(data = contour_df, aes(x = x, y = y, fill = dataset), alpha = polygonalpha) +  
      scale_color_manual(values = color, labels = c("HE" = "HE", "ssDNA" = "ssDNA")) + 
      scale_fill_manual(values = color, labels = c("HE" = "HE", "ssDNA" = "ssDNA")) +  
      guides(color = guide_legend(title = ""), fill = guide_legend(title = "")) +
      theme(panel.background = element_blank(),
            panel.border = element_blank(), 
            legend.background = element_blank(),
            legend.box = "horizontal",  # Changed from element_blank() to valid value		
            legend.text = element_text(size = legendtextsize),
            legend.title = element_text(size = legendtitlesize),				
            legend.position = "right",
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            plot.title = element_text(hjust = 1)) +  
      ggtitle(paste0('Index: ',index,"\n",'Pixel Area Ratio Difference Sum: ',pards,"\n",'Registration Coefficient: ',ratio,"\n",'spot proportion: ',freq, '%',"\nConcentric Circle Radius Interval: ",interval,"\nDeflection Angle Interval: ",intervalAngle)) +  
      guides(color=guide_legend(title = 'Clusters', override.aes = list(size = guidelegendsize), ncol = 2)) +
      geom_vline(xintercept = point_x, color = "yellow", size = 0.2, alpha = 0.8) +  
      geom_hline(yintercept = point_y, color = "yellow", size = 0.2, alpha = 0.8) +
      geom_text(aes(x = point_x, y = point_y + radii[length(radii)-1], label = "up"), vjust = -2, color = "black") +  
      geom_text(aes(x = point_x, y = point_y - radii[length(radii)-1], label = "down"), vjust = 2.5, color = "black") + 
      geom_text(aes(x = point_x - radii[length(radii)-1], y = point_y, label = "left"), hjust = 2, color = "black") +  
      geom_text(aes(x = point_x + radii[length(radii)-1], y = point_y, label = "right"), hjust = -0.5, color = "black") + 
      scale_x_continuous(limits = xlim_range) +  
      scale_y_continuous(limits = ylim_range) +
      geom_path(data = circle_data, aes(x = x, y = y, group = group), 
                color = "yellow", linetype = "dashed", size = 0.8, inherit.aes = FALSE)
  }  
  
  angles <- seq(0, 180, by = intervalAngle) 
  slopes <- tan(angles * pi / 180)
  intercepts <- point_y - slopes * point_x
  line_data <- data.frame(slope = slopes, intercept = intercepts)
  p <- p + geom_abline(data = line_data, aes(slope = slope, intercept = intercept), color = "yellow", linetype = "dashed", size = 0.2, alpha = 0.8)
  
  return(p)
}


