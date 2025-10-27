autoReplot_old <- function(df, clusters, color, labelfontsize, axistitlefontsize, 
                           spatialpointsize, pointalpha, pointshape, legendtextsize, 
                           legendtitlesize, guidelegendsize, ratio, freq, index = 1,
                           pards = 0, point_x, point_y, interval = 0.05, intervalAngle = 30)
{
  suppressMessages(library(ggforce))  
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
  
  p <- ggplot() +
    geom_point(df, mapping=aes(x = x, y = y, color = eval(str2lang(clusters))), size = spatialpointsize, alpha = pointalpha, shape = pointshape) + 
    scale_color_manual(values = color) +
    theme(panel.background = element_blank(),
          panel.border = element_blank(), 
          legend.background = element_blank(),
          legend.box = element_blank(),			
          legend.text = element_text(size = legendtextsize),
          legend.title = element_text(size = legendtitlesize),				
          legend.position = "right",
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          plot.title = element_text(hjust = 1)) +  
    ggtitle(paste0('Index: ',index,"\n",'Pixel Area Ratio Difference Sum: ',pards,"\n",'Registration Coefficient: ',ratio,"\n",'spot proportion: ',freq, '%',"\nConcentric Circle Radius Interval: ",interval)) +  
    guides(color=guide_legend(title = 'Clusters', override.aes = list(size = guidelegendsize), ncol = 2)) +
    geom_vline(xintercept = point_x, color = "yellow") +  
    geom_hline(yintercept = point_y, color = "yellow") +
    geom_text(aes(x = point_x, y = point_y + radii[length(radii)-1], label = "up"), vjust = -2, color = "black") +  
    geom_text(aes(x = point_x, y = point_y - radii[length(radii)-1], label = "down"), vjust = 2.5, color = "black") +  
    geom_text(aes(x = point_x - radii[length(radii)-1], y = point_y, label = "left"), hjust = 2, color = "black") +  
    geom_text(aes(x = point_x + radii[length(radii)-1], y = point_y, label = "right"), hjust = -0.5, color = "black") +
    scale_x_continuous(limits = xlim_range) + 
    scale_y_continuous(limits = ylim_range) 
  
  p <- p + geom_circle(aes(x0 = point_x, y0 = point_y, r = radii), color = "yellow", linetype = "dashed", size = 0.8, inherit.aes = FALSE)
  
  return(p)
}


autoReplot_V0 <- function(df, clusters, color, labelfontsize, axistitlefontsize, 
                          spatialpointsize, pointalpha, pointshape, legendtextsize, 
                          legendtitlesize, guidelegendsize, ratio, freq, index = 1,
                          pards = 0, point_x, point_y, interval = 0.05, intervalAngle = 30)
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
  
  p <- ggplot() +
    geom_point(df, mapping=aes(x = x, y = y, color = eval(str2lang(clusters))), size = spatialpointsize, alpha = pointalpha, shape = pointshape) + 
    scale_color_manual(values = color) +
    theme(panel.background = element_blank(),
          panel.border = element_blank(), 
          legend.background = element_blank(),
          legend.box = element_blank(),			
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
  
  angles <- seq(0, 180, by = intervalAngle) 
  slopes <- tan(angles * pi / 180)
  intercepts <- point_y - slopes * point_x
  line_data <- data.frame(slope = slopes, intercept = intercepts)
  
  p <- p + geom_abline(data = line_data, aes(slope = slope, intercept = intercept), color = "yellow", linetype = "dashed", size = 0.2, alpha = 0.8)
  
  return(p)
}

autoReplot <- function(df, clusters, color, labelfontsize, axistitlefontsize, 
                       spatialpointsize, pointalpha, pointshape, legendtextsize, 
                       legendtitlesize, guidelegendsize, ratio, freq, index = 1,
                       pards = 0, point_x, point_y, interval = 0.05, intervalAngle = 30)
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
  
  p <- ggplot() +
    geom_point(df, mapping=aes(x = x, y = y, color = eval(str2lang(clusters))), size = spatialpointsize, alpha = pointalpha, shape = pointshape) + 
    scale_color_manual(values = color) +
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
  
  angles <- seq(0, 180, by = intervalAngle) 
  slopes <- tan(angles * pi / 180)
  intercepts <- point_y - slopes * point_x
  line_data <- data.frame(slope = slopes, intercept = intercepts)
  
  p <- p + geom_abline(data = line_data, aes(slope = slope, intercept = intercept), color = "yellow", linetype = "dashed", size = 0.2, alpha = 0.8)
  
  return(p)
}
