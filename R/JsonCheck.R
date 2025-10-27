JsonCheck = function(json)
{
	suppressMessages(library(jsonlite))
	jsonData = fromJSON(json)  
	absolute_path <- dirname(normalizePath(json))
	png(file = paste0(absolute_path, "/JsonPolygonCheck.png"), width = 500, height = 500)
	i = which(jsonData$shapes$label == 'Wai')
	j = which(jsonData$shapes$label == 'Nei')
	jsonData$shapes$points[[i]]
	n = nrow(jsonData$shapes$points[[j]])
	a = rbind(jsonData$shapes$points[[i]], jsonData$shapes$points[[j]][n:1,], jsonData$shapes$points[[i]][1, ])
	plot(a, col = "gray", cex = 0.5, pch = 16, axes = FALSE, xlab = "", ylab = "", frame.plot = FALSE)
	polygon(a, col = "red", border = NA)
	dev.off()
}