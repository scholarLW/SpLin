generate_grid_numbers <- function(x_grid, y_grid)
{
  grid_rows <- length(y_grid) - 1
  grid_cols <- length(x_grid) - 1
  grid_numbers <- matrix(0, nrow = grid_rows, ncol = grid_cols)
  
  output_matrix <- matrix(0, nrow = grid_rows * grid_cols, ncol = 3)
  output_index <- 1  
  
  number <- 1
  for (j in 1:grid_cols) {  
    if (j %% 2 == 1) {  
      for (i in grid_rows:1) {
        grid_numbers[i, j] <- number
        output_matrix[output_index, ] <- c(j, grid_rows - i + 1, number)
        number <- number + 1
        output_index <- output_index + 1
      }
    } else {  
      for (i in 1:grid_rows) {
        grid_numbers[i, j] <- number
        output_matrix[output_index, ] <- c(j, grid_rows - i + 1, number)
        number <- number + 1
        output_index <- output_index + 1
      }
    }
  }
  colnames(output_matrix) <- c("x_grid_id", "y_grid_id", "window")
  return(as.data.frame(output_matrix))
}

wcellDueGrid = function(rds, Assay = 'Spatial', grid_density = 30)
{
  suppressMessages(library(Seurat))
  options(stringsAsFactors = FALSE)
  obj = readRDS(rds)
  absolute_path <- paste0(dirname(normalizePath(rds)),"/IGSI/winCellGrid")
  dir.create(absolute_path, showWarnings = FALSE, recursive = TRUE, mode = "0777")
  
  metadata = obj@meta.data
  if(grid_density < 10 || grid_density > 1000 || !isTRUE(all.equal(round(grid_density), grid_density)) || is.na(grid_density))
  {
    stop("Grid density for splitting the signal matrix must be a positive integer between 10 and 1000. Please choose a value in this range.")
  }
  
  x_range <- range(metadata$x)
  y_range <- range(metadata$y)
  grid_length <- ceiling((x_range[2] - x_range[1]) / grid_density)
  grid_width <- ceiling((y_range[2] - y_range[1]) / grid_density)
  
  x_grid <- seq(x_range[1], x_range[1] + grid_density * grid_length, by = grid_length)
  y_grid <- seq(y_range[1], y_range[1] + grid_density * grid_width, by = grid_width)
  window_labels = generate_grid_numbers(x_grid, y_grid)
  
  assign_to_grid <- function(x, y, x_grid, y_grid, window_labels) 
  {
    x_index <- findInterval(x, x_grid, rightmost.closed = TRUE)
    y_index <- findInterval(y, y_grid, rightmost.closed = TRUE)
    ix <- which(window_labels$x_grid_id == x_index & window_labels$y_grid_id == y_index)
    if (length(ix) > 0) 
    {
      ixn <- window_labels[ix[1], ]  
      return(c(ixn$x_grid_id, ixn$y_grid_id, ixn$window))  
    } else {
      return(c(NA, NA, NA))
    }
  }
  res <- mapply(assign_to_grid, 
                metadata$x, metadata$y, 
                MoreArgs = list(x_grid = x_grid, y_grid = y_grid, 
                                window_labels = window_labels))
  
  res_matrix <- t(simplify2array(res))								 
  colnames(res_matrix) = c('x_grid_id','y_grid_id','window')							 
  metadata = as.data.frame(cbind(metadata, res_matrix))
  ix = is.na(metadata$window)
  metadata = metadata[!ix,]
  if(sum(table(metadata$window)>0)>100000)
  {
    stop("Grid density is too high.")
  }
  if(sum(table(metadata$window)>0)<20)
  {
    stop("Grid density is too small.")
  }	
  
  metadata$xnew = metadata$x - grid_length * (metadata$x_grid_id - 1) + grid_length * (metadata$window - 1)
  metadata$ynew = metadata$y - grid_width * (metadata$y_grid_id - 1)
  print(head(metadata))	
  
  obj = subset(obj, cells = rownames(metadata))
  obj@meta.data = metadata[rownames(obj@meta.data), ]
  count_matrix <- GetAssayData(obj, assay = Assay, slot = "counts")
  count_matrix = count_matrix[, rownames(metadata)]
  data = apply(count_matrix, 1, calculate_all_GFai_2, colnames(count_matrix), metadata, grid_density*grid_density, maxWindow = 100*100)	
  rownames(data) = paste("Window", sort(unique(metadata$window)), sep = "")
  data = t(data)
  objwd <- CreateSeuratObject(counts = data, assay = "Spatial", min.cells = 1, min.features = 1)
  saveRDS(objwd, file = paste0(absolute_path, "/wcellDueGrid.RDS"))
  saveRDS(obj, file = paste0(absolute_path, "/SpLin_EXP_output_withWindowGrid.RDS"))
}
