merge_intersecting_lists <- function(setlist)
{
  merged_list <- list()
  
  merge_to_existing <- function(current_set, merged_list) {
    for (j in seq_along(merged_list)) {
      if (length(intersect(current_set, merged_list[[j]])) > 0) {
        merged_list[[j]] <- unique(c(merged_list[[j]], current_set))
        return(merged_list)
      }
    }
    return(merged_list)
  }
  
  for (i in seq_along(setlist)) {
    current_set <- setlist[[i]]
    merged_list <- merge_to_existing(current_set, merged_list)
    
    if (length(merged_list) == 0 || all(sapply(merged_list, function(x) length(intersect(x, current_set)) == 0))) {
      merged_list <- c(merged_list, list(current_set))
    }
  }
  
  repeat {
    new_merged_list <- list()
    merged <- FALSE
    
    for (j in seq_along(merged_list)) {
      current_set <- merged_list[[j]]
      new_merged_list <- merge_to_existing(current_set, new_merged_list)
      
      if (length(new_merged_list) == 0 || all(sapply(new_merged_list, function(x) length(intersect(x, current_set)) == 0))) {
        new_merged_list <- c(new_merged_list, list(current_set))
      }
    }
    
    if (identical(merged_list, new_merged_list)) {
      break
    }
    
    merged_list <- new_merged_list
  }
  
  return(merged_list)
}

process_regionsGrid <- function(data, gridMatrix)
{
  listW = list()
  allN = unique(data$region)
  ymax = max(gridMatrix$y_grid_id)
  xmax = max(gridMatrix$x_grid_id)
  
  
  for(i in 1:length(allN))
  {
    ix = data$region == allN[i]
    hotW = unique(data$window[ix])
    tmp = hotW
    for(j in 1:length(hotW))
    {
      ix = gridMatrix$window == hotW[j]
      iid = gridMatrix$x_grid_id[ix]
      jid = gridMatrix$y_grid_id[ix]
      if(jid < ymax)
      {
        tmp = c(tmp, gridMatrix$window[gridMatrix$x_grid_id == iid & gridMatrix$y_grid_id == jid + 1])
      }		
      if(jid > 1)
      {
        tmp = c(tmp, gridMatrix$window[gridMatrix$x_grid_id == iid & gridMatrix$y_grid_id == jid - 1])
      }		
      if(iid > 1)
      {
        tmp = c(tmp, gridMatrix$window[gridMatrix$x_grid_id == iid - 1 & gridMatrix$y_grid_id == jid])
      }		
      if(iid < xmax)
      {
        tmp = c(tmp, gridMatrix$window[gridMatrix$x_grid_id == iid + 1 & gridMatrix$y_grid_id == jid])
      }		
      if(iid > 1 & jid > 1)
      {
        tmp = c(tmp, gridMatrix$window[gridMatrix$x_grid_id == iid - 1 & gridMatrix$y_grid_id == jid - 1])
      }			
      if(iid > 1 & jid < ymax)
      {
        tmp = c(tmp, gridMatrix$window[gridMatrix$x_grid_id == iid - 1 & gridMatrix$y_grid_id == jid + 1])
      }			
      if(iid < xmax & jid > 1)
      {
        tmp = c(tmp, gridMatrix$window[gridMatrix$x_grid_id == iid + 1 & gridMatrix$y_grid_id == jid - 1])
      }			
      if(iid < xmax & jid < ymax)
      {
        tmp = c(tmp, gridMatrix$window[gridMatrix$x_grid_id == iid + 1 & gridMatrix$y_grid_id == jid + 1])
      }
    }
    listW <- append(listW, list(unique(tmp)))
  }	
  setlist <- merge_intersecting_lists(listW)
  
  return(setlist)
}