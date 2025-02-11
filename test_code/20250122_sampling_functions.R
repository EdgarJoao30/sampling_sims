#################################################### Functions for sampling every month ###############################################################
library(sf)
library(sp)
library(raster)
library(tidyverse)

# wd <- '~/OneDrive - University of Glasgow/PhD/0_simulations'
# roi <- st_read(paste0(wd, '/data/20240312_ROI_4326.shp')) |> 
#   st_transform(crs = 32650)


############## 
############## Function to sample randomly and fixed them for the whole year 
############## 
fixed_random_sample_per_month <- function(roi, sample_size = 15, points_per_sample = 3) {
  
  all_grids <- list()
  all_samples <- list()
  
  sampled_indices <- sample(1:nrow(roi), sample_size)
  sampled_grids <- roi[sampled_indices, ]
  
  for (month in 1:12) {
    
    month_samples <- list()
    
    for (i in 1:nrow(sampled_grids)) {
      
      pts <- spsample(as(sampled_grids[i, ], "Spatial"), n = points_per_sample, type = "random")
      pts_df <- as.data.frame(pts)
      pts_df$month <- month
      pts_df$sample_id <- i
      month_samples[[i]] <- pts_df
    }
    
    grid_samples <- sampled_grids
    grid_samples$month <- month
    all_grids[[month]] <- grid_samples
    
    all_samples[[month]] <- do.call(rbind, month_samples)
  }
  
  return(list(do.call(rbind, all_samples), do.call(rbind, all_grids)))
}

# test <- fixed_random_sample_per_month(roi)
# 
# test_samples <- st_as_sf(test[[1]], coords = c('x', 'y'), crs = 32650)
# test_grids <- test[[2]]
# 
# ggplot() +
#   #geom_sf(data = roi) +
#   geom_sf(data = test_grids) +
#   geom_sf(data = test_samples, aes(color = month))

############## 
############## Function to sample randomly from the ROI for each month
############## 
random_sample_per_month <- function(roi, sample_size = 15, points_per_sample = 3) {
  
  all_grids <- list()
  all_samples <- list()
  
  for (month in 1:12) {
    sampled_indices <- sample(1:nrow(roi), sample_size)
    sampled_grids <- roi[sampled_indices, ]
    month_samples <- list()
    
    for (i in 1:nrow(sampled_grids)) {

      pts <- spsample(as(sampled_grids[i, ], "Spatial"), n = points_per_sample, type = "random")
      pts_df <- as.data.frame(pts)
      pts_df$month <- month
      pts_df$sample_id <- i
      month_samples[[i]] <- pts_df
    }
    
    sampled_grids$month <- month
    all_grids[[month]] <- sampled_grids
    
    all_samples[[month]] <- do.call(rbind, month_samples)
  }
  
  return(list(do.call(rbind, all_samples), do.call(rbind, all_grids)))
}

############## 
############## Function to make a stratified sample using the categories for each month
############## 

stratified_sample_per_month <- function(roi, sample_size = 15, points_per_sample = 3) {
  
  all_grids <- list()
  all_samples <- list()
  
  for (month in 1:12) {
    categories <- unique(roi$cat)
    samples_per_category <- ceiling(sample_size / length(categories))
    
    month_samples <- list()
    month_grids <- list()
    
    for (cat in categories) {
      cat_shp <- roi[roi$cat == cat, ]
      sampled_indices <- sample(1:nrow(cat_shp), min(samples_per_category, nrow(cat_shp)))
      cat_samples <- cat_shp[sampled_indices, ]
      
      for (i in 1:nrow(cat_samples)) {
        pts <- spsample(as(cat_samples[i, ], "Spatial"), n = points_per_sample, type = "random")
        pts_df <- as.data.frame(pts)
        pts_df$month <- month
        pts_df$sample_id <- paste(cat, i, sep = "_")
        month_samples[[length(month_samples) + 1]] <- pts_df
        
      }
      
      cat_samples$month <- month
      month_grids[[length(month_grids) + 1]] <- cat_samples
      #month_grids$month <- month
    }
    
    all_grids[[month]] <- do.call(rbind, month_grids)
    all_samples[[month]] <- do.call(rbind, month_samples)
  }
  
  return(list(do.call(rbind, all_samples), do.call(rbind, all_grids)))
}

############## 
############## Function to get samples equally spaced in the overall geometry of all the polygons together for each month
############## 

equally_spaced_sample_per_month <- function(roi, sample_size = 15, points_per_sample = 3) {

  all_samples <- list()
  
  dissolved <- roi %>% st_union() 
  
  repeat{
    samples_loc <- st_sample(dissolved, size = sample_size, type = 'regular')
    if(length(samples_loc) == sample_size){
      break
    }
  }
  
  joined <- st_join(roi, st_as_sf(samples_loc), join = st_intersects, left = F)
  
  grid_samples <- list()
  for (month in 1:12) {
    month_samples <- list()
    for (i in 1:nrow(joined)) {
      pts <- spsample(as(joined[i, ], "Spatial"), n = points_per_sample, type = "random")
      pts_df <- as.data.frame(pts)
      pts_df$month <- month
      pts_df$sample_id <- i
      month_samples[[i]] <- pts_df
    }
    
    # Add a new column with the iteration number (month)
    grid_month <- joined
    grid_month$month <- month
    
    # Append the copy to the list
    grid_samples[[month]] <- grid_month
    
    all_samples[[month]] <- do.call(rbind, month_samples)
  }
  
  
  return(list(do.call(rbind, all_samples), do.call(rbind, grid_samples)))
}
