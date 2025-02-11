#################################################### Extract values from habitat and simulations ###############################################################
library(sf)
library(tidyverse)
library(terra)
library(tidyterra)
library(tmap)
library(raster)
wd <- '~/OneDrive - University of Glasgow/PhD/0_simulations'

df <- st_read(paste0(wd, '/data/20250129_points_sampling_scenarios_nobuffer.geojson'))
boundary <- st_read(paste0(wd, '/data/20240312_ROI_4326.shp')) |> st_union() |> st_transform(crs = 32650) 
boundary_sp <- as(boundary, "Spatial")
r <- rast(paste0(wd, '/data/20250202_sim_raster001.tif'))
landcover <- rast(paste0(wd, '/data/Landcover_AllClass.tif')) %>% project(crs(boundary)) %>% raster::crop(boundary)
boundary_extent <- extent(boundary_sp)
empty_raster <- raster::raster(boundary_extent, res = 500, crs = crs(boundary_sp))
aligned_landcover <- terra::resample(landcover, rast(empty_raster), method = "near")
mask <- rasterize(vect(boundary), aligned_landcover)
hab_mask <- terra::mask(aligned_landcover, mask)
hab_mask <- hab_mask %>% filter(Landcover_AllClass < 99)
hab_mask <- round(hab_mask)

# g <- st_read(paste0(wd, '/data/20250129_grids_sampling_scenarios_nobuffer.geojson'))
# g <- g %>% dplyr::filter(month == m, iteration == 1, sample_size == 5, scenario == 'a') 
# 
# dfs <- list()
# m = 1
# df_m <- df %>% dplyr::filter(month == m, iteration == 1, sample_size == 5, scenario == 'a') 
# df_m <- df_m %>% mutate(ID = 1:nrow(df_m))
# df_h <- terra::extract(hab_mask, df_m, xy = T, na.rm = T) 
# df_h <- st_as_sf(df_h, coords = c('x', 'y'), crs = st_crs(32650))
# 
# r_m <- r[[m]]
# df_e <- terra::extract(r_m, df_m, xy = T)
# 
# category_colors <- c("0" = "blue", "1" = "green", "2" = "yellow", "3" = "orange", "4" = "purple")
# 
# ggplot() +
#   #geom_spatraster(data = r[[1]]) +
#   geom_spatraster(data = as.factor(hab_mask)) +
#   scale_fill_manual(values = category_colors, na.value = "transparent") +
#   geom_sf(data = boundary, col = 'red', alpha = .5) +
#   geom_sf(data = g, col = 'red', alpha = .5) +
#   geom_sf(data = df_m) +
#   geom_sf(data = df_h, col = 'red') +
#   geom_text(data = df_h, aes(label = ID, geometry = geometry), stat = "sf_coordinates", size = 5, color = "black") +
#   geom_text(data = df_m, aes(label = ID, geometry = geometry), stat = "sf_coordinates", size = 5, color = "red") 

df_h <- terra::extract(hab_mask, df) 
df$cat_500m <- df_h$Landcover_AllClass

dfs <- list()

for (m in 1:12) {
  df_m <- df %>% dplyr::filter(month == m)
  r_m <- r[[m]]
  
  df_e <- terra::extract(r_m, df_m)
  colnames(df_e) <- c('ID', 'sim_anoph')
  df_m$sim_anoph <- df_e$sim_anoph
  dfs[[m]] <- df_m
}

df_all <- do.call(rbind, dfs)

st_write(df_all, paste0(wd, '/data/20250205_points_sampling_scenarios_alldata.geojson'), append = F)
