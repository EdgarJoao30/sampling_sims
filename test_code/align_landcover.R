library(sf)
library(terra)
library(tidyverse)
library(tidyterra)
library(stars)
library(raster)
wd <- '~/OneDrive - University of Glasgow/PhD/0_simulations'
boundary <- st_read(paste0(wd, '/data/20240312_ROI_4326.shp')) |> 
  st_union() |> 
  st_transform(crs = 32650)
boundary_sp <- as(boundary, "Spatial")
# Land cover
landcover <- rast(paste0(wd, '/data/Landcover_AllClass.tif')) %>% project(crs(boundary)) %>% raster::crop(boundary)
boundary_extent <- extent(boundary_sp)
empty_raster <- raster(boundary_extent, res = 500, crs = crs(boundary_sp))
aligned_landcover <- terra::resample(landcover, rast(empty_raster), method = "near")
mask <- rasterize(vect(boundary), aligned_landcover)
hab_mask <- terra::mask(aligned_landcover, mask)
hab_mask <- hab_mask %>% filter(Landcover_AllClass < 99)
hab_mask <- round(hab_mask)

writeRaster(hab_mask, paste0(wd, '/data/aligned_landcover.tif'), overwrite=TRUE)
