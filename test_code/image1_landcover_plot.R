library(tidyverse)
library(terra)
library(sf)
library(tidyterra)
wd <- '~/OneDrive - University of Glasgow/PhD/0_simulations'
boundary <- st_read(paste0(wd, '/data/20240312_ROI_4326.shp')) |> 
  st_union() |> 
  st_transform(crs = 32650)
boundary_sp <- as(boundary, "Spatial")
landcover <- rast(paste0(wd, '/data/Landcover_AllClass.tif')) %>% project(crs(boundary)) %>% raster::crop(boundary)
boundary_extent <- extent(boundary_sp)
empty_raster <- raster(boundary_extent, res = 500, crs = crs(boundary_sp))
aligned_landcover <- terra::resample(landcover, rast(empty_raster), method = "near")
mask <- rasterize(vect(boundary), aligned_landcover)
hab_mask <- terra::mask(aligned_landcover, mask)
hab_mask <- hab_mask %>% filter(Landcover_AllClass < 99)
hab_mask <- round(hab_mask)

(lcplot <- ggplot() +
  geom_spatraster(data = as.factor(hab_mask)) +
  scale_fill_manual(
    values = c('#9fa86a', 'green', 'darkgreen', '#def016', 'gray'),
    labels = c("0" = "Oil palm plantations", "1" = "Secondary forest", "2" = "Primary forest", "3" = "Other plantations", "4" = "Built-up"),
    na.value = NA,
    na.translate = FALSE,
    guide = guide_legend(ncol = 2)
  ) +
  labs(fill = 'Land cover', x='', y='')+
  theme_minimal()+
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      legend.position = "bottom",
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 16)
    ) +
  ggspatial::annotation_scale(location = "bl", line_width = 0.6, tick_height = 0.45, text_cex = 0.6,
                   pad_y=unit(2,"mm")) +
  ggspatial::annotation_north_arrow(location = "bl", height = unit(1.2,"cm"),width = unit(1.2,"cm"),
                         pad_y=unit(6,"mm"),
                         style = ggspatial::north_arrow_nautical(line_width = 0.6, text_size = 6)) 
)
ggsave('~/Documents/GitHub/sampling_sims/images/image1_landcover.jpeg', lcplot, dpi = 300)
