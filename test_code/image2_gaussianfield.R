library(INLA)
library(inlabru)
library(fmesher)
library(sf)
library(terra)
library(tidyverse)
library(tidyterra)
library(stars)
library(raster)
library(ggridges)

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
df <- as.data.frame(hab_mask, xy=T)
df$class <- factor(df$Landcover_AllClass)
df$class <- plyr::revalue(df$class, c("0"="Oil", "1"="Secondary", "2"="Primary", "3"="Plantation", "4"="Built"))
df$class <- relevel(factor(df$class), ref = "Primary")
df_sf <- st_as_sf(df, coords = c('x', 'y'), crs = 32650)


k = 12
rho = 0.3
beta <- c(0, 0.4567584, 1.5648494, 1.0986123, 1.8282377)
sigma <- 0.5
variance <- sigma^2
alpha <- 2
range <- 3000
n <- 12
kappa <- sqrt(8 * (alpha - 1)) / range
theta <- c(-0.5 * log(4 * pi * variance * kappa^2), log(kappa))
seed <- 1234
sd_mu <- 0.001

boundary_sp <- as(boundary, "Spatial")
mask <- rasterize(vect(boundary), aligned_landcover)
# Mesh and true surface, units = meters
mesh <- fm_mesh_2d(df_sf, max.edge = c(2500, 5000), cutoff = 1000)
spde <- inla.spde2.matern(mesh, alpha = alpha)
Q <- inla.spde2.precision(spde, theta = theta)

true_field <- inla.qsample(n, Q, seed = seed)

df_sf$field <- fm_evaluate(mesh, loc = df_sf, field = true_field)

df_sf$field_AR1 <- df_sf$field
for (j in 2:k) {
  df_sf$field_AR1[, j] <- rho * df_sf$field_AR1[, j - 1] + sqrt(1 - rho^2) * df_sf$field[, j]
}  

colsc <- function(...) {
  scale_fill_gradientn(
    colours = rev(RColorBrewer::brewer.pal(11, "RdYlBu")),
    limits = range(..., na.rm = TRUE)
  )
}

df2 <- as.data.frame(df_sf$field_AR1)
df2 <- cbind(df_sf$class, df2, df_sf$geometry) %>% st_as_sf()
colnames(df2)[1] <- 'class'
df2 <- gather(df2, key = 'month', value = 'value', -c(class, geometry))
df2$month <- as.integer(gsub("sample:", "", df2$month))


csc <- colsc(df_sf$field[, 1])

(g <- ggplot() +
    gg(df_sf, aes(fill = field[, 1]), geom = "tile") +
    labs(fill = 'Value', x='', y='')+
    theme_minimal()+
    theme(axis.text.x = element_blank(), axis.text.y = element_blank())+
    csc +
    ggspatial::annotation_scale(location = "bl", line_width = 0.6, tick_height = 0.45, text_cex = 0.6,
                                pad_y=unit(2,"mm")) +
    ggspatial::annotation_north_arrow(location = "bl", height = unit(1.2,"cm"),width = unit(1.2,"cm"),
                                      pad_y=unit(6,"mm"),
                                      style = ggspatial::north_arrow_nautical(line_width = 0.6, text_size = 6)) 
)

csc2 <- colsc(df2$value)
(g2 <- ggplot() +
    gg(df2, aes(fill = value), geom = "tile") +
    labs(fill = 'Value', x='', y='')+
    facet_wrap(~month, ncol = 3)+
    theme_minimal()+
    theme(
      axis.text.x = element_blank(), 
      axis.text.y = element_blank(),
      strip.text = element_text(size = 12, face = "bold"))+
    csc2
)

ggsave('~/Documents/GitHub/sampling_sims/images/image2_gaussianfield.jpeg', g, dpi = 300)
ggsave('~/Documents/GitHub/sampling_sims/images/image3_gaussianfield_ar1.jpeg', g2, dpi = 300)


# Add regression covariates
ccov <- factor(replicate(k, df_sf$class))
n <- nrow(df)
mu <- beta[unclass(ccov)] + df_sf$field_AR1 + rnorm(n * k, 0, sd_mu)
df_sf$mu <- exp(mu)

generate_nbinomial <- function(x) {
  rnbinom(mu = x, n = 1, size = 10)
}

nbinomial_sample <- apply(df_sf$mu, c(1, 2), generate_nbinomial)

df_sf$mosq <- nbinomial_sample

df3 <- as.data.frame(df_sf$mosq)
df3 <- cbind(df_sf$class, df3, df_sf$geometry) %>% st_as_sf()
colnames(df3)[1] <- 'class'
df3 <- gather(df3, key = 'month', value = 'value', -c(class, geometry))
df3$month <- as.integer(gsub("sample:", "", df2$month))

csc3 <- colsc(df3$value)
(g3 <- ggplot() +
    gg(df3, aes(fill = value), geom = "tile") +
    labs(fill = 'Biting rate', x='', y='')+
    facet_wrap(~month, ncol = 3)+
    theme_minimal()+
    theme(
      axis.text.x = element_blank(), 
      axis.text.y = element_blank(),
      strip.text = element_text(size = 12, face = "bold"))+
    csc3
)
ggsave('~/Documents/GitHub/sampling_sims/images/image4_sim_mosq.jpeg', g3, dpi = 300)

(ridge_plot <- ggplot(df3, aes(x = value, y = class, fill = class)) +
    geom_density_ridges(scale = 1, alpha = 0.7) +
    scale_fill_manual(values = c('darkgreen', '#9fa86a', 'green', '#def016', 'gray')) +
    labs(
         x = "Biting Rate",
         y = "") +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      plot.title = element_text(size = 16, face = "bold")
    )
)

ggsave('~/Documents/GitHub/sampling_sims/images/image5_sim_density.jpeg', ridge_plot, dpi = 300)
