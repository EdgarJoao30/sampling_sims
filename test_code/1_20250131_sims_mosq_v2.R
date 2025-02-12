library(INLA)
library(inlabru)
library(fmesher)
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

generate_data <- function(boundary,
                          aligned_landcover,
                          df_sf, 
                          seed = seed, 
                          alpha = alpha, 
                          theta = theta, 
                          rho = rho, 
                          beta = beta, 
                          k = k, 
                          sd_mu = sd_mu) {
  boundary_sp <- as(boundary, "Spatial")
  mask <- rasterize(vect(boundary), aligned_landcover)
  # Mesh and true surface, units = meters
  mesh <- fm_mesh_2d(df_sf, max.edge = c(2500, 5000), cutoff = 1000)
  spde <- inla.spde2.matern(mesh, alpha = alpha)
  Q <- inla.spde2.precision(spde, theta = theta)
  
  true_field <- inla.qsample(n, Q, seed = seed)
  
  df_sf$field <- fm_evaluate(mesh, loc = df_sf, field = true_field)
  
  # Compute AR1
  df_sf$field_AR1 <- df_sf$field 
  for (j in 2:n) {
    df_sf$field_AR1[, j] <- rho * df_sf$field_AR1[, j - 1] + sqrt(1 - rho^2) * df_sf$field[, j]  
  }
  
  # Add regression covariates
  ccov <- factor(replicate(k, df_sf$class))
  n <- nrow(df)
  mu <- beta[unclass(ccov)] + df_sf$field_AR1 + rnorm(n * k, 0, sd_mu)
  df_sf$mu <- exp(mu)
  
  generate_binomial <- function(x) {
    rnbinom(mu = x, n = 1, size = 10)
  }
  
  binomial_sample <- apply(df_sf$mu, c(1, 2), generate_binomial)
  
  df_sf$mosq <- binomial_sample
  # df_sf$mosq[is.na(df_sf$mosq)] <- 0
  
  
  df2 <- cbind(df_sf, as.data.frame(df_sf$mosq)) %>% dplyr::select(sample.1:sample.12)
  # Convert df_sf to raster
  df_raster <- st_rasterize(df2)
  df_raster_terra <- rast(df_raster)
  df_mask <- terra::mask(df_raster_terra, mask) 
  names(df_mask) <- c('jan', 'feb', 'mar', 'apr', 'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec')
  
  return(list(df = df_sf, raster = df_mask))
}

# Example usage:
result <- generate_data(boundary,
                        aligned_landcover,
                        df_sf, 
                        seed = seed, 
                        alpha = alpha, 
                        theta = theta, 
                        rho = rho, 
                        beta = beta, 
                        k = k, 
                        sd_mu = sd_mu) 
df_sf <- result$df
df_raster <- result$raster
writeRaster(df_raster, paste0(wd, '/data/20250212_sim_raster001.tif'), overwrite=TRUE)
