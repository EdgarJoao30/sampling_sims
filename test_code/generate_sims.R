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
# Land cover
landcover <- rast(paste0(wd, '/data/aligned_landcover.tif')) 


k = 12
rho = 0.3
beta <- c(0, 0.4567584, 1.5648494, 1.0986123, 1.8282377)
sigma <- 0.5
variance <- sigma^2
alpha <- 2
range <- 3000
kappa <- sqrt(8 * (alpha - 1)) / range
theta <- c(-0.5 * log(4 * pi * variance * kappa^2), log(kappa))
seed <- 1234
sd_mu <- 0.001
phi <- 10 # overdispersion

generate_sims <- function(boundary = boundary,
                         landcover = landcover,
                         alpha = 2, 
                         theta = theta, 
                         rho = rho, 
                         beta = beta, 
                         k = k, # number of samples
                         sd_mu = sd_mu,
                         phi = phi,
                         seed = seed) {
  
  points <- as.data.frame(landcover, xy=T)
  # colnames(points)[3] <- 'class'
  # points$class <- factor(points$class)
  # points$class <- relevel(factor(points$class), ref = 2)
  points$class <- factor(points$Landcover_AllClass)
  points$class <- plyr::revalue(points$class, c("0"="Oil", "1"="Secondary", "2"="Primary", "3"="Plantation", "4"="Built"))
  points$class <- relevel(factor(points$class), ref = "Primary")
  points <- st_as_sf(points, coords = c('x', 'y'), crs = 32650)
  mask <- rasterize(vect(boundary), landcover)
  # Mesh and true surface, units = meters
  mesh <- fm_mesh_2d(points, max.edge = c(2500, 5000), cutoff = 1000)
  spde <- inla.spde2.matern(mesh, alpha = alpha)
  Q <- inla.spde2.precision(spde, theta = theta)
  true_field <- inla.qsample(k, Q, seed = seed)
  points$field <- fm_evaluate(mesh, loc = points, field = true_field)
  # Compute AR1
  points$field_AR1 <- points$field
  for (j in 2:k) {
    points$field_AR1[, j] <- rho * points$field_AR1[, j - 1] + sqrt(1 - rho^2) * points$field[, j]
  }
  # Add regression covariates
  ccov <- factor(replicate(k, points$class))
  n <- nrow(points)
  mu <- beta[unclass(ccov)] + points$field_AR1 + rnorm(n * k, 0, sd_mu)
  points$mu <- exp(mu)
  # Draw from nbinomial distribution
  generate_nbinomial <- function(x) {
    rnbinom(mu = x, n = 1, size = phi)
  }
  set.seed(seed)
  nbinomial_sample <- apply(points$mu, c(1, 2), generate_nbinomial)
  points$mosq <- nbinomial_sample
  points <- cbind(points, as.data.frame(points$mosq)) %>% dplyr::select(sample.1:sample.12)
  # Convert points to raster
  surface <- st_rasterize(points) %>% rast() %>% terra::mask(mask)
  names(surface) <- c('jan', 'feb', 'mar', 'apr', 'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec')
  
  return(list(df = points, raster = df_mask))
}


# Example usage:
result <- generate_data(boundary,
                        landcover,
                        seed = 0, 
                        alpha = alpha, 
                        theta = theta, 
                        rho = 0.3, 
                        beta = beta, 
                        k = k, 
                        sd_mu = sd_mu) 
df_sf <- result$df
df_raster <- result$raster
rasters <- c(df_raster, df_raster_2)
writeRaster(df_raster, paste0(wd, '/data/20250212_sim_raster001.tif'), overwrite=TRUE)

