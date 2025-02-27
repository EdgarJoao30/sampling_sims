library(sf)
library(terra)
library(tidyverse)
library(tidyterra)
library(stars)
library(raster)
library(INLA)
library(inlabru)
library(fmesher)
wd <- '~/OneDrive - University of Glasgow/PhD/0_simulations'
boundary <- st_read(paste0(wd, '/data/20240312_ROI_4326.shp')) |> 
  st_union() |> 
  st_transform(crs = 32650)
# Land cover
landcover <- rast(paste0(wd, '/data/aligned_landcover.tif')) 


# Define mosquito class
setClass(
  "Mosquito",
  slots = list(
    lc_coefficients = "numeric",
    spatial_range = "numeric",
    rho = 'numeric',
    theta = 'numeric',
    land_cover = 'SpatRaster',
    simulated_surface = 'SpatRaster',
    seed = 'numeric'
  )
)

setGeneric("_get_theta", function(object) {
  standardGeneric("_get_theta")
})

setMethod("_get_theta", "Mosquito", function(object) {
  alpha <- 2
  sigma <- 0.5
  variance <- sigma^2
  range <-object@spatial_range
  n <- 12
  kappa <- sqrt(8 * (alpha - 1)) / range
  object@theta <- c(-0.5 * log(4 * pi * variance * kappa^2), log(kappa))
  return(object)  # Return the modified object 
})

setGeneric("_get_sims", function(object) {
  standardGeneric("_get_sims")
})

setMethod("_get_sims", "Mosquito", function(object) {
  result <- generate_sims(boundary,
                          object@land_cover,
                          theta = object@theta, 
                          rho = object@rho, 
                          beta = object@lc_coefficients, 
                          k = k, 
                          sd_mu = sd_mu,
                          seed = seed) 
  
  object@simulated_surface <- result$raster
  return(object)  # Return the modified object 
})

# Constructor function for the Mosquito class
Mosquito <- function(lc_coefficients, spatial_range, rho, land_cover) {
  obj <- new("Mosquito", lc_coefficients = lc_coefficients, 
             spatial_range = spatial_range,
             rho = rho,
             land_cover = land_cover)
  
  obj <- `_get_theta`(obj)  # Automatically apply the _get_theta method
  obj <- `_get_sims`(obj) 
  return(obj)
}


beta <- c(0, 0.4567584, 1.5648494, 1.0986123, 1.8282377)
range <- 3000
rho <- 0.3

anopheles <- Mosquito(lc_coefficients = beta, 
                      spatial_range = range,
                      rho = rho,
                      land_cover = landcover)

plot(anopheles@simulated_surface)
