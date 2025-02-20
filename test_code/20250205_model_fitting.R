#################################################### Model fitting ###############################################################
library(sf)
library(tidyverse)
library(terra)
library(tidyterra)
library(tmap)
library(raster)
library(inlabru)
library(fmesher)
library(INLA)
wd <- '~/OneDrive - University of Glasgow/PhD/0_simulations'
# Landcover
boundary <- st_read(paste0(wd, '/data/20240312_ROI_4326.shp')) |> st_union() |> st_transform(crs = 32650) 
boundary_sp <- as(boundary, "Spatial")
landcover <- rast(paste0(wd, '/data/Landcover_AllClass.tif')) %>% project(crs(boundary)) %>% raster::crop(boundary)
boundary_extent <- extent(boundary_sp)
empty_raster <- raster::raster(boundary_extent, res = 500, crs = crs(boundary_sp))
aligned_landcover <- terra::resample(landcover, rast(empty_raster), method = "near")
mask <- rasterize(vect(boundary), aligned_landcover)
landcover_mask <- terra::mask(aligned_landcover, mask)
landcover_mask <- landcover_mask %>% filter(Landcover_AllClass < 99)
# Simulations
simulations <- rast(paste0(wd, '/data/20250212_sim_raster001.tif'))
simulations <- st_as_sf(as.data.frame(simulations, xy=T), crs = 32650, coords = c('x', 'y'))
lc_values <- terra::extract(landcover_mask, simulations, xy = T)
simulations$cat_500m <- round(lc_values$Landcover_AllClass)
simulations$cat_500m <- plyr::revalue(as.character(simulations$cat_500m), c("0"="Oil", "1"="Secondary", "2"="Primary", "3"="Plantation", "4"="Built"))
simulations$cat_500m <- relevel(factor(simulations$cat_500m), ref = "Primary")
land_category_dummies <- model.matrix(~ cat_500m - 1, data = simulations)
simulations <- cbind(simulations, land_category_dummies)

sim_long <- gather(simulations, key = 'month', value = 'sim', -c(cat_500m:geometry))
sim_long <- sim_long %>% 
  mutate(month = replace(month, month == 'jan', 1),
         month = replace(month, month == 'feb', 2),
         month = replace(month, month == 'mar', 3),
         month = replace(month, month == 'apr', 4),
         month = replace(month, month == 'may', 5),
         month = replace(month, month == 'jun', 6),
         month = replace(month, month == 'jul', 7),
         month = replace(month, month == 'aug', 8),
         month = replace(month, month == 'sep', 9),
         month = replace(month, month == 'oct', 10),
         month = replace(month, month == 'nov', 11),
         month = replace(month, month == 'dec', 12)
         )
sim_long$month <- as.numeric(sim_long$month)
# Samples
samples <- st_read(paste0(wd, '/data/20250219_points_sampling_scenarios_alldata.geojson'))
samples$cat_500m <- plyr::revalue(as.character(samples$cat_500m), c("0"="Oil", "1"="Secondary", "2"="Primary", "3"="Plantation", "4"="Built"))
samples$cat_500m <- relevel(factor(samples$cat_500m), ref = "Primary")
land_category_dummies <- model.matrix(~ cat_500m - 1, data = samples)
samples <- cbind(samples, land_category_dummies)

test <- samples %>% dplyr::filter(iteration == 1, scenario == 'c', sample_size == 15)

mesh <- fm_mesh_2d(simulations, max.edge = c(2500, 5000), cutoff = 1000)
matern <- inla.spde2.matern(mesh, alpha = 2, constr = T)

# model <- sim_anoph ~ land_cover(cat_500m, model = 'factor_contrast') +
#   #Intercept(1) + 
#   field(geometry, model = matern) +
#   time(month, model = "ar1")

model <- sim_anoph ~ -1 + 
  land_cover(cat_500m, model = 'factor_full') +
  field(geometry, model = matern) +
  time(month, model = "ar1")

fit <- bru(model, test, family = "nbinomial",
           options = list(control.family = list(link = "log"), 
                          control.compute = list(dic = TRUE, cpo = TRUE, config=T, dic = TRUE, waic = TRUE)
                          ))
pred <- predict(
  fit, sim_long,
  ~  exp(land_cover + field + time )
)

samp <- generate(fit, sim_long,
                 ~  exp(land_cover + field + time),
                 n.samples = 1
)

pred$sample <- samp[, 1]

summary(pred$sample)

pl_truth <- ggplot() +
  gg(pred %>% filter(month == 1), aes(fill = sim), geom = "tile") +
  facet_wrap( ~ month, nrow = 3) +
  gg(boundary_sp, alpha = 0) +
  ggtitle("Simulated")

pl_posterior_mean <- ggplot() +
  gg(pred %>% filter(month == 1), aes(fill = mean), geom = "tile") +
  facet_wrap( ~ month, nrow = 3) +
  gg(boundary_sp, alpha = 0) +
  ggtitle("Posterior mean")

pl_posterior_sample <- ggplot() +
  gg(pred %>% filter(month == 1), aes(fill = sample), geom = "tile") +
  facet_wrap( ~ month, nrow = 3) +
  gg(boundary_sp, alpha = 0) +
  ggtitle("Posterior sample")

# Common colour scale for the truth and estimate:
colsc <- function(...) {
  scale_fill_gradientn(
    colours = rev(RColorBrewer::brewer.pal(11, "RdYlBu")),
    limits = range(..., na.rm = TRUE)
  )
}

csc <- colsc(pred$sim, pred$mean, pred$sample)

multiplot(pl_truth + csc,
          pl_posterior_mean + csc,
          pl_posterior_sample + csc,
          cols = 3
)

spde.range <- spde.posterior(fit, "field", what = "range")
spde.logvar <- spde.posterior(fit, "field", what = "log.variance")
range.plot <- plot(spde.range)
var.plot <- plot(spde.logvar)

multiplot(range.plot, var.plot)

corplot <- plot(spde.posterior(fit, "field", what = "matern.correlation"))
covplot <- plot(spde.posterior(fit, "field", what = "matern.covariance"))
multiplot(covplot, corplot)


flist <- vector("list", NROW(fit$summary.random$land_cover))
for (i in seq_along(flist)) flist[[i]] <- plot(fit, "land_cover", index = i)
multiplot(plotlist = flist, cols = 3)

# test delete
library(MASS)
model_test <- glm.nb(sim_anoph ~ cat_500m, data = test)
predictions <- predict(model_test, newdata = sim_long, type = "response")

