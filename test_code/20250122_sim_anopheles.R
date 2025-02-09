library(INLA)
library(inlabru)
library(fmesher)
library(sf)
library(terra)
library(tidyverse)
library(mgcv)
library(tidyterra)
library(patchwork)
# data
wd <- '~/OneDrive - University of Glasgow/PhD/0_simulations'
source(paste0(wd, '/code/spde-book-functions.R'))
boundary <- st_read(paste0(wd, '/data/20240312_ROI_4326.shp')) |> 
  st_union() |> 
  st_transform(crs = 32650)
# Land cover
hab <- rast("~/OneDrive - University of Glasgow/PhD/0_simulations/data/Landcover_AllClass.tif") %>% project(crs(boundary)) %>% raster::crop(boundary)
mask <- rasterize(vect(boundary), hab)
hab_mask <- terra::mask(hab, mask)
hab_mask <- hab_mask %>% filter(Landcover_AllClass < 99)
## Create new raster
r <- rast(ext(hab_mask), resolution = 500)
r <- terra::resample(hab_mask, r)
r <- round(r)
df <- as.data.frame(r, xy=T)
df$class <- factor(df$Landcover_AllClass)
df$class <- plyr::revalue(df$class, c("0"="Oil", "1"="Secondary", "2"="Primary", "3"="Plantation", "4"="Built"))
df$class <- relevel(factor(df$class), ref = "Primary")
df_sf <- st_as_sf(df, coords = c('x', 'y'), crs = 32650)

### mesh and true surface, units = meters
mesh <- fm_mesh_2d(df_sf, 
                    max.edge = c(2500, 5000), 
                    cutoff = 1000)

g <- ggplot() +
  gg(mesh) +
  geom_sf(data = boundary, col = 'green', alpha = .5) +
  labs(x = '', y = '') +
  theme_minimal()+
  theme(axis.text.x = element_blank(), axis.text.y = element_blank())
# 
# ggsave('~/OneDrive - University of Glasgow/PhD/0_simulations/figures/mesh.jpg', g, dpi = 300)


sigma=0.5
variance=sigma^2
alpha=2
range = 1
n = 12
kappa = sqrt(8*(alpha-1))/range
theta <- c(-0.5*log(4*pi*variance*kappa^2), log(kappa))

spde <- inla.spde2.matern(mesh, alpha=alpha)
Q <- inla.spde2.precision(spde, theta=theta)

true_field <- inla.qsample(n, Q, seed = 1234)

df_sf$field <- fm_evaluate(
  mesh,
  loc = df_sf,
  field = true_field
)

## Set rho (value for temporal correlation)
rho <- 0.3
df_sf$field_AR1 <- NULL
df_sf$field_AR1 <- df_sf$field
## Compute AR1
for (j in 2:n) {
  df_sf$field_AR1[, j] <- rho * df_sf$field_AR1[, j - 1] + sqrt(1 - rho^2) * df_sf$field[, j]
}

ggplot() +
  gg(df_sf, aes(fill = field_AR1), geom = "tile") +
  coord_equal()
## Add regression covariates (from glm) *** in right order for sense 
beta <- c(0, 0.4567584, 1.5648494, 1.0986123, 1.8282377) # anoph - add in right order for categories
## Set number of time points
k <- 12
## Repeat habitat data per year
ccov <- factor(replicate(12, df_sf$class))
## Set number of cells
n <- nrow(df)
## Calculate mu from fixed effects + spatial effects
sd.mu <- 0.001
mu <- beta[unclass(ccov)] + df_sf$field_AR1 + rnorm(n * k, 0, sd.mu)

df_sf$mu <- exp(mu)

generate_binomial <- function(x) {
  rnbinom(mu= x, n=1, size=10) #### FLAG 2 - Change to poisson, to prevent competition with previous terms
}

# Apply the function to each element of the data frame
binomial_sample <- apply(df_sf$mu, c(1, 2), generate_binomial)

df_sf$mosq <- binomial_sample
df_sf$mosq[is.na(df_sf$mosq)] <- 0

# ggplot() +
#   gg(df_sf, aes(fill = mosq[,1]), geom = "tile") +
#   ggtitle("True field")+
#   theme_minimal()

plots <- list()
months <- c('January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December')

for (i in 1:12) {
  g <- ggplot() +
    gg(df_sf, aes(fill = mosq[,i]), geom = "tile") +
    ggtitle(months[i])+
    labs(fill = 'Counts', x='', y='')+
    theme_minimal()+
    scale_fill_gradient(low = "white", high = "red")+
    theme(axis.text.x = element_blank(), axis.text.y = element_blank())
  plots[[i]] <- g
}

(combined_plot <- wrap_plots(plots, ncol = 3, nrow = 4))

(g <- ggplot() +
  gg(df_sf, aes(fill = field_AR1[,12]), geom = "tile") +
  ggtitle('Anopheles count')+
  labs(fill = 'Value', x='', y='')+
  theme_minimal()+
  scale_fill_gradient(low = "white", high = "red")+
  theme(axis.text.x = element_blank(), axis.text.y = element_blank())
)
  # 
# (g <- ggplot() +
#   geom_spatraster(data = hab_mask) +
#   scale_fill_manual(values = c('red', 'blue', 'white', 'gray')) +
#   ggtitle('Land cover')+
#   labs(fill = 'Value', x='', y='')+
#   theme_minimal()+
#   theme(axis.text.x = element_blank(), axis.text.y = element_blank())
# )
# ggsave('~/OneDrive - University of Glasgow/PhD/0_simulations/figures/anopheles.jpg', g, dpi = 300)
