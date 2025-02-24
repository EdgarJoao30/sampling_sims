######################################################################################################################################################################################################################
######################################################################################################################################################################################################################
##########################################################################################################################################################

## POSTERIOR SAMPLES

##########################################################################
################################ ANOPHELES ###############################

# create raster and make shapefile to get coordinates back
roi  <- readOGR(dsn="C:/Users/Kimberly Fornace/OneDrive - University of Glasgow/SENSOR/R/SAFE_ento/", layer="ROI_extent")
roi <- spTransform(roi, CRS("+proj=utm +zone=50 +datum=WGS84 +units=km +no_defs +ellps=WGS84 +towgs84=0,0,0"))
rr <- raster(ext = extent(roi), resolution = 1)
values(rr) <- 1 # set all vaues to one
df_rr <- as.data.frame(rr, xy=T) # get df from raster for spatial points
xmin <- min(df_rr$x)
ymin <- min(df_rr$y)

#simulated
pred <- read.csv("C:/Users/Kimberly Fornace/OneDrive - University of Glasgow/SENSOR/R/SAFE_ento/simulated_anoph_500m.csv") # sim_mos
pred$hab <- relevel(factor(pred$hab), ref = "Primary")
pred$hab <- plyr::revalue(pred$hab, c("Oil"="oil", "Secondary"="secondary", "Primary"="primary", "Built"="built", "Plantation"="plantation"))
pred$xcoo <- pred$xcoo + xmin
pred$ycoo <- pred$ycoo + ymin

# vector of actual ('real' simulated) values 
pmos <- pred$mos

## Simulate 100 draws from the posterior negative binomial distribution
# Kim code 2023
pred.nb <- function(l, nsample=nrow(pred)){
  lp <- l$latent[1:nsample] ## extract the linear predictor in log scale
  lambda  <-  exp(lp)   ## apply the link function
  y <- rnbinom(nsample,mu=lambda, size=l$hyperpar[1]) ## posterior sample of non-zero inflated data
}

### Hab-strat

iteration <- c(1:74,76:100)

for(j in 1:length(iteration)){ 
  
  ### 15 sites ###
  load(paste0("C:/R/SENSOR/EJ_SAFE/outputs/INLA_Anopheles/Models_An/spde_ss4_it",iteration[j],"_n15.RData"))
  
  ## Posterior mean total number of mosquito bites per location
  Post.Res <- inla.posterior.sample(100, res)
  
  # Simulated draws from the negative binomial distribution
  samples.df <- sapply(Post.Res, pred.nb)
  
  # Summary values for posterior samples per row (per cell)
  sum.post <- apply(samples.df, 1, summary)
  sum.post <- as.data.frame(sum.post)
  sum.post <- t(sum.post)
  
  # SD values for posterior samples per row (per cell) 
  sd.post <- apply(samples.df, 1, sd)
  sd.post <- as.data.frame(round(sd.post, digits=4))
  
  # Join 
  output.post <- cbind(sum.post, sd.post)
  rm(sd.post); rm(sum.post)
  
  # Cell-wise (row-wise) RMSE of posteriors (per raster grid across space/time)
  for(i in 1:nrow(samples.df)) {
    
    predicted <- samples.df[i,]
    actual <- pmos[i]
    d <- (predicted - actual)^2
    rsd <- sqrt(sd(d))
    rmse <- rmse(actual, predicted)
    rmse.rsd <- cbind(rmse, rsd)
    #   rmse.rsd <- round(rmse.rsd, digits=4)
    
    if(exists("rmse.samp")){
      rmse.samp <- rbind(rmse.samp, rmse.rsd)
      
    }
    
    # If merged dataset does not exist, create dataset
    if(!exists("rmse.samp")){
      rmse.samp <- rmse.rsd
    }
    rm(rmse); rm(rsd); rm(rmse.rsd)
    
  } 
  
  output <- cbind(output.post, rmse.samp)
  rm(output.post);  rm(rmse.samp)
  
  
  # rbind vertical join 
  output$type <- "hab_fixed"
  output$num <- 15
  output$iter <- iteration[j]
  
  
  if(exists("output.all")){
    output.all <- rbind(output.all, output)
    
  }
  
  # If merged dataset does not exist, create dataset
  if(!exists("output.all")){
    output.all <- output
  }
  
  rm(output)
  
  
  # Col-wise overall RMSE per posterior sample 100
  for(i in 1:ncol(samples.df)) {
    
    predicted <- samples.df[,i]
    actual <- pmos
    rmse <- rmse(actual, predicted)
    
    if(exists("rmse.mod")){
      rmse.mod <- rbind(rmse.mod, rmse)
      
    }
    
    # If merged dataset does not exist, create dataset
    if(!exists("rmse.mod")){
      rmse.mod <- rmse
    }
    rm(rmse)
  } 
  
  rmse.df <- as.data.frame(rmse.mod)
  rmse.df <- t(rmse.df)
  rmse.df <- as.data.frame(rmse.df)
  rmse.df$iter <- iteration[j]
  rmse.df$type <- "hab_fixed"
  rmse.df$num <- 15
  
  if(exists("rmse.all")){
    rmse.all <- rbind(rmse.all, rmse.df)
    
  }
  
  # If merged dataset does not exist, create dataset
  if(!exists("rmse.all")){
    rmse.all <- rmse.df
  }
  rm(rmse.df); rm(rmse.mod)
  print(iteration[j])
  
}

write.csv(output.all,"C:/R/SENSOR/EJ_SAFE/outputs/INLA_Anopheles/Posteriors_CSV/An_output_posterior_ss4_n15.csv")
write.csv(rmse.all,"C:/R/SENSOR/EJ_SAFE/outputs/INLA_Anopheles/Posteriors_CSV/An_rmse_posterior_ss4_n15.csv")
#
rm(output.all); rm(rmse.all) 




# ss5 - geospatial #######################  

iteration <- 1:100

for(j in 1:length(iteration)){ 
  
  ### 5 sites ###
  load(paste0("C:/R/SENSOR/EJ_SAFE/outputs/INLA_Anopheles/Models_An/spde_ss5_it",iteration[j],"_n5.RData"))
  
  ## Posterior mean total number of mosquito bites per location
  Post.Res <- inla.posterior.sample(100, res)
  
  # Simulated draws from the negative binomial distribution
  samples.df <- sapply(Post.Res, pred.nb)
  
  # Summary values for posterior samples per row (per cell)
  sum.post <- apply(samples.df, 1, summary)
  sum.post <- as.data.frame(sum.post)
  sum.post <- t(sum.post)
  
  # SD values for posterior samples per row (per cell) 
  sd.post <- apply(samples.df, 1, sd)
  sd.post <- as.data.frame(round(sd.post, digits=4))
  
  # Join 
  output.post <- cbind(sum.post, sd.post)
  rm(sd.post); rm(sum.post)
  
  # Cell-wise (row-wise) RMSE of posteriors (per raster grid across space/time)
  for(i in 1:nrow(samples.df)) {
    
    predicted <- samples.df[i,]
    actual <- pmos[i]
    d <- (predicted - actual)^2
    rsd <- sqrt(sd(d))
    rmse <- rmse(actual, predicted)
    rmse.rsd <- cbind(rmse, rsd)
    #   rmse.rsd <- round(rmse.rsd, digits=4)
    
    if(exists("rmse.samp")){
      rmse.samp <- rbind(rmse.samp, rmse.rsd)
      
    }
    
    # If merged dataset does not exist, create dataset
    if(!exists("rmse.samp")){
      rmse.samp <- rmse.rsd
    }
    rm(rmse); rm(rsd); rm(rmse.rsd)
    
  } 
  
  output <- cbind(output.post, rmse.samp)
  rm(output.post);  rm(rmse.samp)
  
  
  # rbind vertical join 
  output$type <- "geospat"
  output$num <- 5
  output$iter <- iteration[j]
  
  
  if(exists("output.all")){
    output.all <- rbind(output.all, output)
    
  }
  
  # If merged dataset does not exist, create dataset
  if(!exists("output.all")){
    output.all <- output
  }
  
  rm(output)
  
  
  # Col-wise overall RMSE per posterior sample 100
  for(i in 1:ncol(samples.df)) {
    
    predicted <- samples.df[,i]
    actual <- pmos
    rmse <- rmse(actual, predicted)
    
    if(exists("rmse.mod")){
      rmse.mod <- rbind(rmse.mod, rmse)
      
    }
    
    # If merged dataset does not exist, create dataset
    if(!exists("rmse.mod")){
      rmse.mod <- rmse
    }
    rm(rmse)
  } 
  
  rmse.df <- as.data.frame(rmse.mod)
  rmse.df <- t(rmse.df)
  rmse.df <- as.data.frame(rmse.df)
  rmse.df$iter <- iteration[j]
  rmse.df$type <- "geospat"
  rmse.df$num <- 5
  
  if(exists("rmse.all")){
    rmse.all <- rbind(rmse.all, rmse.df)
    
  }
  
  # If merged dataset does not exist, create dataset
  if(!exists("rmse.all")){
    rmse.all <- rmse.df
  }
  rm(rmse.df); rm(rmse.mod)
  print(iteration[j])
  
} 

write.csv(output.all,"C:/R/SENSOR/EJ_SAFE/outputs/INLA_Anopheles/Posteriors_CSV/An_output_posterior_ss5_n5.csv")
write.csv(rmse.all,"C:/R/SENSOR/EJ_SAFE/outputs/INLA_Anopheles/Posteriors_CSV/An_rmse_posterior_ss5_n5.csv")
rm(output.all); rm(rmse.all)




for(j in 1:length(iteration)){ 
  
  ### 10 sites ###
  load(paste0("C:/R/SENSOR/EJ_SAFE/outputs/INLA_Anopheles/Models_An/spde_ss5_it",iteration[j],"_n10.RData"))
  
  ## Posterior mean total number of mosquito bites per location
  Post.Res <- inla.posterior.sample(100, res)
  
  # Simulated draws from the negative binomial distribution
  samples.df <- sapply(Post.Res, pred.nb)
  
  # Summary values for posterior samples per row (per cell)
  sum.post <- apply(samples.df, 1, summary)
  sum.post <- as.data.frame(sum.post)
  sum.post <- t(sum.post)
  
  # SD values for posterior samples per row (per cell) 
  sd.post <- apply(samples.df, 1, sd)
  sd.post <- as.data.frame(round(sd.post, digits=4))
  
  # Join 
  output.post <- cbind(sum.post, sd.post)
  rm(sd.post); rm(sum.post)
  
  # Cell-wise (row-wise) RMSE of posteriors (per raster grid across space/time)
  for(i in 1:nrow(samples.df)) {
    
    predicted <- samples.df[i,]
    actual <- pmos[i]
    d <- (predicted - actual)^2
    rsd <- sqrt(sd(d))
    rmse <- rmse(actual, predicted)
    rmse.rsd <- cbind(rmse, rsd)
    #   rmse.rsd <- round(rmse.rsd, digits=4)
    
    if(exists("rmse.samp")){
      rmse.samp <- rbind(rmse.samp, rmse.rsd)
      
    }
    
    # If merged dataset does not exist, create dataset
    if(!exists("rmse.samp")){
      rmse.samp <- rmse.rsd
    }
    rm(rmse); rm(rsd); rm(rmse.rsd)
    
  } 
  
  output <- cbind(output.post, rmse.samp)
  rm(output.post);  rm(rmse.samp)
  
  
  # rbind vertical join 
  output$type <- "geospat"
  output$num <- 10
  output$iter <- iteration[j]
  
  
  if(exists("output.all")){
    output.all <- rbind(output.all, output)
    
  }
  
  # If merged dataset does not exist, create dataset
  if(!exists("output.all")){
    output.all <- output
  }
  
  rm(output)
  
  
  # Col-wise overall RMSE per posterior sample 100
  for(i in 1:ncol(samples.df)) {
    
    predicted <- samples.df[,i]
    actual <- pmos
    rmse <- rmse(actual, predicted)
    
    if(exists("rmse.mod")){
      rmse.mod <- rbind(rmse.mod, rmse)
      
    }
    
    # If merged dataset does not exist, create dataset
    if(!exists("rmse.mod")){
      rmse.mod <- rmse
    }
    rm(rmse)
  } 
  
  rmse.df <- as.data.frame(rmse.mod)
  rmse.df <- t(rmse.df)
  rmse.df <- as.data.frame(rmse.df)
  rmse.df$iter <- iteration[j]
  rmse.df$type <- "geospat"
  rmse.df$num <- 10
  
  if(exists("rmse.all")){
    rmse.all <- rbind(rmse.all, rmse.df)
    
  }
  
  # If merged dataset does not exist, create dataset
  if(!exists("rmse.all")){
    rmse.all <- rmse.df
  }
  rm(rmse.df); rm(rmse.mod)
  print(iteration[j])
  
} 

write.csv(output.all,"C:/R/SENSOR/EJ_SAFE/outputs/INLA_Anopheles/Posteriors_CSV/An_output_posterior_ss5_n10.csv")
write.csv(rmse.all,"C:/R/SENSOR/EJ_SAFE/outputs/INLA_Anopheles/Posteriors_CSV/An_rmse_posterior_ss5_n10.csv")
rm(output.all); rm(rmse.all) 


for(j in 1:length(iteration)){ 
  
  
  ### 15 sites ###
  load(paste0("C:/R/SENSOR/EJ_SAFE/outputs/INLA_Anopheles/Models_An/spde_ss5_it",iteration[j],"_n15.RData"))
  
  ## Posterior mean total number of mosquito bites per location
  Post.Res <- inla.posterior.sample(100, res)
  
  # Simulated draws from the negative binomial distribution
  samples.df <- sapply(Post.Res, pred.nb)
  
  # Summary values for posterior samples per row (per cell)
  sum.post <- apply(samples.df, 1, summary)
  sum.post <- as.data.frame(sum.post)
  sum.post <- t(sum.post)
  
  # SD values for posterior samples per row (per cell) 
  sd.post <- apply(samples.df, 1, sd)
  sd.post <- as.data.frame(round(sd.post, digits=4))
  
  # Join 
  output.post <- cbind(sum.post, sd.post)
  rm(sd.post); rm(sum.post)
  
  # Cell-wise (row-wise) RMSE of posteriors (per raster grid across space/time)
  for(i in 1:nrow(samples.df)) {
    
    predicted <- samples.df[i,]
    actual <- pmos[i]
    d <- (predicted - actual)^2
    rsd <- sqrt(sd(d))
    rmse <- rmse(actual, predicted)
    rmse.rsd <- cbind(rmse, rsd)
    #   rmse.rsd <- round(rmse.rsd, digits=4)
    
    if(exists("rmse.samp")){
      rmse.samp <- rbind(rmse.samp, rmse.rsd)
      
    }
    
    # If merged dataset does not exist, create dataset
    if(!exists("rmse.samp")){
      rmse.samp <- rmse.rsd
    }
    rm(rmse); rm(rsd); rm(rmse.rsd)
    
  } 
  
  output <- cbind(output.post, rmse.samp)
  rm(output.post);  rm(rmse.samp)
  
  
  # rbind vertical join 
  output$type <- "geospat"
  output$num <- 15
  output$iter <- iteration[j]
  
  
  if(exists("output.all")){
    output.all <- rbind(output.all, output)
    
  }
  
  # If merged dataset does not exist, create dataset
  if(!exists("output.all")){
    output.all <- output
  }
  
  rm(output)
  
  
  # Col-wise overall RMSE per posterior sample 100
  for(i in 1:ncol(samples.df)) {
    
    predicted <- samples.df[,i]
    actual <- pmos
    rmse <- rmse(actual, predicted)
    
    if(exists("rmse.mod")){
      rmse.mod <- rbind(rmse.mod, rmse)
      
    }
    
    # If merged dataset does not exist, create dataset
    if(!exists("rmse.mod")){
      rmse.mod <- rmse
    }
    rm(rmse)
  } 
  
  rmse.df <- as.data.frame(rmse.mod)
  rmse.df <- t(rmse.df)
  rmse.df <- as.data.frame(rmse.df)
  rmse.df$iter <- iteration[j]
  rmse.df$type <- "geospat"
  rmse.df$num <- 15
  
  if(exists("rmse.all")){
    rmse.all <- rbind(rmse.all, rmse.df)
    
  }
  
  # If merged dataset does not exist, create dataset
  if(!exists("rmse.all")){
    rmse.all <- rmse.df
  }
  rm(rmse.df); rm(rmse.mod)
  print(iteration[j])
  
} 

write.csv(output.all,"C:/R/SENSOR/EJ_SAFE/outputs/INLA_Anopheles/Posteriors_CSV/An_output_posterior_ss5_n15.csv")
write.csv(rmse.all,"C:/R/SENSOR/EJ_SAFE/outputs/INLA_Anopheles/Posteriors_CSV/An_rmse_posterior_ss5_n15.csv")
rm(output.all); rm(rmse.all) 




##########################################################################
################################ AEDES ###################################
##########################################################################

# create raster and make shapefile to get coordinates back
roi  <- readOGR(dsn="C:/Users/Kimberly Fornace/OneDrive - University of Glasgow/SENSOR/R/SAFE_ento/", layer="ROI_extent")
roi <- spTransform(roi, CRS("+proj=utm +zone=50 +datum=WGS84 +units=km +no_defs +ellps=WGS84 +towgs84=0,0,0"))
rr <- raster(ext = extent(roi), resolution = 1)
values(rr) <- 1 # set all vaues to one
df_rr <- as.data.frame(rr, xy=T) # get df from raster for spatial points
xmin <- min(df_rr$x)
ymin <- min(df_rr$y)

#simulated
pred <- read.csv("C:/Users/Kimberly Fornace/OneDrive - University of Glasgow/SENSOR/R/SAFE_ento/simulated_aedes_500m.csv") # sim_mos
pred$hab <- relevel(factor(pred$hab), ref = "Primary")
pred$hab <- plyr::revalue(pred$hab, c("Oil"="oil", "Secondary"="secondary", "Primary"="primary", "Built"="built", "Plantation"="plantation"))
pred$xcoo <- pred$xcoo + xmin
pred$ycoo <- pred$ycoo + ymin

# vector of actual ('real' simulated) values 
pmos <- pred$mos

## Simulate 100 draws from the posterior negative binomial distribution
# Kim code 2023
pred.nb <- function(l, nsample=nrow(pred)){
  lp <- l$latent[1:nsample] ## extract the linear predictor in log scale
  lambda  <-  exp(lp)   ## apply the link function
  y <- rnbinom(nsample,mu=lambda, size=l$hyperpar[1]) ## posterior sample of non-zero inflated data
}

setwd("C:/R/SENSOR/EJ_SAFE/outputs/INLA_Aedes/Posteriors_Ae")

# ss1 - random ###########################################################

iteration <- 1:100

for(j in 1:length(iteration)){ 
  
  ### 5 sites ###
  load(paste0("C:/R/SENSOR/EJ_SAFE/outputs/INLA_Aedes/Models_Ae/spde_ss1_it",iteration[j],"_n5.RData"))
  
  ## Posterior mean total number of mosquito bites per location
  Post.Res <- inla.posterior.sample(100, res)
  
  # Simulated draws from the negative binomial distribution
  samples.df <- sapply(Post.Res, pred.nb)
  
  # Summary values for posterior samples per row (per cell)
  sum.post <- apply(samples.df, 1, summary)
  sum.post <- as.data.frame(sum.post)
  sum.post <- t(sum.post)
  
  # SD values for posterior samples per row (per cell) 
  sd.post <- apply(samples.df, 1, sd)
  sd.post <- as.data.frame(round(sd.post, digits=4))
  
  # Join 
  output.post <- cbind(sum.post, sd.post)
  rm(sd.post); rm(sum.post)
  
  # Cell-wise (row-wise) RMSE of posteriors (per raster grid across space/time)
  for(i in 1:nrow(samples.df)) {
    
    predicted <- samples.df[i,]
    actual <- pmos[i]
    d <- (predicted - actual)^2
    rsd <- sqrt(sd(d))
    rmse <- rmse(actual, predicted)
    rmse.rsd <- cbind(rmse, rsd)
    #   rmse.rsd <- round(rmse.rsd, digits=4)
    
    if(exists("rmse.samp")){
      rmse.samp <- rbind(rmse.samp, rmse.rsd)
      
    }
    
    # If merged dataset does not exist, create dataset
    if(!exists("rmse.samp")){
      rmse.samp <- rmse.rsd
    }
    rm(rmse); rm(rsd); rm(rmse.rsd)
    
  } 
  
  output <- cbind(output.post, rmse.samp)
  rm(output.post);  rm(rmse.samp)
  
  
  # rbind vertical join 
  output$type <- "random"
  output$num <- 5
  output$iter <- iteration[j]
  
  
  if(exists("output.all")){
    output.all <- rbind(output.all, output)
    
  }
  
  # If merged dataset does not exist, create dataset
  if(!exists("output.all")){
    output.all <- output
  }
  
  rm(output)
  
  
  # Col-wise overall RMSE per posterior sample 100
  for(i in 1:ncol(samples.df)) {
    
    predicted <- samples.df[,i]
    actual <- pmos
    rmse <- rmse(actual, predicted)
    
    if(exists("rmse.mod")){
      rmse.mod <- rbind(rmse.mod, rmse)
      
    }
    
    # If merged dataset does not exist, create dataset
    if(!exists("rmse.mod")){
      rmse.mod <- rmse
    }
    rm(rmse)
  } 
  
  rmse.df <- as.data.frame(rmse.mod)
  rmse.df <- t(rmse.df)
  rmse.df <- as.data.frame(rmse.df)
  rmse.df$iter <- iteration[j]
  rmse.df$type <- "random"
  rmse.df$num <- 5
  
  if(exists("rmse.all")){
    rmse.all <- rbind(rmse.all, rmse.df)
    
  }
  
  # If merged dataset does not exist, create dataset
  if(!exists("rmse.all")){
    rmse.all <- rmse.df
  }
  rm(rmse.df); rm(rmse.mod)
  print(iteration[j])
  
} 

write.csv(output.all,"C:/R/SENSOR/EJ_SAFE/outputs/INLA_Aedes/Posteriors_CSV/Ae_output_posterior_ss1_n5.csv")
write.csv(rmse.all,"C:/R/SENSOR/EJ_SAFE/outputs/INLA_Aedes/Posteriors_CSV/Ae_rmse_posterior_ss1_n5.csv")
rm(output.all); rm(rmse.all)

for(j in 1:length(iteration)){ 
  
  ### 10 sites ###
  load(paste0("C:/R/SENSOR/EJ_SAFE/outputs/INLA_Aedes/Models_Ae/spde_ss1_it",iteration[j],"_n10.RData"))
  
  ## Posterior mean total number of mosquito bites per location
  Post.Res <- inla.posterior.sample(100, res)
  
  # Simulated draws from the negative binomial distribution
  samples.df <- sapply(Post.Res, pred.nb)
  
  # Summary values for posterior samples per row (per cell)
  sum.post <- apply(samples.df, 1, summary)
  sum.post <- as.data.frame(sum.post)
  sum.post <- t(sum.post)
  
  # SD values for posterior samples per row (per cell) 
  sd.post <- apply(samples.df, 1, sd)
  sd.post <- as.data.frame(round(sd.post, digits=4))
  
  # Join 
  output.post <- cbind(sum.post, sd.post)
  rm(sd.post); rm(sum.post)
  
  # Cell-wise (row-wise) RMSE of posteriors (per raster grid across space/time)
  for(i in 1:nrow(samples.df)) {
    
    predicted <- samples.df[i,]
    actual <- pmos[i]
    d <- (predicted - actual)^2
    rsd <- sqrt(sd(d))
    rmse <- rmse(actual, predicted)
    rmse.rsd <- cbind(rmse, rsd)
    #   rmse.rsd <- round(rmse.rsd, digits=4)
    
    if(exists("rmse.samp")){
      rmse.samp <- rbind(rmse.samp, rmse.rsd)
      
    }
    
    # If merged dataset does not exist, create dataset
    if(!exists("rmse.samp")){
      rmse.samp <- rmse.rsd
    }
    rm(rmse); rm(rsd); rm(rmse.rsd)
    
  } 
  
  output <- cbind(output.post, rmse.samp)
  rm(output.post);  rm(rmse.samp)
  
  
  # rbind vertical join 
  output$type <- "random"
  output$num <- 10
  output$iter <- iteration[j]
  
  
  if(exists("output.all")){
    output.all <- rbind(output.all, output)
    
  }
  
  # If merged dataset does not exist, create dataset
  if(!exists("output.all")){
    output.all <- output
  }
  
  rm(output)
  
  
  # Col-wise overall RMSE per posterior sample 100
  for(i in 1:ncol(samples.df)) {
    
    predicted <- samples.df[,i]
    actual <- pmos
    rmse <- rmse(actual, predicted)
    
    if(exists("rmse.mod")){
      rmse.mod <- rbind(rmse.mod, rmse)
      
    }
    
    # If merged dataset does not exist, create dataset
    if(!exists("rmse.mod")){
      rmse.mod <- rmse
    }
    rm(rmse)
  } 
  
  rmse.df <- as.data.frame(rmse.mod)
  rmse.df <- t(rmse.df)
  rmse.df <- as.data.frame(rmse.df)
  rmse.df$iter <- iteration[j]
  rmse.df$type <- "random"
  rmse.df$num <- 10
  
  if(exists("rmse.all")){
    rmse.all <- rbind(rmse.all, rmse.df)
    
  }
  
  # If merged dataset does not exist, create dataset
  if(!exists("rmse.all")){
    rmse.all <- rmse.df
  }
  rm(rmse.df); rm(rmse.mod)
  print(iteration[j])
  
} 

write.csv(output.all,"C:/R/SENSOR/EJ_SAFE/outputs/INLA_Aedes/Posteriors_CSV/Ae_output_posterior_ss1_n10.csv")
write.csv(rmse.all,"C:/R/SENSOR/EJ_SAFE/outputs/INLA_Aedes/Posteriors_CSV/Ae_rmse_posterior_ss1_n10.csv")
rm(output.all); rm(rmse.all)



for(j in 1:length(iteration)){ 
  
  ### 15 sites ###
  load(paste0("C:/R/SENSOR/EJ_SAFE/outputs/INLA_Aedes/Models_Ae/spde_ss1_it",iteration[j],"_n15.RData"))
  
  ## Posterior mean total number of mosquito bites per location
  Post.Res <- inla.posterior.sample(100, res)
  
  # Simulated draws from the negative binomial distribution
  samples.df <- sapply(Post.Res, pred.nb)
  
  # Summary values for posterior samples per row (per cell)
  sum.post <- apply(samples.df, 1, summary)
  sum.post <- as.data.frame(sum.post)
  sum.post <- t(sum.post)
  
  # SD values for posterior samples per row (per cell) 
  sd.post <- apply(samples.df, 1, sd)
  sd.post <- as.data.frame(round(sd.post, digits=4))
  
  # Join 
  output.post <- cbind(sum.post, sd.post)
  rm(sd.post); rm(sum.post)
  
  # Cell-wise (row-wise) RMSE of posteriors (per raster grid across space/time)
  for(i in 1:nrow(samples.df)) {
    
    predicted <- samples.df[i,]
    actual <- pmos[i]
    d <- (predicted - actual)^2
    rsd <- sqrt(sd(d))
    rmse <- rmse(actual, predicted)
    rmse.rsd <- cbind(rmse, rsd)
    #   rmse.rsd <- round(rmse.rsd, digits=4)
    
    if(exists("rmse.samp")){
      rmse.samp <- rbind(rmse.samp, rmse.rsd)
      
    }
    
    # If merged dataset does not exist, create dataset
    if(!exists("rmse.samp")){
      rmse.samp <- rmse.rsd
    }
    rm(rmse); rm(rsd); rm(rmse.rsd)
    
  } 
  
  output <- cbind(output.post, rmse.samp)
  rm(output.post);  rm(rmse.samp)
  
  
  # rbind vertical join 
  output$type <- "random"
  output$num <- 15
  output$iter <- iteration[j]
  
  
  if(exists("output.all")){
    output.all <- rbind(output.all, output)
    
  }
  
  # If merged dataset does not exist, create dataset
  if(!exists("output.all")){
    output.all <- output
  }
  
  rm(output)
  
  
  # Col-wise overall RMSE per posterior sample 100
  for(i in 1:ncol(samples.df)) {
    
    predicted <- samples.df[,i]
    actual <- pmos
    rmse <- rmse(actual, predicted)
    
    if(exists("rmse.mod")){
      rmse.mod <- rbind(rmse.mod, rmse)
      
    }
    
    # If merged dataset does not exist, create dataset
    if(!exists("rmse.mod")){
      rmse.mod <- rmse
    }
    rm(rmse)
  } 
  
  rmse.df <- as.data.frame(rmse.mod)
  rmse.df <- t(rmse.df)
  rmse.df <- as.data.frame(rmse.df)
  rmse.df$iter <- iteration[j]
  rmse.df$type <- "random"
  rmse.df$num <- 15
  
  if(exists("rmse.all")){
    rmse.all <- rbind(rmse.all, rmse.df)
    
  }
  
  # If merged dataset does not exist, create dataset
  if(!exists("rmse.all")){
    rmse.all <- rmse.df
  }
  rm(rmse.df); rm(rmse.mod)
  print(iteration[j])
  
} 

write.csv(output.all,"C:/R/SENSOR/EJ_SAFE/outputs/INLA_Aedes/Posteriors_CSV/Ae_output_posterior_ss1_n15.csv")
write.csv(rmse.all,"C:/R/SENSOR/EJ_SAFE/outputs/INLA_Aedes/Posteriors_CSV/Ae_rmse_posterior_ss1_n15.csv")
rm(output.all); rm(rmse.all)




# ss2 - fixed  ############################

iteration <- c(1:71,73:100)

for(j in 1:length(iteration)){ 
  
  ### 5 sites ###
  load(paste0("C:/R/SENSOR/EJ_SAFE/outputs/INLA_Aedes/Models_Ae/spde_ss2_it",iteration[j],"_n5.RData"))
  
  ## Posterior mean total number of mosquito bites per location
  Post.Res <- inla.posterior.sample(100, res)
  
  # Simulated draws from the negative binomial distribution
  samples.df <- sapply(Post.Res, pred.nb)
  
  # Summary values for posterior samples per row (per cell)
  sum.post <- apply(samples.df, 1, summary)
  sum.post <- as.data.frame(sum.post)
  sum.post <- t(sum.post)
  
  # SD values for posterior samples per row (per cell) 
  sd.post <- apply(samples.df, 1, sd)
  sd.post <- as.data.frame(round(sd.post, digits=4))
  
  # Join 
  output.post <- cbind(sum.post, sd.post)
  rm(sd.post); rm(sum.post)
  
  # Cell-wise (row-wise) RMSE of posteriors (per raster grid across space/time)
  for(i in 1:nrow(samples.df)) {
    
    predicted <- samples.df[i,]
    actual <- pmos[i]
    d <- (predicted - actual)^2
    rsd <- sqrt(sd(d))
    rmse <- rmse(actual, predicted)
    rmse.rsd <- cbind(rmse, rsd)
    #   rmse.rsd <- round(rmse.rsd, digits=4)
    
    if(exists("rmse.samp")){
      rmse.samp <- rbind(rmse.samp, rmse.rsd)
      
    }
    
    # If merged dataset does not exist, create dataset
    if(!exists("rmse.samp")){
      rmse.samp <- rmse.rsd
    }
    rm(rmse); rm(rsd); rm(rmse.rsd)
    
  } 
  
  output <- cbind(output.post, rmse.samp)
  rm(output.post);  rm(rmse.samp)
  
  
  # rbind vertical join 
  output$type <- "fixed"
  output$num <- 5
  output$iter <- iteration[j]
  
  
  if(exists("output.all")){
    output.all <- rbind(output.all, output)
    
  }
  
  # If merged dataset does not exist, create dataset
  if(!exists("output.all")){
    output.all <- output
  }
  
  rm(output)
  
  
  # Col-wise overall RMSE per posterior sample 100
  for(i in 1:ncol(samples.df)) {
    
    predicted <- samples.df[,i]
    actual <- pmos
    rmse <- rmse(actual, predicted)
    
    if(exists("rmse.mod")){
      rmse.mod <- rbind(rmse.mod, rmse)
      
    }
    
    # If merged dataset does not exist, create dataset
    if(!exists("rmse.mod")){
      rmse.mod <- rmse
    }
    rm(rmse)
  } 
  
  rmse.df <- as.data.frame(rmse.mod)
  rmse.df <- t(rmse.df)
  rmse.df <- as.data.frame(rmse.df)
  rmse.df$iter <- iteration[j]
  rmse.df$type <- "fixed"
  rmse.df$num <- 5
  
  if(exists("rmse.all")){
    rmse.all <- rbind(rmse.all, rmse.df)
    
  }
  
  # If merged dataset does not exist, create dataset
  if(!exists("rmse.all")){
    rmse.all <- rmse.df
  }
  rm(rmse.df); rm(rmse.mod)
  print(iteration[j])
  
}

write.csv(output.all,"C:/R/SENSOR/EJ_SAFE/outputs/INLA_Aedes/Posteriors_CSV/Ae_output_posterior_ss2_n5.csv")
write.csv(rmse.all,"C:/R/SENSOR/EJ_SAFE/outputs/INLA_Aedes/Posteriors_CSV/Ae_rmse_posterior_ss2_n5.csv")
rm(output.all); rm(rmse.all)



iteration <- 1:100

for(j in 1:length(iteration)){ 
  
  ### 10 sites ###
  load(paste0("C:/R/SENSOR/EJ_SAFE/outputs/INLA_Aedes/Models_Ae/spde_ss2_it",iteration[j],"_n10.RData"))
  
  ## Posterior mean total number of mosquito bites per location
  Post.Res <- inla.posterior.sample(100, res)
  
  # Simulated draws from the negative binomial distribution
  samples.df <- sapply(Post.Res, pred.nb)
  
  # Summary values for posterior samples per row (per cell)
  sum.post <- apply(samples.df, 1, summary)
  sum.post <- as.data.frame(sum.post)
  sum.post <- t(sum.post)
  
  # SD values for posterior samples per row (per cell) 
  sd.post <- apply(samples.df, 1, sd)
  sd.post <- as.data.frame(round(sd.post, digits=4))
  
  # Join 
  output.post <- cbind(sum.post, sd.post)
  rm(sd.post); rm(sum.post)
  
  # Cell-wise (row-wise) RMSE of posteriors (per raster grid across space/time)
  for(i in 1:nrow(samples.df)) {
    
    predicted <- samples.df[i,]
    actual <- pmos[i]
    d <- (predicted - actual)^2
    rsd <- sqrt(sd(d))
    rmse <- rmse(actual, predicted)
    rmse.rsd <- cbind(rmse, rsd)
    #   rmse.rsd <- round(rmse.rsd, digits=4)
    
    if(exists("rmse.samp")){
      rmse.samp <- rbind(rmse.samp, rmse.rsd)
      
    }
    
    # If merged dataset does not exist, create dataset
    if(!exists("rmse.samp")){
      rmse.samp <- rmse.rsd
    }
    rm(rmse); rm(rsd); rm(rmse.rsd)
    
  } 
  
  output <- cbind(output.post, rmse.samp)
  rm(output.post);  rm(rmse.samp)
  
  
  # rbind vertical join 
  output$type <- "fixed"
  output$num <- 10
  output$iter <- iteration[j]
  
  
  if(exists("output.all")){
    output.all <- rbind(output.all, output)
    
  }
  
  # If merged dataset does not exist, create dataset
  if(!exists("output.all")){
    output.all <- output
  }
  
  rm(output)
  
  
  # Col-wise overall RMSE per posterior sample 100
  for(i in 1:ncol(samples.df)) {
    
    predicted <- samples.df[,i]
    actual <- pmos
    rmse <- rmse(actual, predicted)
    
    if(exists("rmse.mod")){
      rmse.mod <- rbind(rmse.mod, rmse)
      
    }
    
    # If merged dataset does not exist, create dataset
    if(!exists("rmse.mod")){
      rmse.mod <- rmse
    }
    rm(rmse)
  } 
  
  rmse.df <- as.data.frame(rmse.mod)
  rmse.df <- t(rmse.df)
  rmse.df <- as.data.frame(rmse.df)
  rmse.df$iter <- iteration[j]
  rmse.df$type <- "fixed"
  rmse.df$num <- 10
  
  if(exists("rmse.all")){
    rmse.all <- rbind(rmse.all, rmse.df)
    
  }
  
  # If merged dataset does not exist, create dataset
  if(!exists("rmse.all")){
    rmse.all <- rmse.df
  }
  rm(rmse.df); rm(rmse.mod)
  print(iteration[j])
  
} 

write.csv(output.all,"C:/R/SENSOR/EJ_SAFE/outputs/INLA_Aedes/Posteriors_CSV/Ae_output_posterior_ss2_n10.csv")
write.csv(rmse.all,"C:/R/SENSOR/EJ_SAFE/outputs/INLA_Aedes/Posteriors_CSV/Ae_rmse_posterior_ss2_n10.csv")
rm(output.all); rm(rmse.all)



for(j in 1:length(iteration)){ 
  
  ### 15 sites ###
  load(paste0("C:/R/SENSOR/EJ_SAFE/outputs/INLA_Aedes/Models_Ae/spde_ss2_it",iteration[j],"_n15.RData"))
  
  ## Posterior mean total number of mosquito bites per location
  Post.Res <- inla.posterior.sample(100, res)
  
  # Simulated draws from the negative binomial distribution
  samples.df <- sapply(Post.Res, pred.nb)
  
  # Summary values for posterior samples per row (per cell)
  sum.post <- apply(samples.df, 1, summary)
  sum.post <- as.data.frame(sum.post)
  sum.post <- t(sum.post)
  
  # SD values for posterior samples per row (per cell) 
  sd.post <- apply(samples.df, 1, sd)
  sd.post <- as.data.frame(round(sd.post, digits=4))
  
  # Join 
  output.post <- cbind(sum.post, sd.post)
  rm(sd.post); rm(sum.post)
  
  # Cell-wise (row-wise) RMSE of posteriors (per raster grid across space/time)
  for(i in 1:nrow(samples.df)) {
    
    predicted <- samples.df[i,]
    actual <- pmos[i]
    d <- (predicted - actual)^2
    rsd <- sqrt(sd(d))
    rmse <- rmse(actual, predicted)
    rmse.rsd <- cbind(rmse, rsd)
    #   rmse.rsd <- round(rmse.rsd, digits=4)
    
    if(exists("rmse.samp")){
      rmse.samp <- rbind(rmse.samp, rmse.rsd)
      
    }
    
    # If merged dataset does not exist, create dataset
    if(!exists("rmse.samp")){
      rmse.samp <- rmse.rsd
    }
    rm(rmse); rm(rsd); rm(rmse.rsd)
    
  } 
  
  output <- cbind(output.post, rmse.samp)
  rm(output.post);  rm(rmse.samp)
  
  
  # rbind vertical join 
  output$type <- "fixed"
  output$num <- 15
  output$iter <- iteration[j]
  
  
  if(exists("output.all")){
    output.all <- rbind(output.all, output)
    
  }
  
  # If merged dataset does not exist, create dataset
  if(!exists("output.all")){
    output.all <- output
  }
  
  rm(output)
  
  
  # Col-wise overall RMSE per posterior sample 100
  for(i in 1:ncol(samples.df)) {
    
    predicted <- samples.df[,i]
    actual <- pmos
    rmse <- rmse(actual, predicted)
    
    if(exists("rmse.mod")){
      rmse.mod <- rbind(rmse.mod, rmse)
      
    }
    
    # If merged dataset does not exist, create dataset
    if(!exists("rmse.mod")){
      rmse.mod <- rmse
    }
    rm(rmse)
  } 
  
  rmse.df <- as.data.frame(rmse.mod)
  rmse.df <- t(rmse.df)
  rmse.df <- as.data.frame(rmse.df)
  rmse.df$iter <- iteration[j]
  rmse.df$type <- "fixed"
  rmse.df$num <- 15
  
  if(exists("rmse.all")){
    rmse.all <- rbind(rmse.all, rmse.df)
    
  }
  
  # If merged dataset does not exist, create dataset
  if(!exists("rmse.all")){
    rmse.all <- rmse.df
  }
  rm(rmse.df); rm(rmse.mod)
  print(iteration[j])
  
} 

write.csv(output.all,"C:/R/SENSOR/EJ_SAFE/outputs/INLA_Aedes/Posteriors_CSV/Ae_output_posterior_ss2_n15.csv")
write.csv(rmse.all,"C:/R/SENSOR/EJ_SAFE/outputs/INLA_Aedes/Posteriors_CSV/Ae_rmse_posterior_ss2_n15.csv")
rm(output.all); rm(rmse.all)



# ss3 - hab #######################  

iteration <- 1:100

for(j in 1:length(iteration)){ 
  
  ### 5 sites ###
  load(paste0("C:/R/SENSOR/EJ_SAFE/outputs/INLA_Aedes/Models_Ae/spde_ss3_it",iteration[j],"_n5.RData"))
  
  ## Posterior mean total number of mosquito bites per location
  Post.Res <- inla.posterior.sample(100, res)
  
  # Simulated draws from the negative binomial distribution
  samples.df <- sapply(Post.Res, pred.nb)
  
  # Summary values for posterior samples per row (per cell)
  sum.post <- apply(samples.df, 1, summary)
  sum.post <- as.data.frame(sum.post)
  sum.post <- t(sum.post)
  
  # SD values for posterior samples per row (per cell) 
  sd.post <- apply(samples.df, 1, sd)
  sd.post <- as.data.frame(round(sd.post, digits=4))
  
  # Join 
  output.post <- cbind(sum.post, sd.post)
  rm(sd.post); rm(sum.post)
  
  # Cell-wise (row-wise) RMSE of posteriors (per raster grid across space/time)
  for(i in 1:nrow(samples.df)) {
    
    predicted <- samples.df[i,]
    actual <- pmos[i]
    d <- (predicted - actual)^2
    rsd <- sqrt(sd(d))
    rmse <- rmse(actual, predicted)
    rmse.rsd <- cbind(rmse, rsd)
    #   rmse.rsd <- round(rmse.rsd, digits=4)
    
    if(exists("rmse.samp")){
      rmse.samp <- rbind(rmse.samp, rmse.rsd)
      
    }
    
    # If merged dataset does not exist, create dataset
    if(!exists("rmse.samp")){
      rmse.samp <- rmse.rsd
    }
    rm(rmse); rm(rsd); rm(rmse.rsd)
    
  } 
  
  output <- cbind(output.post, rmse.samp)
  rm(output.post);  rm(rmse.samp)
  
  
  # rbind vertical join 
  output$type <- "hab_strat"
  output$num <- 5
  output$iter <- iteration[j]
  
  
  if(exists("output.all")){
    output.all <- rbind(output.all, output)
    
  }
  
  # If merged dataset does not exist, create dataset
  if(!exists("output.all")){
    output.all <- output
  }
  
  rm(output)
  
  
  # Col-wise overall RMSE per posterior sample 100
  for(i in 1:ncol(samples.df)) {
    
    predicted <- samples.df[,i]
    actual <- pmos
    rmse <- rmse(actual, predicted)
    
    if(exists("rmse.mod")){
      rmse.mod <- rbind(rmse.mod, rmse)
      
    }
    
    # If merged dataset does not exist, create dataset
    if(!exists("rmse.mod")){
      rmse.mod <- rmse
    }
    rm(rmse)
  } 
  
  rmse.df <- as.data.frame(rmse.mod)
  rmse.df <- t(rmse.df)
  rmse.df <- as.data.frame(rmse.df)
  rmse.df$iter <- iteration[j]
  rmse.df$type <- "hab_strat"
  rmse.df$num <- 5
  
  if(exists("rmse.all")){
    rmse.all <- rbind(rmse.all, rmse.df)
    
  }
  
  # If merged dataset does not exist, create dataset
  if(!exists("rmse.all")){
    rmse.all <- rmse.df
  }
  rm(rmse.df); rm(rmse.mod)
  print(iteration[j])
  
} 

write.csv(output.all,"C:/R/SENSOR/EJ_SAFE/outputs/INLA_Aedes/Posteriors_CSV/Ae_output_posterior_ss3_n5.csv")
write.csv(rmse.all,"C:/R/SENSOR/EJ_SAFE/outputs/INLA_Aedes/Posteriors_CSV/Ae_rmse_posterior_ss3_n5.csv")
rm(output.all); rm(rmse.all)



iteration <- c(1:68,70:100)

for(j in 1:length(iteration)){ 
  
  ### 10 sites ###
  load(paste0("C:/R/SENSOR/EJ_SAFE/outputs/INLA_Aedes/Models_Ae/spde_ss3_it",iteration[j],"_n10.RData"))
  
  ## Posterior mean total number of mosquito bites per location
  Post.Res <- inla.posterior.sample(100, res)
  
  # Simulated draws from the negative binomial distribution
  samples.df <- sapply(Post.Res, pred.nb)
  
  # Summary values for posterior samples per row (per cell)
  sum.post <- apply(samples.df, 1, summary)
  sum.post <- as.data.frame(sum.post)
  sum.post <- t(sum.post)
  
  # SD values for posterior samples per row (per cell) 
  sd.post <- apply(samples.df, 1, sd)
  sd.post <- as.data.frame(round(sd.post, digits=4))
  
  # Join 
  output.post <- cbind(sum.post, sd.post)
  rm(sd.post); rm(sum.post)
  
  # Cell-wise (row-wise) RMSE of posteriors (per raster grid across space/time)
  for(i in 1:nrow(samples.df)) {
    
    predicted <- samples.df[i,]
    actual <- pmos[i]
    d <- (predicted - actual)^2
    rsd <- sqrt(sd(d))
    rmse <- rmse(actual, predicted)
    rmse.rsd <- cbind(rmse, rsd)
    #   rmse.rsd <- round(rmse.rsd, digits=4)
    
    if(exists("rmse.samp")){
      rmse.samp <- rbind(rmse.samp, rmse.rsd)
      
    }
    
    # If merged dataset does not exist, create dataset
    if(!exists("rmse.samp")){
      rmse.samp <- rmse.rsd
    }
    rm(rmse); rm(rsd); rm(rmse.rsd)
    
  } 
  
  output <- cbind(output.post, rmse.samp)
  rm(output.post);  rm(rmse.samp)
  
  
  # rbind vertical join 
  output$type <- "hab_strat"
  output$num <- 10
  output$iter <- iteration[j]
  
  
  if(exists("output.all")){
    output.all <- rbind(output.all, output)
    
  }
  
  # If merged dataset does not exist, create dataset
  if(!exists("output.all")){
    output.all <- output
  }
  
  rm(output)
  
  
  # Col-wise overall RMSE per posterior sample 100
  for(i in 1:ncol(samples.df)) {
    
    predicted <- samples.df[,i]
    actual <- pmos
    rmse <- rmse(actual, predicted)
    
    if(exists("rmse.mod")){
      rmse.mod <- rbind(rmse.mod, rmse)
      
    }
    
    # If merged dataset does not exist, create dataset
    if(!exists("rmse.mod")){
      rmse.mod <- rmse
    }
    rm(rmse)
  } 
  
  rmse.df <- as.data.frame(rmse.mod)
  rmse.df <- t(rmse.df)
  rmse.df <- as.data.frame(rmse.df)
  rmse.df$iter <- iteration[j]
  rmse.df$type <- "hab_strat"
  rmse.df$num <- 10
  
  if(exists("rmse.all")){
    rmse.all <- rbind(rmse.all, rmse.df)
    
  }
  
  # If merged dataset does not exist, create dataset
  if(!exists("rmse.all")){
    rmse.all <- rmse.df
  }
  rm(rmse.df); rm(rmse.mod)
  print(iteration[j])
  
} 

write.csv(output.all,"C:/R/SENSOR/EJ_SAFE/outputs/INLA_Aedes/Posteriors_CSV/Ae_output_posterior_ss3_n10.csv")
write.csv(rmse.all,"C:/R/SENSOR/EJ_SAFE/outputs/INLA_Aedes/Posteriors_CSV/Ae_rmse_posterior_ss3_n10.csv")
rm(output.all); rm(rmse.all)



iteration <- 1:100

for(j in 1:length(iteration)){ 
  
  ### 15 sites ###
  load(paste0("C:/R/SENSOR/EJ_SAFE/outputs/INLA_Aedes/Models_Ae/spde_ss3_it",iteration[j],"_n15.RData"))
  
  ## Posterior mean total number of mosquito bites per location
  Post.Res <- inla.posterior.sample(100, res)
  
  # Simulated draws from the negative binomial distribution
  samples.df <- sapply(Post.Res, pred.nb)
  
  # Summary values for posterior samples per row (per cell)
  sum.post <- apply(samples.df, 1, summary)
  sum.post <- as.data.frame(sum.post)
  sum.post <- t(sum.post)
  
  # SD values for posterior samples per row (per cell) 
  sd.post <- apply(samples.df, 1, sd)
  sd.post <- as.data.frame(round(sd.post, digits=4))
  
  # Join 
  output.post <- cbind(sum.post, sd.post)
  rm(sd.post); rm(sum.post)
  
  # Cell-wise (row-wise) RMSE of posteriors (per raster grid across space/time)
  for(i in 1:nrow(samples.df)) {
    
    predicted <- samples.df[i,]
    actual <- pmos[i]
    d <- (predicted - actual)^2
    rsd <- sqrt(sd(d))
    rmse <- rmse(actual, predicted)
    rmse.rsd <- cbind(rmse, rsd)
    #   rmse.rsd <- round(rmse.rsd, digits=4)
    
    if(exists("rmse.samp")){
      rmse.samp <- rbind(rmse.samp, rmse.rsd)
      
    }
    
    # If merged dataset does not exist, create dataset
    if(!exists("rmse.samp")){
      rmse.samp <- rmse.rsd
    }
    rm(rmse); rm(rsd); rm(rmse.rsd)
    
  } 
  
  output <- cbind(output.post, rmse.samp)
  rm(output.post);  rm(rmse.samp)
  
  
  # rbind vertical join 
  output$type <- "hab_strat"
  output$num <- 15
  output$iter <- iteration[j]
  
  
  if(exists("output.all")){
    output.all <- rbind(output.all, output)
    
  }
  
  # If merged dataset does not exist, create dataset
  if(!exists("output.all")){
    output.all <- output
  }
  
  rm(output)
  
  
  # Col-wise overall RMSE per posterior sample 100
  for(i in 1:ncol(samples.df)) {
    
    predicted <- samples.df[,i]
    actual <- pmos
    rmse <- rmse(actual, predicted)
    
    if(exists("rmse.mod")){
      rmse.mod <- rbind(rmse.mod, rmse)
      
    }
    
    # If merged dataset does not exist, create dataset
    if(!exists("rmse.mod")){
      rmse.mod <- rmse
    }
    rm(rmse)
  } 
  
  rmse.df <- as.data.frame(rmse.mod)
  rmse.df <- t(rmse.df)
  rmse.df <- as.data.frame(rmse.df)
  rmse.df$iter <- iteration[j]
  rmse.df$type <- "hab_strat"
  rmse.df$num <- 15
  
  if(exists("rmse.all")){
    rmse.all <- rbind(rmse.all, rmse.df)
    
  }
  
  # If merged dataset does not exist, create dataset
  if(!exists("rmse.all")){
    rmse.all <- rmse.df
  }
  rm(rmse.df); rm(rmse.mod)
  print(iteration[j])
  
} 

write.csv(output.all,"C:/R/SENSOR/EJ_SAFE/outputs/INLA_Aedes/Posteriors_CSV/Ae_output_posterior_ss3_n15.csv")
write.csv(rmse.all,"C:/R/SENSOR/EJ_SAFE/outputs/INLA_Aedes/Posteriors_CSV/Ae_rmse_posterior_ss3_n15.csv")
rm(output.all); rm(rmse.all)



# ss4 - hab fixed #######################  

iteration <- 1:100

for(j in 1:length(iteration)){ 
  
  
  ### 5 sites ###
  load(paste0("C:/R/SENSOR/EJ_SAFE/outputs/INLA_Aedes/Models_Ae/spde_ss4_it",iteration[j],"_n5.RData"))
  
  ## Posterior mean total number of mosquito bites per location
  Post.Res <- inla.posterior.sample(100, res)
  
  # Simulated draws from the negative binomial distribution
  samples.df <- sapply(Post.Res, pred.nb)
  
  # Summary values for posterior samples per row (per cell)
  sum.post <- apply(samples.df, 1, summary)
  sum.post <- as.data.frame(sum.post)
  sum.post <- t(sum.post)
  
  # SD values for posterior samples per row (per cell) 
  sd.post <- apply(samples.df, 1, sd)
  sd.post <- as.data.frame(round(sd.post, digits=4))
  
  # Join 
  output.post <- cbind(sum.post, sd.post)
  rm(sd.post); rm(sum.post)
  
  # Cell-wise (row-wise) RMSE of posteriors (per raster grid across space/time)
  for(i in 1:nrow(samples.df)) {
    
    predicted <- samples.df[i,]
    actual <- pmos[i]
    d <- (predicted - actual)^2
    rsd <- sqrt(sd(d))
    rmse <- rmse(actual, predicted)
    rmse.rsd <- cbind(rmse, rsd)
    #   rmse.rsd <- round(rmse.rsd, digits=4)
    
    if(exists("rmse.samp")){
      rmse.samp <- rbind(rmse.samp, rmse.rsd)
      
    }
    
    # If merged dataset does not exist, create dataset
    if(!exists("rmse.samp")){
      rmse.samp <- rmse.rsd
    }
    rm(rmse); rm(rsd); rm(rmse.rsd)
    
  } 
  
  output <- cbind(output.post, rmse.samp)
  rm(output.post);  rm(rmse.samp)
  
  
  # rbind vertical join 
  output$type <- "hab_fixed"
  output$num <- 5
  output$iter <- iteration[j]
  
  
  if(exists("output.all")){
    output.all <- rbind(output.all, output)
    
  }
  
  # If merged dataset does not exist, create dataset
  if(!exists("output.all")){
    output.all <- output
  }
  
  rm(output)
  
  
  # Col-wise overall RMSE per posterior sample 100
  for(i in 1:ncol(samples.df)) {
    
    predicted <- samples.df[,i]
    actual <- pmos
    rmse <- rmse(actual, predicted)
    
    if(exists("rmse.mod")){
      rmse.mod <- rbind(rmse.mod, rmse)
      
    }
    
    # If merged dataset does not exist, create dataset
    if(!exists("rmse.mod")){
      rmse.mod <- rmse
    }
    rm(rmse)
  } 
  
  rmse.df <- as.data.frame(rmse.mod)
  rmse.df <- t(rmse.df)
  rmse.df <- as.data.frame(rmse.df)
  rmse.df$iter <- iteration[j]
  rmse.df$type <- "hab_fixed"
  rmse.df$num <- 5
  
  if(exists("rmse.all")){
    rmse.all <- rbind(rmse.all, rmse.df)
    
  }
  
  # If merged dataset does not exist, create dataset
  if(!exists("rmse.all")){
    rmse.all <- rmse.df
  }
  rm(rmse.df); rm(rmse.mod)
  print(iteration[j])
  
} 

write.csv(output.all,"C:/R/SENSOR/EJ_SAFE/outputs/INLA_Aedes/Posteriors_CSV/Ae_output_posterior_ss4_n5.csv")
write.csv(rmse.all,"C:/R/SENSOR/EJ_SAFE/outputs/INLA_Aedes/Posteriors_CSV/Ae_rmse_posterior_ss4_n5.csv")
rm(output.all); rm(rmse.all)


for(j in 1:length(iteration)){ 
  
  ### 10 sites ###
  load(paste0("C:/R/SENSOR/EJ_SAFE/outputs/INLA_Aedes/Models_Ae/spde_ss4_it",iteration[j],"_n10.RData"))
  
  ## Posterior mean total number of mosquito bites per location
  Post.Res <- inla.posterior.sample(100, res)
  
  # Simulated draws from the negative binomial distribution
  samples.df <- sapply(Post.Res, pred.nb)
  
  # Summary values for posterior samples per row (per cell)
  sum.post <- apply(samples.df, 1, summary)
  sum.post <- as.data.frame(sum.post)
  sum.post <- t(sum.post)
  
  # SD values for posterior samples per row (per cell) 
  sd.post <- apply(samples.df, 1, sd)
  sd.post <- as.data.frame(round(sd.post, digits=4))
  
  # Join 
  output.post <- cbind(sum.post, sd.post)
  rm(sd.post); rm(sum.post)
  
  # Cell-wise (row-wise) RMSE of posteriors (per raster grid across space/time)
  for(i in 1:nrow(samples.df)) {
    
    predicted <- samples.df[i,]
    actual <- pmos[i]
    d <- (predicted - actual)^2
    rsd <- sqrt(sd(d))
    rmse <- rmse(actual, predicted)
    rmse.rsd <- cbind(rmse, rsd)
    #   rmse.rsd <- round(rmse.rsd, digits=4)
    
    if(exists("rmse.samp")){
      rmse.samp <- rbind(rmse.samp, rmse.rsd)
      
    }
    
    # If merged dataset does not exist, create dataset
    if(!exists("rmse.samp")){
      rmse.samp <- rmse.rsd
    }
    rm(rmse); rm(rsd); rm(rmse.rsd)
    
  } 
  
  output <- cbind(output.post, rmse.samp)
  rm(output.post);  rm(rmse.samp)
  
  
  # rbind vertical join 
  output$type <- "hab_fixed"
  output$num <- 10
  output$iter <- iteration[j]
  
  
  if(exists("output.all")){
    output.all <- rbind(output.all, output)
    
  }
  
  # If merged dataset does not exist, create dataset
  if(!exists("output.all")){
    output.all <- output
  }
  
  rm(output)
  
  
  # Col-wise overall RMSE per posterior sample 100
  for(i in 1:ncol(samples.df)) {
    
    predicted <- samples.df[,i]
    actual <- pmos
    rmse <- rmse(actual, predicted)
    
    if(exists("rmse.mod")){
      rmse.mod <- rbind(rmse.mod, rmse)
      
    }
    
    # If merged dataset does not exist, create dataset
    if(!exists("rmse.mod")){
      rmse.mod <- rmse
    }
    rm(rmse)
  } 
  
  rmse.df <- as.data.frame(rmse.mod)
  rmse.df <- t(rmse.df)
  rmse.df <- as.data.frame(rmse.df)
  rmse.df$iter <- iteration[j]
  rmse.df$type <- "hab_fixed"
  rmse.df$num <- 10
  
  if(exists("rmse.all")){
    rmse.all <- rbind(rmse.all, rmse.df)
    
  }
  
  # If merged dataset does not exist, create dataset
  if(!exists("rmse.all")){
    rmse.all <- rmse.df
  }
  rm(rmse.df); rm(rmse.mod)
  print(iteration[j])
  
} 


write.csv(output.all,"C:/R/SENSOR/EJ_SAFE/outputs/INLA_Aedes/Posteriors_CSV/Ae_output_posterior_ss4_n10.csv")
write.csv(rmse.all,"C:/R/SENSOR/EJ_SAFE/outputs/INLA_Aedes/Posteriors_CSV/Ae_rmse_posterior_ss4_n10.csv")
rm(output.all); rm(rmse.all)


######################################################




# create raster and make shapefile to get coordinates back
roi  <- readOGR(dsn="C:/Users/Kimberly Fornace/OneDrive - University of Glasgow/SENSOR/R/SAFE_ento/", layer="ROI_extent")
roi <- spTransform(roi, CRS("+proj=utm +zone=50 +datum=WGS84 +units=km +no_defs +ellps=WGS84 +towgs84=0,0,0"))
rr <- raster(ext = extent(roi), resolution = 1)
values(rr) <- 1 # set all vaues to one
df_rr <- as.data.frame(rr, xy=T) # get df from raster for spatial points
xmin <- min(df_rr$x)
ymin <- min(df_rr$y)

#simulated
pred <- read.csv("C:/Users/Kimberly Fornace/OneDrive - University of Glasgow/SENSOR/R/SAFE_ento/simulated_anoph_500m.csv") # sim_mos
pred$hab <- relevel(factor(pred$hab), ref = "Primary")
pred$hab <- plyr::revalue(pred$hab, c("Oil"="oil", "Secondary"="secondary", "Primary"="primary", "Built"="built", "Plantation"="plantation"))
pred$xcoo <- pred$xcoo + xmin
pred$ycoo <- pred$ycoo + ymin

# vector of actual ('real' simulated) values 
pmos <- pred$mos

#random
b <- read.csv("C:/Users/Kimberly Fornace/OneDrive - University of Glasgow/Desktop/Mosq Modelling/Output_Sims/Random_Grids_Monthly_100.csv")
b$hab <- relevel(factor(b$hab), ref = "2") # ref = Primary (cat=2)
b$hab <- plyr::revalue(b$hab, c("0"="oil", "1"="secondary", "2"="primary", "3"="plantation", "4"="built"))
b$num <- b$loc
b$xcoo <- b$x + xmin
b$ycoo <- b$y + ymin

#random-fixed
c <- read.csv("C:/Users/Kimberly Fornace/OneDrive - University of Glasgow/Desktop/Mosq Modelling/Output_Sims/Random_Fixed_Monthly_100.csv")
c$hab <- relevel(factor(c$hab), ref = "2") # ref = Primary (cat=2)
c$hab <- plyr::revalue(c$hab, c("0"="oil", "1"="secondary", "2"="primary", "3"="plantation", "4"="built"))
c$num <- c$loc
c$xcoo <- c$x + xmin
c$ycoo <- c$y + ymin

#stratified
d <- read.csv("C:/Users/Kimberly Fornace/OneDrive - University of Glasgow/Desktop/Mosq Modelling/Output_Sims/Stratified_Random_Monthly_100.csv")
d$hab <- relevel(factor(d$hab), ref = "2") # ref = Primary (cat=2)
d$hab <- plyr::revalue(d$hab, c("0"="oil", "1"="secondary", "2"="primary", "3"="plantation", "4"="built"))
d$num <- d$loc
d$xcoo <- d$x + xmin
d$ycoo <- d$y + ymin

#strat-fixed
e <- read.csv("C:/Users/Kimberly Fornace/OneDrive - University of Glasgow/Desktop/Mosq Modelling/Output_Sims/Fixed_Stratified_Monthly_100.csv")
e$hab <- relevel(factor(e$hab), ref = "2") # ref = Primary (cat=2)
e$hab <- plyr::revalue(e$hab, c("0"="oil", "1"="secondary", "2"="primary", "3"="plantation", "4"="built"))
e$num <- e$loc
e$xcoo <- e$x + xmin
e$ycoo <- e$y + ymin

#geospatial
f <- read.csv("C:/Users/Kimberly Fornace/OneDrive - University of Glasgow/Desktop/Mosq Modelling/Output_Sims/Geospatial_Fixed_100.csv")
f$hab <- relevel(factor(f$hab), ref = "2") # ref = Primary (cat=2)
f$hab <- plyr::revalue(f$hab, c("0"="oil", "1"="secondary", "2"="primary", "3"="plantation", "4"="built"))
f$num <- f$loc
f$xcoo <- f$x + xmin
f$ycoo <- f$y + ymin


# 5 sites ------------------

iteration <- 1:100

# Null  model
dat <- e
dat <- subset(dat, dat$num == 5)

for(j in 1:length(iteration)){ 
  
  # Subset iteration 
  it <- subset(dat, dat$iter == iteration[j])
  
  # inla
  coo <- cbind(it$xcoo, it$ycoo)
  mesh <- inla.mesh.2d(loc=coo, max.edge = c(2, 4), cutoff = 1)
  spde.a <- inla.spde2.matern(mesh = mesh, alpha=2, constr = T)
  indexs <- inla.spde.make.index("s", spde.a$n.spde)
  A <- inla.spde.make.A(mesh = mesh, loc = coo) 
  coop <- cbind(pred$xcoo, pred$ycoo)
  Ap <- inla.spde.make.A(mesh = mesh, loc = coop)
  stke <- inla.stack(tag="est", data=list(mos = it$mos),A = list(1, A), effects = list(data.frame(hab = it$hab, time=it$mon), s = indexs))
  stkp <- inla.stack(tag = "pred", data = list(mos = NA),A = list(1, Ap), effects = list(data.frame( hab = pred$hab, time=pred$time), s = indexs)) # b0 = rep(1, nrow(sim_mos_pred)),
  stkfull <- inla.stack(stke, stkp)
  res <- inla(mos ~ 1 + f(s, model = spde.a)  + f(time, model = "ar1"), 
              family = "nbinomial", control.family = list(link = "log"), data = inla.stack.data(stkfull), 
              control.predictor = list(compute = TRUE, link=1, A = inla.stack.A(stkfull)),
              control.compute = list(config=T, dic = TRUE, waic = TRUE)) # sim_mos
  rm(coo);rm(mesh);rm(indexs);rm(A);rm(coop);rm(Ap);rm(stke);rm(stkp)
  
  # Save coefficients
  inla_coefs_spde_ar <- plyr::rbind.fill(res$summary.fixed, res$summary.hyperpar)
  rownames(inla_coefs_spde_ar) <- c("(Intercept)" , "size for the nbinomial observations (1/overdispersion)", "Theta1 for s", "Theta2 for s", "Precision for time", "Rho for time")
  
  inla_coefs_spde_ar$type <- "hab_fixed"
  inla_coefs_spde_ar$num <- 5
  inla_coefs_spde_ar$iter <- iteration[j]
  inla_coefs_spde_ar$DIC <- res$dic$dic
  inla_coefs_spde_ar$WAIC <- res$waic$waic
  
  if(exists("coefs_df_spde_ar")){
    coefs_df_spde_ar <- rbind(coefs_df_spde_ar, inla_coefs_spde_ar)
    
  }
  
  # If merged dataset does not exist, create dataset
  if(!exists("coefs_df_spde_ar")){
    coefs_df_spde_ar <- inla_coefs_spde_ar
  }
  
  rm(inla_coefs_spde_ar)
  
  ### Calculate spatial range ----------------------------------
  spde.est <- inla.spde2.result(inla = res, name = "s", spde = spde.a, do.transf = TRUE)
  sr <- as.data.frame(inla.zmarginal(spde.est$marginals.range.nominal[[1]]))
  sr$type <- "hab_fixed"
  sr$num <- 5
  sr$iter <- iteration[j]
  
  if(exists("INLA_spatial_range_spde_ar")){
    INLA_spatial_range_spde_ar <- rbind(INLA_spatial_range_spde_ar, sr)
    
  }
  
  # If merged dataset does not exist, create dataset
  if(!exists("INLA_spatial_range_spde_ar")){
    INLA_spatial_range_spde_ar <- sr
  }
  rm(sr)
  
  
  
  print(j)
  
}

# 10 sites ------------------
iteration <- 47:100

# Null  model
dat <- e
dat <- subset(dat, dat$num == 10)

for(j in 1:length(iteration)){ 
  
  # Subset iteration 
  it <- subset(dat, dat$iter == iteration[j])
  
  # inla
  coo <- cbind(it$xcoo, it$ycoo)
  mesh <- inla.mesh.2d(loc=coo, max.edge = c(2, 4), cutoff = 1)
  spde.a <- inla.spde2.matern(mesh = mesh, alpha=2, constr = T)
  indexs <- inla.spde.make.index("s", spde.a$n.spde)
  A <- inla.spde.make.A(mesh = mesh, loc = coo) 
  coop <- cbind(pred$xcoo, pred$ycoo)
  Ap <- inla.spde.make.A(mesh = mesh, loc = coop)
  stke <- inla.stack(tag="est", data=list(mos = it$mos),A = list(1, A), effects = list(data.frame(hab = it$hab, time=it$mon), s = indexs))
  stkp <- inla.stack(tag = "pred", data = list(mos = NA),A = list(1, Ap), effects = list(data.frame( hab = pred$hab, time=pred$time), s = indexs)) # b0 = rep(1, nrow(sim_mos_pred)),
  stkfull <- inla.stack(stke, stkp)
  res <- inla(mos ~ 1 + f(s, model = spde.a)  + f(time, model = "ar1"), 
              family = "nbinomial", control.family = list(link = "log"), data = inla.stack.data(stkfull), 
              control.predictor = list(compute = TRUE, link=1, A = inla.stack.A(stkfull)),
              control.compute = list(config=T, dic = TRUE, waic = TRUE)) # sim_mos
  rm(coo);rm(mesh);rm(indexs);rm(A);rm(coop);rm(Ap);rm(stke);rm(stkp)
  
  # Save coefficients
  inla_coefs_spde_ar <- plyr::rbind.fill(res$summary.fixed, res$summary.hyperpar)
  rownames(inla_coefs_spde_ar) <- c("(Intercept)" , "size for the nbinomial observations (1/overdispersion)", "Theta1 for s", "Theta2 for s", "Precision for time", "Rho for time")
  
  inla_coefs_spde_ar$type <- "hab_fixed"
  inla_coefs_spde_ar$num <- 10
  inla_coefs_spde_ar$iter <- iteration[j]
  inla_coefs_spde_ar$DIC <- res$dic$dic
  inla_coefs_spde_ar$WAIC <- res$waic$waic
  
  if(exists("coefs_df_spde_ar")){
    coefs_df_spde_ar <- rbind(coefs_df_spde_ar, inla_coefs_spde_ar)
    
  }
  
  # If merged dataset does not exist, create dataset
  if(!exists("coefs_df_spde_ar")){
    coefs_df_spde_ar <- inla_coefs_spde_ar
  }
  
  rm(inla_coefs_spde_ar)
  
  ### Calculate spatial range ----------------------------------
  spde.est <- inla.spde2.result(inla = res, name = "s", spde = spde.a, do.transf = TRUE)
  sr <- as.data.frame(inla.zmarginal(spde.est$marginals.range.nominal[[1]]))
  sr$type <- "hab_fixed"
  sr$num <- 10
  sr$iter <- iteration[j]
  
  if(exists("INLA_spatial_range_spde_ar")){
    INLA_spatial_range_spde_ar <- rbind(INLA_spatial_range_spde_ar, sr)
    
  }
  
  # If merged dataset does not exist, create dataset
  if(!exists("INLA_spatial_range_spde_ar")){
    INLA_spatial_range_spde_ar <- sr
  }
  rm(sr)
  
  
  print(j)
  
}


print(paste0("10done"))


# 15 sites ------------------
iteration <- 1:100

# Null  model
dat <- e
dat <- subset(dat, dat$num == 15)

for(j in 1:length(iteration)){ 
  
  # Subset iteration 
  it <- subset(dat, dat$iter == iteration[j])
  
  # inla
  coo <- cbind(it$xcoo, it$ycoo)
  mesh <- inla.mesh.2d(loc=coo, max.edge = c(2, 4), cutoff = 1)
  spde.a <- inla.spde2.matern(mesh = mesh, alpha=2, constr = T)
  indexs <- inla.spde.make.index("s", spde.a$n.spde)
  A <- inla.spde.make.A(mesh = mesh, loc = coo) 
  coop <- cbind(pred$xcoo, pred$ycoo)
  Ap <- inla.spde.make.A(mesh = mesh, loc = coop)
  stke <- inla.stack(tag="est", data=list(mos = it$mos),A = list(1, A), effects = list(data.frame(hab = it$hab, time = it$mon), s = indexs))
  stkp <- inla.stack(tag = "pred", data = list(mos = NA),A = list(1, Ap), effects = list(data.frame( hab = pred$hab, time=pred$time), s = indexs)) # b0 = rep(1, nrow(sim_mos_pred)),
  stkfull <- inla.stack(stke, stkp)
  res <- inla(mos ~ 1 + f(s, model = spde.a)  + f(time, model = "ar1"), 
              family = "nbinomial", control.family = list(link = "log"), data = inla.stack.data(stkfull), 
              control.predictor = list(compute = TRUE, link=1, A = inla.stack.A(stkfull)),
              control.compute = list(config=T, dic = TRUE, waic = TRUE)) # sim_mos
  rm(coo);rm(mesh);rm(indexs);rm(A);rm(coop);rm(Ap);rm(stke);rm(stkp)
  
  # Save coefficients
  inla_coefs_spde_ar <- plyr::rbind.fill(res$summary.fixed, res$summary.hyperpar)
  rownames(inla_coefs_spde_ar) <- c("(Intercept)" , "size for the nbinomial observations (1/overdispersion)", "Theta1 for s", "Theta2 for s", "Precision for time", "Rho for time")
  
  inla_coefs_spde_ar$type <- "hab_fixed"
  inla_coefs_spde_ar$num <- 15
  inla_coefs_spde_ar$iter <- iteration[j]
  inla_coefs_spde_ar$DIC <- res$dic$dic
  inla_coefs_spde_ar$WAIC <- res$waic$waic
  
  if(exists("coefs_df_spde_ar")){
    coefs_df_spde_ar <- rbind(coefs_df_spde_ar, inla_coefs_spde_ar)
    
  }
  
  # If merged dataset does not exist, create dataset
  if(!exists("coefs_df_spde_ar")){
    coefs_df_spde_ar <- inla_coefs_spde_ar
  }
  
  rm(inla_coefs_spde_ar)
  
  ### Calculate spatial range ----------------------------------
  spde.est <- inla.spde2.result(inla = res, name = "s", spde = spde.a, do.transf = TRUE)
  sr <- as.data.frame(inla.zmarginal(spde.est$marginals.range.nominal[[1]]))
  sr$type <- "hab_fixed"
  sr$num <- 15
  sr$iter <- iteration[j]
  
  if(exists("INLA_spatial_range_spde_ar")){
    INLA_spatial_range_spde_ar <- rbind(INLA_spatial_range_spde_ar, sr)
    
  }
  
  # If merged dataset does not exist, create dataset
  if(!exists("INLA_spatial_range_spde_ar")){
    INLA_spatial_range_spde_ar <- sr
  }
  rm(sr)
  
  print(iteration[j])
  
}
### 

write.csv(coefs_df_spde_ar,"C:/R/SENSOR/EJ_SAFE/outputs/INLA_Anopheles/NoFixed_An/An_coefs_nofix_ss4.csv")
write.csv(INLA_spatial_range_spde_ar,"C:/R/SENSOR/EJ_SAFE/outputs/INLA_Anopheles/NoFixed_An/An_nofix_spatial_range_ss4.csv")
rm(coefs_df_spde_ar)
rm(INLA_spatial_range_spde_ar)

print(paste0("15done"))   
