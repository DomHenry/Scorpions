library(raster)
library(boot)
library(coda)
library(tidyverse)
library(sp)
library(rgdal)
library(here)
library(lattice)
library(rasterVis)

# Read in basic JAGS model workspace --------------------------------------
load(here("data output","JAGSout_B_P1.RData"))
mod <- outB

# Read in all pentad covs, scale to create predictor ----------------------
pencovs_all_sc <- pencovs_all %>% 
  mutate_at(vars(map:ndvi),funs(scale)) %>% 
  mutate_at(vars(map:ndvi),as.numeric)

# Create array ------------------------------------------------------------
(n_samp <- nrow(mod))                      # 27000 is too many!
select_samp <- sort(sample(1:n_samp, 100)) # Chose random sample of 100 posterior draws
(n_samp <- length(select_samp))    
n_all_pen <- nrow(pencovs_all_sc)

## Create posterior predictive distribution for Z for all pentads
str(zPen <- array(NA, dim = c(n_all_pen,n_spp, n_samp)) ) 

# Acccess node names and extract posterior samples ------------------------
att <- attributes(mod)$dimnames[2][[1]][]
coeffs <- c("lpsi","beta1","beta2","beta3","beta4")
coeffOut <- list()

for(i in 1:length(coeffs)){
  ref <- which(grepl(coeffs[i],att) & 
        !grepl(c("sd"),att) & !grepl(c("mu"),att)) # Exclude sd and mu
  coeffOut[[i]] <- mod[,ref]
}
names(coeffOut) <- c("LPSI","BETA1","BETA2","BETA3","BETA4")

# Generate occupancy predictions for 100 samples --------------------------
for(i in 1:n_all_pen){             
  print.vec <- seq(0,1906, by = 100)
  if(i %in% print.vec){cat(paste("\nPentad ", i))}

  for(u in 1:length(select_samp)){
    psi <- plogis(coeffOut[["LPSI"]][select_samp[u],] + 
                    coeffOut[["BETA1"]][select_samp[u],] * pencovs_all_sc$ndvi[i] + 
                    coeffOut[["BETA2"]][select_samp[u],] * pencovs_all_sc$map_ctn[i] +
                    coeffOut[["BETA3"]][select_samp[u],] * pencovs_all_sc$elev[i] + 
                    coeffOut[["BETA4"]][select_samp[u],] * pencovs_all_sc$elev_range[i])
    zPen[i,,u] <- rbinom(n_spp, 1, psi)
    
  }
}

# Compute posterior distribution of SR - collapse z array -----------------
SR <- apply(zPen, c(1,3), sum, na.rm= TRUE)   # posterior distribution
pmSR <- apply(SR, 1, mean, na.rm = TRUE)      # posterior mean
sdSR <- apply(SR, 1, sd, na.rm = TRUE)        # posterior standard deviation
pencovs_all$spec_rich_mean <- pmSR;
pencovs_all$spec_rich_sd <- sdSR

# Join covs, species richness and SD with spatial df ----------------------
pentads_sp <- readOGR(here("data input","BioGapsAllPens.shp"))
towns <- readOGR(here("data input","karoo_towns.shp"))
pentads_sp <- merge(pentads_sp, pencovs_all, by.x = "PENTADE",by.y="pentad")
head(pentads_sp)

ras_ref <- c("spec_rich_mean","spec_rich_sd")
main_lab <- c("Scorpion species richness","Scorpion species richness - SD")

# Plot species richness maps ----------------------------------------------
col_mean <- rev(terrain.colors(255)) 
col_sd = rev(heat.colors(n = 8))
col_list <- list(col_mean,col_sd)

pdf(here("data output","Scorpions_SR_raster.pdf"), width = 16, height = 9)
par(mfrow = c(1,1))

for(i in 1:length(ras_ref)){
  ras <- raster(ncol = 150, nrow = 150)
  extent(ras) <- extent(pentads_sp)+0.5 #Add 0.5 buffer
  res(ras) <- c(0.0833,0.0833)
  ras <- raster::rasterize(pentads_sp,ras,ras_ref[i])
  plot(ras, main = main_lab[i], col = col_list[[i]])
  plot(towns, add = T, pch = 20,cex = 0.7)
  maptools::pointLabel(coordinates(towns),labels=towns$field_2, cex = 0.7)
}

dev.off()

# Notes -------------------------------------------------------------------

## Take the elements of the columns of the covs.all and subtract the mean from covs data set
#trans <- sweep(covs.all[,-1], 2, d.mean,"-") 
## Take the elements of the columns of the trans1 and divide by the sd from the covs data set
#covs.pred <- sweep(trans, 2, d.sd,"/")
## Do it in one line
covs.pred <- sweep(sweep(covs.all[,-1], 2, d.mean,"-"), 2, d.sd,"/") %>% 
  mutate(PENTAD = covs.all$PENTAD) %>%  as_tibble()
covs.pred
## See these link to further understand sweep() 
## https://stackoverflow.com/questions/3444889/how-to-use-the-r-function-sweep
## https://magesblog.com/post/2012-04-14-sweeping-through-data-in-r/

