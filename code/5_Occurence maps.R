library(raster)
library(boot)
library(coda)
library(tidyverse)
library(sp)
library(rgdal)
library(here)
library(lattice)
library(rasterVis)
## Occurence maps for each species

# Load basic JAGS workspace -----------------------------------------------
load(here("data output","JAGSout_B_P1.RData"))
mod <- outB

# Scale pentad covariates for study area ----------------------------------
pencovs_all_sc <- pencovs_all %>% 
  mutate_at(vars(map:ndvi),funs(scale)) %>% 
  mutate_at(vars(map:ndvi),as.numeric)

# Create posterior predictive distribution for Z for all pentads ----------
(n_samp <- nrow(mod))                      
select_samp <- sort(sample(1:n_samp, 100)) # Chose random sample of 100 posterior draws
(n_samp <- length(select_samp))    
n_pen_all <- nrow(pencovs_all_sc)
str(psi_pen <- array(NA, dim = c(n_pen_all,n_samp,n_spp))) 

# Acccess node names and extract posterior samples ------------------------
att <- attributes(mod)$dimnames[2][[1]][]
coeffs <- c("lpsi","beta1","beta2","beta3","beta4")
coeff_out <- list()

for(i in 1:length(coeffs)){
  ref <- which(grepl(coeffs[i],att) & 
                 !grepl(c("sd"),att) & !grepl(c("mu"),att)) # Exclude sd and mu
  coeff_out[[i]] <- mod[,ref]
}
names(coeff_out) <- c("LPSI","BETA1","BETA2","BETA3","BETA4")
coeff_out$LPSI[1:5,1:5]
coeff_out$BETA4[1:5,1:5]

# Occupancy prediction ----------------------------------------------------
# beta1[k] ~ dnorm(mu.beta1, tau.beta1)                # ndvi
# beta2[k] ~ dnorm(mu.beta2, tau.beta2)                # map_ctn
# beta3[k] ~ dnorm(mu.beta3, tau.beta3)                # elev   
# beta4[k] ~ dnorm(mu.beta4, tau.beta4)                # elev_range 

for (s in 1:n_spp){
  for(i in 1:n_pen_all){
    cat(paste("\nSpecies", s, "\nPentad", i, "\n"))
    for(u in 1:length(select_samp)){
      psi <- plogis(coeff_out[["LPSI"]][select_samp[u],s] + 
                      coeff_out[["BETA1"]][select_samp[u],s] * pencovs_all_sc$ndvi[i]+ 
                      coeff_out[["BETA2"]][select_samp[u],s] * pencovs_all_sc$map_ctn[i]+
                      coeff_out[["BETA3"]][select_samp[u],s] * pencovs_all_sc$elev[i] +
                      coeff_out[["BETA4"]][select_samp[u],s] * pencovs_all_sc$elev_range[i])
      psi_pen[i,u,s] <- psi
    }
  }
}
pm_psi <- apply(psi_pen, c(1,3), mean, na.rm = TRUE)      # posterior mean
sd_psi <- apply(psi_pen, c(1,3), sd, na.rm = TRUE)   

# Save workspace ----------------------------------------------------------

## First remove large objects 
save(list = c(ls()[!ls() %in% c("mod","coeff_out")]),file = here("data output","Occurence maps.RData"))
## Import "Occurence maps.RData" if need to alter the appearence of maps

# Plot species mean occupancy ---------------------------------------------
pentads_sp <- readOGR(here("data input","BioGapsAllPens.shp"))
towns <- readOGR(here("data input","karoo_towns.shp"))
df_psi <- as.data.frame(pm_psi)
covs_psi <- cbind(pencovs_all,df_psi)
pentads <- subset(pentads_sp, pentads_sp$PENTADE %in% covs_psi$pentad)
pentads <- merge(pentads, covs_psi, by.x = "PENTADE",by.y="pentad")
head(pentads)

brks <- seq(0, 1, by=0.1) 
nb <- length(brks)-1 
cols <- rev(terrain.colors(nb))
sppref <- paste0("V",c(1:n_spp))

pdf(here("data output","Scorpions_SppOccupancy_mean.pdf"), width = 16, height = 9)
par(mfrow = c(1,1))
for(i in seq_along(sppref)){
  PSIras <- raster(ncol = 150, nrow = 150)
  extent(PSIras) <- extent(pentads)+0.5 #Add 0.5 buffer
  res(PSIras) <- c(0.0833,0.0833)
  PSIras <- rasterize(pentads,PSIras,sppref[i])
  plot(PSIras,breaks=brks, col=cols, lab.breaks=brks,
       zlim=c(0,1), main = spp_ref[i])
  plot(towns, add = T, pch = 20,cex = 0.7)
  maptools::pointLabel(coordinates(towns),labels=towns$field_2, cex = 0.7)
}
dev.off()

# Plot species SD occupancy -----------------------------------------------
pentads_sp <- readOGR(here("data input","BioGapsAllPens.shp"))
df_psi <- as.data.frame(sd_psi)
covs_psi <- cbind(pencovs_all,df_psi)
pentads <- subset(pentads_sp, pentads_sp$PENTADE %in% covs_psi$pentad)
pentads <- merge(pentads, covs_psi, by.x = "PENTADE",by.y="pentad")
head(pentads)

brks <- seq(0, 0.5, by=0.05) 
nb <- length(brks)-1 
cols <- rev(heat.colors(nb))
sppref <- paste0("V",c(1:n_spp))

pdf(here("data output","Scorpions_SppOccupancy_standev.pdf"), width = 16, height = 9)

par(mfrow = c(1,1))
for(i in seq_along(sppref)){
  PSIras <- raster(ncol = 150, nrow = 150)
  extent(PSIras) <- extent(pentads)+0.5 #Add 0.5 buffer
  res(PSIras) <- c(0.0833,0.0833)
  PSIras <- rasterize(pentads,PSIras,sppref[i])
  plot(PSIras,breaks=brks, #lab.breaks=brks,
       col = cols, zlim=c(0,0.5), main = spp_ref[i]) # Change this zlim in other plots
  plot(towns, add = T, pch = 20,cex = 0.7)
  maptools::pointLabel(coordinates(towns),labels=towns$field_2, cex = 0.7)
  
}
dev.off()



# End ---------------------------------------------------------------------


