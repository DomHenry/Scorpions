library(boot)
library(coda)
library(tidyverse)
library(jagsUI)
library(here)
library(gridExtra)
library(ggpubr)
## Plotting of parameters, predictions and diagnostics of bird TTD models

# Import JAGS model workspaces --------------------------------------------
(wslist <- dir(here("data output"),pattern = "F_"))
load(here("data output",wslist[1]))

# Load functions ----------------------------------------------------------
source(here("code","0_Analysis functions.R"))

# Create plotting dataframe -----------------------------------------------
plotdf <- as_tibble(mod.out.F$summary[,1:7]) %>% 
  mutate(param = row.names(mod.out.F$summary)) %>% 
  select(param,mean,`2.5%`,`97.5%`) %>% 
  rename(lower = `2.5%`, upper = `97.5%`)

plotdf;plotdf$param

# Psi & p table -----------------------------------------------------------
pp_means<- plotdf %>% 
  filter(param %in% c("mu.psi","mu.p"))

pp_sds<- plotdf %>%
  filter(param %in% c("sd.lpsi","sd.lp")) %>% 
  mutate_at(vars(mean:upper),funs(inv.logit)) %>% 
  mutate(param = str_replace(param, "l",""))

pp <- bind_rows(pp_means,pp_sds)
pp


# Open PDF ----------------------------------------------------------------
pdf(here("data output","Scorpion_JAGS_model_plots.pdf"), width = 16, height = 9)

# Plot :: Species richness ------------------------------------------------
NsiteR <- plotdf %>% filter(str_detect(param,"Nsite")) %>% 
  mutate(param = as.factor(pen_ref))

ggplot(NsiteR, aes(x = fct_reorder(param,mean), y = mean, ymin = lower, ymax = upper)) +
  geom_pointrange() + theme_bw() + scale_y_continuous(limits =  c(0,max(NsiteR$upper)+2))+
  coord_flip() + geom_hline(yintercept = 0, linetype = "dotted") +
  ylab("Species richness") + xlab("")+
  theme(axis.text.y=element_text(size=8),
        axis.text.x=element_text(size=10),
        axis.title=element_text(size=17,face="bold"))

# Plot :: Detection coefficients ------------------------------------------

# alpha1 ~ dnorm(0, 0.1)  # Air temperature  
# alpha2 ~ dnorm(0, 0.1)  # Effect of observer

alphaR <- plotdf %>% filter(str_detect(param,"alpha")) %>% 
  mutate(param = c("Air temp","Observer"))

ggplot(alphaR, aes(x = fct_reorder(param,mean), y = mean, ymin = lower, ymax = upper)) +
  geom_pointrange() + theme_bw() + scale_y_continuous(limits = c(-1.1,0.6))+
  coord_flip() + geom_hline(yintercept = 0, linetype = "dotted") +
  ylab("Alpha coefficients - detection") + xlab("")+
  theme(axis.text.x=element_text(size=18),
        axis.text.y=element_text(size=18),
        axis.title=element_text(size=18,face="bold"))

# Plot :: Detection covariate predictions ---------------------------------
plot_coeff_curves(airtemp_arr,c(15,30,0.1),"alpha1",NA,"mu.p","Air temperature","Detection probability",c(0,0.8))

# Plot :: Species detection probability -----------------------------------
p.spR <- plotdf %>% filter(str_detect(param, "lp\\[")) %>% 
  mutate(param = spp_ref) %>% 
  mutate_at(vars(mean:upper),funs(inv.logit))

ggplot(p.spR, aes(x = fct_reorder(param,mean), y = mean, ymin = lower, ymax = upper)) +
  geom_pointrange() + theme_bw() + scale_y_continuous(limits = c(0,1))+
  coord_flip() + geom_hline(yintercept = 0, linetype = "dotted") +
  ylab("Mean detection probability") + xlab("")+
  theme(axis.text.y=element_text(size=6),
        axis.text.x=element_text(size=16),
        axis.title=element_text(size=16,face="bold"))

# Plot :: Occupancy coefficients ------------------------------------------

# beta1[k] ~ dnorm(mu.beta1, tau.beta1)                # ndvi
# beta2[k] ~ dnorm(mu.beta2, tau.beta2)                # map_ctn
# beta3[k] ~ dnorm(mu.beta3, tau.beta3)                # elev   
# beta4[k] ~ dnorm(mu.beta4, tau.beta4)                # tri   

betaR <- plotdf %>% filter(str_detect(param, "mu.beta")) %>% 
  mutate(param = as.factor(c("NDVI","RainConc", "Elev","tri")))

ggplot(betaR, aes(x = fct_reorder(param,mean), y = mean, ymin = lower, ymax = upper)) +
  geom_pointrange() + theme_bw() + scale_y_continuous(limits = c(-1,1.5))+
  coord_flip() + geom_hline(yintercept = 0, linetype = "dotted") +
  ylab("beta coefficients - detection") + xlab("")+
  theme(axis.text.x=element_text(size=18),
        axis.text.y=element_text(size=18),
        axis.title=element_text(size=18,face="bold"))


# Plot :: Occupancy predictions -------------------------------------------
covs <- list(pencovs[["ndvi"]],pencovs[["map_ctn"]],pencovs[["elev"]],pencovs[["tri_med"]])
ranges <- list(c(0,0.5,0.001),c(30,50,1),c(300,1700,5),c(10,700,5))
coeffs1 <- list("mu.beta1","mu.beta2","mu.beta3","mu.beta4")
coeffs2 <- NA
intname <- "mu.psi"
xlabs <- list("NDVI","Rainfall concentration","Elevation","TRI")
ylabs <- list(rep("Occupancy probability",4))
ylim <- list(c(0,1))

# AMAZING! 
par(mfrow = c(2,2))

plots_out <- pmap(list(cov = covs, covRange = ranges,
          coeffname1 = coeffs1,coeffname2 = coeffs2,
          intname = intname,xlab = xlabs,
          ylab = ylabs, ylim = ylim),
     plot_coeff_curves)

plots_out

## Display plots
grid.arrange(grobs = plots_out, ncol = 2) 

## Alternatively 
ggarrange(plots_out[[1]],
          plots_out[[2]] + rremove("y.text") + rremove("ylab"),
          plots_out[[3]],
          plots_out[[4]] + rremove("y.text") + rremove("ylab"),
          ncol = 2, nrow = 2)

# Plot :: Species occupancy probability -----------------------------------
psiR <- plotdf %>% filter(str_detect(param, "lpsi")) %>% 
  filter(!str_detect(param,"sd|mu")) %>% mutate(param = spp_ref) %>% 
  mutate_at(vars(mean:upper),inv.logit)

ggplot(psiR, aes(x = fct_reorder(param,mean), y = mean, ymin = lower, ymax = upper)) +
  geom_pointrange() + theme_bw() + scale_y_continuous(limits = c(0,1))+
  coord_flip() + geom_hline(yintercept = 0, linetype = "dotted") +
  ylab("Occupancy probability for each species") + xlab("")+
  theme(axis.text.y=element_text(size=7),
        axis.text.x=element_text(size=16),
        axis.title=element_text(size=16,face="bold"))

# Plot :: Number of occupied sites  ---------------------------------------
occ.fsR <- plotdf %>% filter(str_detect(param, "occ.fs")) %>% 
  mutate(param = spp_ref) 

ggplot(occ.fsR, aes(x = fct_reorder(param,mean), y = mean, ymin = lower, ymax = upper)) +
  geom_pointrange() + theme_bw() + scale_y_continuous(limits = c(0,30))+
  coord_flip() + geom_hline(yintercept = 0, linetype = "dotted") +
  ylab("Est no of sites occupied by each species") + xlab("")+
  theme(axis.text.y=element_text(size=7),
        axis.text.x=element_text(size=16),
        axis.title=element_text(size=16,face="bold"))

# Plot:: Random effects ---------------------------------------------------

coeffs <- list("beta1","beta2","beta3","beta4")
xlims <- list(c(-5,5),c(-2,2),c(-2,2),c(-2,2))
mainlabs <- list("NDVI","MAP_CTN","ELEV","TRI")

par(mfrow = c(1,2), cex.lab = 1.3, cex.axis = 1.3)

pwalk(list(coeff_name = coeffs, main_lab = mainlabs, xlims = xlims, n = n_spp),
     random_effect_plots)


# Plot :: Posterior distributions (all parameters) ------------------------
str(mod.out.F$sims.list)
names(mod.out.F$sims.list)
par(mfrow= c(4,3))

for(i in 1:14){
  id <- names(mod.out.F$sims.list)[i]
  hist(mod.out.F$sims.list[[i]][], col = "grey", breaks = 60, main = "", xlab = id,freq = F)
  abline(v = mod.out.F$summary[id,1], col = "blue", lwd = 3) #mean
  abline(v = mod.out.F$summary[id,c(3,7)], col = "red", lwd = 3)   # add 95% CRI
}

for(i in 1:n_pen){
  ref <- paste0("Nsite[",i,"]")
  xlab <- paste0("No. of species: ",dimnames(Y)[[1]][i])
  hist(mod.out.F$sims.list$Nsite[,i], col = "grey", breaks = 60, xlim = c(0,15), main = "", xlab = xlab, freq = F)
  abline(v = mod.out.F$mean$Nsite[i] , col = "blue", lwd = 3) #mean
  abline(v = mod.out.F$summary[ref,c(3,7)], col = "red", lwd = 3)  
  
}

for(i in 1:n_spp){
  ref <- paste0("Nocc.fs[",i,"]")
  xlab <- paste0("No.occ sites: ",dimnames(Y)[[3]][i])
  hist(mod.out.F$sims.list$Nocc.fs[,i], col = "grey", breaks = 60, xlim = c(0, 30), main = "", xlab = xlab, freq = F)
  abline(v = mod.out.F$mean$Nocc.fs[i] , col = "blue", lwd = 3) #mean
  abline(v = mod.out.F$summary[ref,c(3,7)], col = "red", lwd = 3)   
}

for(i in 1:n_spp){
  ref <- paste0("lp[",i,"]")
  xlab <- paste0("Detection prob: ",dimnames(Y)[[3]][i])
  hist(inv.logit(mod.out.F$sims.list$lp[,i]), col = "grey", breaks = 60, xlim = c(0,1), main = "", xlab = xlab, freq = F)
  abline(v = inv.logit(mod.out.F$mean$lp[i]), col = "blue", lwd = 3) #mean
  abline(v = inv.logit(mod.out.F$summary[ref,c(3,7)]), col = "red", lwd = 3)   
}

for(i in 1:n_spp){
  ref <- paste0("lpsi[",i,"]")
  xlab <- paste0("Occupancy prob: ",dimnames(Y)[[3]][i])
  hist(inv.logit(mod.out.F$sims.list$lpsi[,i]), col = "grey", breaks = 60, xlim = c(0,1), main = "", xlab = xlab, freq = F)
  abline(v = inv.logit(mod.out.F$mean$lpsi[i]), col = "blue", lwd = 3) #mean
  abline(v = inv.logit(mod.out.F$summary[ref,c(3,7)]), col = "red", lwd = 3)   
}

for(i in 1:n_spp){
  ref <- paste0("beta1[",i,"]")
  xlab <- paste0("NDVI: ",dimnames(Y)[[3]][i])
  hist(mod.out.F$sims.list$beta1[,i], col = "grey", breaks = 60, xlim = c(-10,10), main = "", xlab = xlab, freq = F)
  abline(v = mod.out.F$mean$beta1[i] , col = "blue", lwd = 3) #mean
  abline(v = mod.out.F$summary[ref,c(3,7)], col = "red", lwd = 3)  
  
}

for(i in 1:n_spp){
  ref <- paste0("beta2[",i,"]")
  xlab <- paste0("MAP_CTN: ",dimnames(Y)[[3]][i])
  hist(mod.out.F$sims.list$beta2[,i], col = "grey", breaks = 60, xlim = c(-5,5), main = "", xlab = xlab, freq = F)
  abline(v = mod.out.F$mean$beta2[i] , col = "blue", lwd = 3) #mean
  abline(v = mod.out.F$summary[ref,c(3,7)], col = "red", lwd = 3)  
  
}

for(i in 1:n_spp){
  ref <- paste0("beta3[",i,"]")
  xlab <- paste0("ELEV: ",dimnames(Y)[[3]][i])
  hist(mod.out.F$sims.list$beta3[,i], col = "grey", breaks = 60, xlim = c(-5,5), main = "", xlab = xlab, freq = F)
  abline(v = mod.out.F$mean$beta3[i] , col = "blue", lwd = 3) #mean
  abline(v = mod.out.F$summary[ref,c(3,7)], col = "red", lwd = 3)  
  
}

for(i in 1:n_spp){
  ref <- paste0("beta4[",i,"]")
  xlab <- paste0("ELEV RANGE: ",dimnames(Y)[[3]][i])
  hist(mod.out.F$sims.list$beta4[,i], col = "grey", breaks = 60, xlim = c(-5,5), main = "", xlab = xlab, freq = F)
  abline(v = mod.out.F$mean$beta4[i] , col = "blue", lwd = 3) #mean
  abline(v = mod.out.F$summary[ref,c(3,7)], col = "red", lwd = 3)  
  
}

hist(mod.out.F$sims.list$deviance, col = "grey", breaks = 60, xlim = c(600,900), main = "", xlab = "Deviance", freq = F)
abline(v = mod.out.F$mean$deviance , col = "blue", lwd = 3) #mean
abline(v = mod.out.F$summary["deviance",c(3,7)], col = "red", lwd = 3)   

# Close PDF ---------------------------------------------------------------
dev.off()

# Plot:: Traceplots -------------------------------------------------------
sum <- as.data.frame(mod.out.F$summary)
rhats <- round(sum$Rhat, digits = 3)
neff <- sum$n.eff
v <- seq(0, nrow(sum), by = 9)
dir.create(here("data output/traceplots"))

for(p in 1:length(v)){
  
  name <- paste0(here("data output/traceplots"),"/traceplot",p,".png")
  png(name, width = 3200, height = 1600)
  par(mfrow=c(3,3))
  
  for(i in 1:9){
    ref <- v[p]+i
    traceplot(mod.out.F$samples[,ref], smooth = FALSE, type = "l", xlab = rownames(sum)[ref], ylab = "",xaxt = "n",
              main = paste0("Rhat ",rhats[ref]," | Neff ",neff[ref]), cex.lab = 3,cex.axis =3, cex.main = 3)
  }
  dev.off()
}
dev.off()

system('CMD /C "ECHO The tracepplots have finsihed writing to disk && PAUSE"', 
       invisible=FALSE, wait=FALSE)
