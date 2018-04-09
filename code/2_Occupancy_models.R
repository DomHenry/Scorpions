library(here)

# Import workspace --------------------------------------------------------
load(here("data output","Arrays.RData"))

# Prepare data ------------------------------------------------------------
dim(Y)
Y[which(is.na(Y))] <- 0
Y[1:10,1:5,1:2]

# Compile data, inits and parameters for JAGS model -----------------------
str(jags_data <- list(Y = Y, n_spp = dim(Y)[3], n_pen = dim(Y)[1], rep_vec = rep_df$n_transects,
                      ndvi = pencovs_sc$ndvi, map_ctn = pencovs_sc$map_ctn, 
                      elev = pencovs_sc$elev, elev_range = pencovs_sc$elev_range,
                      airtemp = airtemp_arr_sc, observer = observer_arr))

zst <- array(1, dim = c(dim(Y)[1],dim(Y)[3])) 
dim(zst)
zst[is.na(zst)] <- 1
inits <- function() list(z = zst)
str(list(z = zst))

# Choose parameters and settings ------------------------------------------
p1 <- c("mu.psi","mu.p","sd.lpsi","sd.lp","alpha1","alpha2",
        "mu.beta1","mu.beta2","mu.beta3","mu.beta4",
        "sd.beta1", "sd.beta2","sd.beta3","sd.beta4",
        "Nsite","Nocc.fs",
        "beta1","beta2","beta3","beta4",
        "lpsi","lp")


p2 <- c("beta1","beta2","beta3")
p3 <- c("z")
p4 <- c("lpsi")
p5 <- c("lp")

## MCMC settings
ni <- 50000
nt <- 2
nb <- 10000
nc <- 3

paramchoose <- p1
pname <- "P1"

# Run JAGS model components -----------------------------------------------
library(jagsUI)
source(here("code","0_JAGS model.R"))

mod.out.F <- jags(jags_data, inits, paramchoose, "OccMod1.txt", n.chains = nc, n.thin = nt,
                n.iter = ni, n.burnin = nb, parallel = T)

print(mod.out.F, dig = 3)

df <- as.data.frame(mod.out.F$summary)
write.csv(df,here("data output","scorpions_modoutF.csv"))
save.image(here("data output","JAGSout_F_P1.RData"))

# Run full JAGS basic model -----------------------------------------------
library(jagsUI)
source(here("code","0_JAGS model.R"))
paramchoose <- c("mu.psi","mu.p","sd.lpsi","sd.lp","alpha1","alpha2",
                 "mu.beta1","mu.beta2","mu.beta3","mu.beta4",
                 "sd.beta1", "sd.beta2","sd.beta3","sd.beta4",
                 "Nsite","Nocc.fs",
                 "beta1","beta2","beta3","beta4",
                 "lpsi","lp")

mod.out.B <- jags.basic(jags_data, inits, paramchoose, "OccMod1.txt", n.chains = nc, n.thin = nt,
                n.iter = ni, n.burnin = nb, parallel = TRUE)

library(coda)
outB <- as.matrix(mod.out.B) # Put output from 3 chains into a matrix
rm(mod.out.B)
save.image(here("data output","JAGSout_B_P1.RData"))

