## Functions for the analysis of the Scorpion data set

# Plot TTD detection covariates  ------------------------------------------
plot_coeff_curves <- function(cov,covRange,coeffname1,coeffname2,intname,xlab,ylim) {
  
  vec <- which(!is.na(cov))
  vec <- cov[vec]
  
  mean <- attr(scale(vec), "scaled:center") # mean 
  sd <- attr(scale(vec), "scaled:scale")  # standard deviation
  orig.pred <- seq(from = covRange[1], to = covRange[2], by = covRange[3])# Predictor (x-values) which are unscaled 
  sc.pred <- (orig.pred - mean)/sd 
  
  predictions <- array(dim = c(length(sc.pred), length(mod.out.F$sims.list[[coeffname1]])))
  dim(predictions)
  
  if (is.na(coeffname2)) {
    
    for(i in 1:length(sc.pred)){
      
      predictions[i,] <- plogis(logit(mod.out.F$sims.list[[intname]]) +
                                  mod.out.F$sims.list[[coeffname1]] * sc.pred[i])
      
    }  
  } else {
    
    for(i in 1:length(sc.pred)){
      predictions[i,] <- plogis(logit(mod.out.F$sims.list[[intname]]) +
                                  mod.out.F$sims.list[[coeffname1]] * sc.pred[i] + 
                                  mod.out.F$sims.list[[coeffname2]] * sc.pred[i]^2)
    }  
    
  }
  
  LPB <-  apply(predictions, 1, quantile, probs = 0.025) # Lower bound
  UPB <-  apply(predictions, 1, quantile, probs = 0.975) # Upper bound
  y <- apply(predictions, 1, mean) # Upper bound
  ylim <- ylim
  
  plot(orig.pred, y, xlab = xlab, ylab = "Detection prob.", type = "n",
       ylim = ylim, col = "red", lwd = 3, frame = F,
       cex.lab = 1.6) #Plot with original labels
  polygon(c(orig.pred, rev(orig.pred)), c(UPB, rev(LPB)), col = "lightsteelblue", border = NA)
  lines(orig.pred, y, lwd = 3, col = "blue")
  
}


# Random effect plots -----------------------------------------------------
random_effect_plots <- function(coeff_name,n,main_lab,xlims){
  
  pm <- apply(mod.out.F$sims.list[[coeff_name]], 2, mean)    # Get posterior means and 95% CRIs
  cri <- apply(mod.out.F$sims.list[[coeff_name]], 2, function(x) quantile(x, prob = c(0.025, 0.975))) # CRIs
  
  plot(pm[1:n], 1:n, xlim = xlims, 
       xlab = "Parameter estimate", ylab = "Species number", 
       main = main_lab, pch = 16, type = "n")
  abline(v = 0, lwd = 2, col = "black")
  segments(cri[1, 1:n], 1:n, cri[2, 1:n], 1:n, col = "grey", lwd = 2)
  sig1 <- (cri[1, 1:n] * cri[2, 1:n]) > 0
  segments(cri[1, 1:n][sig1 == 1], (1:n)[sig1 == 1], cri[2, 1:n][sig1 == 1], (1:n)[sig1 == 1], col = "blue", lwd = 2)
  points(pm[1:n], 1:n, pch = 16)
  
  abline(v = mod.out.F$summary[paste0("mu.",coeff_name),1], col = "red", lwd = 3)
  abline(v = mod.out.F$summary[paste0("mu.",coeff_name),c(3,7)], col = "red", lwd = 3, lty = 2)
}
