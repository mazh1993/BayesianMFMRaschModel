
#### model ####
library(nimble)
library(coda)
piFun2 <- nimbleFunction(
  run = function(r = double(1)) {
    rlength <- length(r)
    rsum <- rep(0, rlength)
    pi <- rep(0, rlength)
    rsum[1] <- pi[1] <- r[1]
    for (i in 2:rlength) {
      rsum[i] <- rsum[i - 1] + r[i]
      if (rsum[i] >= 1) {
        pi[i] <- 1 - rsum[i - 1]
      }
      else {pi[i] <- r[i]}
    }
    for (i in 1:rlength) {
      if (pi[i] < 0) {pi[i] <- 0}
    }
    returnType(double(1))
    return(pi)
  }
)

HeterIRTCode <- nimbleCode({
  for (i in 1:N_ind) {
    for (j in 1:N_item) {
      y[i, j] ~ dbern(py[i, j])
      logit(py[i, j]) <- theta[i] - b[j]
      
      logY[i, j] <- y[i, j] * log(py[i, j]) + (1 - y[i, j]) * log(1 - py[i, j])
    }
    DevY0[i] <- sum(logY[i, 1:N_item])
    
    theta[i] <- thetam[latent_t[i]]
    latent_t[i] ~ dcat(pt[1:M])
  }
  
  for (j in 1:N_item) {
    b[j] <- bm[latent_b[j]]
    latent_b[j] ~ dcat(pb[1:M])
  }
  
  for (l in 1:M) {
    rb[l] ~ dexp(rate = lambda_b)
    rt[l] ~ dexp(rate = lambda_t)
  }
  
  pb[1:M] <- piFun2(r = rb[1:M])
  pt[1:M] <- piFun2(r = rt[1:M])
  
  lambda_b ~ dlnorm(0, varlog = 1)
  lambda_t ~ dlnorm(0, varlog = 1)

  for (h in 1:M) {
    bm[h] ~ dnorm(0, 100)
    thetam[h] ~ dnorm(0, psi)
  }
  psi ~ dgamma(100, 1)
})

ydata <- as.matrix(bookid1_data[,c(10:49)])

run <- function(seed) {
  HeterIRTdata <- list(y = ydata)
  HeterIRTconsts <- list(M = 10, N_ind = nrow(HeterIRTdata$y), N_item = ncol(HeterIRTdata$y))
  HeterIRTinits <- list(
    thetam = rep(0.1, HeterIRTconsts$M), 
    bm = rep(0.5, HeterIRTconsts$M),
    lambda_b = 1, lambda_t = 1, psi = 1,
    rb = rep(0.2, HeterIRTconsts$M),
    rt = rep(0.2, HeterIRTconsts$M), 
    latent_b = rep(1, HeterIRTconsts$N_item), 
    latent_t = rep(1, HeterIRTconsts$N_ind))
  HeterIRT <- nimbleModel(code = HeterIRTCode, name = "HeterIRT", constants = HeterIRTconsts,
                          data = HeterIRTdata, inits = HeterIRTinits)
  cHeterIRT <- compileNimble(HeterIRT)
  HeterIRTconf <- configureMCMC(HeterIRT, print = FALSE, enableWAIC = TRUE)
  HeterIRTconf$addMonitors(c("b","theta", "latent_b", "latent_t",
                             "lambda_b","lambda_t",
                             "DevY0", "py", "psi",
                             "bm", "thetam"))
  HeterIRTmcmc <- buildMCMC(HeterIRTconf)
  cHeterIRTmcmc <- compileNimble(HeterIRTmcmc, project = HeterIRT)
  cHeterIRT$setInits(HeterIRTinits)
  mcmc.out <- runMCMC(cHeterIRTmcmc, niter = 60000, setSeed = seed, thin = 2, WAIC = TRUE)
  pos_mcmc <- as.mcmc(mcmc.out$samples[-c(1:20000), ])
  WAIC <- mcmc.out$WAIC
  return(list(WAIC,pos_mcmc))
}

