# various utility functions 

# density function for Wald (unshifted)
dwald = function(x,gamma,alpha){
  return((alpha/(sqrt(2*pi*x^3)))*exp(-(alpha-gamma*x)^2/(2*x)))
}

# function to generate random shifted Wald data
# adapted from pp. 79-80, Dagpunar, J. (1988). Principles of Random Variate Generation. Clarendon Press, Oxford.
# code modified from Heathcote (2004)
rwald = function(n, gamma, alpha, theta) {
  y2 = rchisq(n,1)
  y2onm = y2/gamma
  u = runif(n)
  r1 = (2*alpha + y2onm - sqrt(y2onm*(4*alpha+y2onm)))/(2*gamma)
  r2 = (alpha/gamma)^2/r1
  ifelse(u < alpha/(alpha+gamma*r1), theta+r1, theta+r2)
}

## Utility function to assess prior probability that all 
## effects are greater than zero 
prior.p.greater = function(M, I, a = alpha, b = beta, sd_mu = sigma_mu){
  
  s2 = MCMCpack::rinvgamma(M, a, b)
  mu = rcauchy(M, 0, sd_mu)
  res = exp(pnorm(0, mu, sqrt(s2), lower.tail = F, log.p = T) * I)
  
  return(mean(res))
}


# define function to build design matrix and define g-priors
prep.models = function(sub, cond){
  I = length(unique(sub))
  R = length(sub)
  X.full = matrix(nrow = R, ncol = 2 * I + 2, 0)
  for (r in 1:R){
    X.full[r, 1] = 1
    X.full[r, sub[r] + 1] = 1
    if (cond[r] == 2) {
      X.full[r, I + 2] <- 1
      X.full[r, I + 2 + sub[r]] <- 1
    }
  }
  
  gMap.full = c(rep(0, I), 1, rep(2, I))
  
  X.one = matrix(nrow = R, ncol = I + 2, 0)
  for (r in 1:R){
    X.one[r, 1] = 1
    X.one[r, sub[r] + 1] = 1
    if (cond[r] == 2) {
      X.one[r, I + 2] = 1
    }
  }
  
  gMap.one = c(rep(0, I), 1)
  
  X.null = matrix(nrow = R, ncol = I + 1, 0)
  for(r in 1:R){
    X.null[r, 1] = 1
    X.null[r, sub[r] + 1] = 1
  }
  
  gMap.null <- rep(0, I)
  
  return(list(X.full = X.full,
              gMap.full = gMap.full,
              X.one = X.one,
              gMap.one = gMap.one,
              X.null = X.null,
              gMap.null= gMap.null,
              R = R,
              I = I))
}


# define function to compute BFs and sample posteriors
makeBF = function(y, meanScale, effectScale, prep, keep=1001:10000){
  bf.full = nWayAOV(y,
                    prep$X.full,
                    prep$gMap.full,
                    rscale = c(1, meanScale, effectScale),
                    posterior = F,
                    method = "auto",
                    iterations = max(keep))
  
  mcmc.full = nWayAOV(y,
                      prep$X.full,
                      prep$gMap.full,
                      rscale = c(1, meanScale, effectScale),
                      posterior = T,
                      method = "auto",
                      iterations = max(keep))
  
  
  # "encompassing approach"
  # find number of evidential iterations conditional on data
  i.delta0 = prep$I + 2
  i.delta = (prep$I + 3):(2*prep$I + 2)
  
  myDelta = mcmc.full[keep, i.delta0] + mcmc.full[keep, i.delta]
  good = myDelta > 0  # returns TRUE if all samples are positive (i.e., evidential)
  all.good = apply(good, 1, mean) # returns proportion of positive sample values
  PostCount = mean(all.good == 1) # returns proportion of "evidential" samples
  
  #Evaluate how often all individuals' theta_i are positive in the prior
  R = max(keep) * 10
  beta = .5 * effectScale^2
  alpha = .5
  mu.theta.sd = meanScale
  PriorCount = prior.p.greater(M = R, I = prep$I, a = alpha, b = beta, sd_mu = mu.theta.sd)
  
  # compute Bayes factors for positive effect over unconstrained
  bf.pu = PostCount/PriorCount
  
  return(list(bf.pu = bf.pu))
}


