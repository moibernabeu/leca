# MCMC sampling functions using JAGS
# Moisès Bernabeu
# Barcelona, May 2022

gmcmcfun <- function(y, nchains = 3, niter = 1000,
                     thin = 1, unifmax = 100, group, prefix, subspl, burnin = 0) {
  require(rjags)

  # Model definition
  model_text <- sprintf('model{
    # Likelihood
    for (i in 1:n) {
      y[i] ~ dgamma(a, b)
    }
    
    m <- a / b
    v <- a / b^2
    mo <- ifelse(a <= 1, 0, (a - 1) / b)
    
    # Prior distributions
    a ~ dunif(0, %d)
    b ~ dunif(0, %d)
  }', unifmax, unifmax)
  
  # Data generation
  N <- length(y)
  
  datlist <- list(y = y, n = N)
  
  # Preparing the model
  model_jags <- jags.model(textConnection(model_text), 
                           data = datlist,
                           n.chains = nchains,
                           n.adapt = round(niter * burnin, digits = 0))
  
  # Running the MCMC sampling
  post <- coda.samples(model_jags, 
                       variable.names = c('a', 'b', 'm', 'v', 'mo'), 
                       n.iter = niter,
                       thin = thin)
  
  # Assigning species name to the posterior sample
  assign(sprintf('%s_%s_gam', group, subspl), post)
  assign(sprintf('%s_%s_y', group, subspl), y)
  
  # Saving the posterior sample and its information into RData file
  save(list = c('N', 'nchains', 'niter', 'thin', 'unifmax',
                sprintf('%s_%s_gam', group, subspl), sprintf('%s_%s_y', group, subspl)),
       file = sprintf('%s_gam_mcmc.RData', prefix))
  
  return(post)
}

nmcmcfun <- function(y, nchains = 3, niter = 1000,
                     thin = 1, unifmax = 100, group, prefix, subspl, burnin = 0) {
  require(rjags)

  # Model definition
  model_text <- sprintf('model{
    # Likelihood
    for (i in 1:n) {
      y[i] ~ dnorm(mu, tau)
    }
    
    tau <- 1 / sig ^ 2
    
    # Prior distributions
    mu ~ dunif(0, %d)
    sig ~ dunif(0, %d)
  }', unifmax, unifmax)
  
  # Data generation
  N <- length(y)
  
  datlist <- list(y = y, n = N)
  
  # Preparing the model
  model_jags <- jags.model(textConnection(model_text), 
                           data = datlist,
                           n.chains = nchains,
                           n.adapt = round(niter * burnin, digits = 0))
  
  # Running the MCMC sampling
  post <- coda.samples(model_jags, 
                       variable.names = c('mu', 'sig'), 
                       n.iter = niter,
                       thin = thin)
  
  # Assigning species name to the posterior sample
  assign(sprintf('%s_%s_nor', group, subspl), post)
  assign(sprintf('%s_%s_y', group, subspl), y)
  
  # Saving the posterior sample and its information into RData file
  save(list = c('N', 'nchains', 'niter', 'thin', 'unifmax',
                sprintf('%s_%s_nor', group, subspl), sprintf('%s_%s_y', group, subspl)),
       file = sprintf('%s_nor_mcmc.RData', prefix))
  
  return(post)
}

lnmcmcfun <- function(y, nchains = 3, niter = 1000,
                      thin = 1, unifmax = 100, group, prefix, subspl, burnin = 0) {
  require(rjags)

  # Model definition
  model_text <- sprintf('model{
    # Likelihood
    for (i in 1:n) {
      y[i] ~ dlnorm(mu, tau)
    }
    
    tau <- 1 / sig ^ 2
    
    m <- exp(mu + (sig ^ 2))
    v <- (exp(sig ^ 2) - 1) * exp(2 * mu + sig ^ 2)
    mo <- exp(mu - sig ^ 2)
    
    # Prior distributions
    mu ~ dunif(0, %d)
    sig ~ dunif(0, %d)
  }', unifmax, unifmax)
  
  # Data generation
  N <- length(y)
  
  datlist <- list(y = y, n = N)
  
  # Preparing the model
  model_jags <- jags.model(textConnection(model_text), 
                           data = datlist,
                           n.chains = nchains,
                           n.adapt = round(niter * burnin, digits = 0))
  
  # Running the MCMC sampling
  post <- coda.samples(model_jags, 
                       variable.names = c('mu', 'sig', 'm', 'v', 'mo'), 
                       n.iter = niter,
                       thin = thin)
  
  # Assigning species name to the posterior sample
  assign(sprintf('%s_%s_nor', group, subspl), post)
  assign(sprintf('%s_%s_y', group, subspl), y)
  
  # Saving the posterior sample and its information into RData file
  save(list = c('N', 'nchains', 'niter', 'thin', 'unifmax',
                sprintf('%s_%s_nor', group, subspl), sprintf('%s_%s_y', group, subspl)),
       file = sprintf('%s_nor_mcmc.RData', prefix))
  
  return(post)
}
