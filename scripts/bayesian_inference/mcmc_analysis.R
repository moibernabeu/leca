# Functions to analyse MCMC
# Moisès Bernabeu
# Barcelona, November 2023

# Requiring packages ----
require(ggmcmc)
require(tidyr)
require(dplyr)
require(cumstats)
require(coda)
require(ggpubr)

theme_set(theme_bw())

# Defining function ----
get_ac_df <- function(mcmc_df, nLags = 50) {
  # Retrieve autocorrelation per parameter and chain
  wc.ac <- mcmc_df %>%
    group_by(Parameter, Chain) %>%
    do(ac(.$value, nLags))
  
  return(wc.ac)
}


get_full_means <- function(mcmc_df) {
  # Estimates of the mean by iteration
  dm.m <- mcmc_df %>%
    group_by(Parameter, Chain) %>%
    summarize(m = mean(value))
  
  return(dm.m)
}

get_running <- function(mcmc_df) {
  # Calculate the running mean
  # Force the object to be sorted by Parameter, and hence avoid 'rm' calculation
  # to be wrong
  dm.rm <- mcmc_df %>%
    arrange(Parameter, Iteration) %>%
    group_by(Parameter, Chain) %>%
    mutate(rm = cumsum(value) / Iteration,
           sd = sqrt(cumvar(value)),
           up = rm + sd,
           do = rm - sd)
  
  return(dm.rm)
}

mcmc_to_df <- function(mcmc) {
  # Converts an MCMC list object into a dataframe to be used in ggplot
  post_list <- mcmc.list(mcmc)
  post_df <- as.data.frame(ggs(post_list))
  post_df$Chain <- as.factor(post_df$Chain)
  
  return(post_df)
}

get_distro_sum <- function(distros) {
  summlist <- c()
  i = 1
  for (distro in distros) {
    dname <- names(distros)[i]
    summ <- summary(distro[['clean']])
    if (dname == 'gamma') {
      summlist[[dname]] <- summ$statistics[1:2, 'Mean']
      names(summlist[[dname]]) <- c('shape', 'rate')
    } else if (dname == 'lnorm') {
      summlist[[dname]] <- summ$statistics[c('mu', 'sig'), 'Mean']
      names(summlist[[dname]]) <- c('shape', 'rate')
    } else {
      summlist[[dname]] <- summ$statistics[, 'Mean']
    }
    i = i + 1
  }
  
  return(summlist)
}

plot_distros <- function(y, distros) {
  require(see)
  
  x <- 0:(max(y, 0.9)*1000)/1000
  distdat <- c()
  for (distro in names(distros)) {
    params <- distros[[distro]]
    print(paste0('d', distro))
    z <- get(paste0('d', distro))(x, params[1], params[2])
    pdat <- data.frame(x, z, distribution = distro)
    distdat <- rbind(distdat, pdat)
  }
  
  p <- ggplot() +
    geom_histogram(aes(x = y, y = after_stat(density)), colour = 'black', alpha = 0.4) +
    geom_line(aes(x, z, colour = distribution, lty = distribution), distdat, size = 1) +
    scale_colour_okabeito() +
    xlab('Distance') +
    ylab('Density') +
    labs(colour = 'Distribution', lty = 'Distribution')
  
  return(p)
}

plot_cdf <- function(y, distros) {
  require(see)

  x <- 0:(max(y, 0.9)*1000)/1000
  distdat <- c()
  for (distro in names(distros)) {
    params <- distros[[distro]]
    print(paste0('p', distro))
    z <- get(paste0('p', distro))(x, params[1], params[2])
    pdat <- data.frame(x, z, distribution = distro)
    distdat <- rbind(distdat, pdat)
  }
  
  p <- ggplot() +
    stat_ecdf(aes(x = y), colour = 'black', geom = 'point') +
    geom_line(aes(x, z, colour = distribution, lty = distribution), distdat, size = 1) +
    scale_colour_okabeito() +
    xlab('Distance') +
    ylab('CDF') +
    labs(colour = 'Distribution', lty = 'Distribution')
  
  return(p)
}

plot_qq <- function(y, distros) {
  require(see)
  
  dat <- data.frame(y = y)
  p <- ggplot(dat, aes(sample = y)) +
    geom_abline(slope = 1) +
    xlab('Theoretical') +
    ylab('Observed')

  for (distro in names(distros)) {
    params <- distros[[distro]]
    names(params) <- NULL
    print(paste0('p', distro))
    p <- p + geom_qq(distribution = get(paste0('q', distro)),
                     dparams = params,
                     colour = okabeito_colors(which(names(distros) == distro)))
  }
  p
  
  return(p)
}

ggtraces <- function(mcmc_df, title = '') {
  # Plots the MCMC traces, checking whether the file is an MCMC, if it is, the
  # function will convert it using mcmc_to_df
  if (class(mcmc_df) %in% c('mcmc.list', 'mcmc')) {
    mcmc_df <- mcmc_to_df(mcmc_df)
  }
  
  gtraces <- ggplot(mcmc_df, aes(x = Iteration, y = value, colour = Chain)) +
    geom_line(size = 0.2, aes(lty = Chain)) +
    facet_wrap(~Parameter, scales = 'free', ncol = 1) +
    labs(title = title)
  
  return(gtraces)
}

ggmcmcdens <- function(mcmc_df, title = '') {
  # Plots the MCMC densities, checking whether the file is an MCMC, if it is, the
  # function will convert it using mcmc_to_df
  if (class(mcmc_df) %in% c('mcmc.list', 'mcmc')) {
    mcmc_df <- mcmc_to_df(mcmc_df)
  }
  
  ghists <- ggplot(mcmc_df, aes(x = value, colour = Chain)) +
    geom_density() +
    facet_wrap(~Parameter, scales = 'free', ncol = 1) +
    labs(title = title)
  
  return(ghists)
}

ggmcmcsummary <- function(mcmc_df, cl_mcmc_df, title = '') {
  # Plots the MCMC densities, checking whether the file is an MCMC, if it is, the
  # function will convert it using mcmc_to_df
  p <- ggarrange(ggtraces(mcmc_df), ggtraces(cl_mcmc_df), ggmcmcdens(cl_mcmc_df),
                 ncol = 3, align = 'hv', common.legend = TRUE,
                 legend = 'right')
  if (title != '') {
    p <- annotate_figure(p, top = text_grob(title, face = "bold", size = 14))
  }
  
  return(p)
}

numeric_diagnostics <- function(mcmc) {
  Rhat <- gelman.diag(mcmc)
  summ <- summary(mcmc)
  
  efs <- try(effectiveSize(mcmc))
  if (class(efs) == 'try-error') {
    efs_chains <- c()
    for (i in 1:length(mcmc)) {
      try(
        efs_chains <- rbind(efs_chains, effectiveSize(mcmc[[i]]))
      )
    }
    efs <- colMeans(efs_chains, na.rm = TRUE)
  }

  summdf <- cbind(summ$statistics[, 1:2],
                  'R' = Rhat$psrf[, 1],
                  'Effective size' = efs)
  
  return(summdf)
}

ggacplot <- function(mcmc_df, title = '') {
  if (class(mcmc_df) %in% c('mcmc.list', 'mcmc')) {
    mcmc_df <- mcmc_to_df(mcmc_df)
  }
  
  post_ac <- get_ac_df(mcmc_df, 100)
  ac <- ggplot(post_ac, aes(x = Lag, y = Autocorrelation, colour = Chain)) +
    geom_line() +
    facet_grid(Parameter ~ Chain) +
    labs(title = title)
  
  return(ac)
}

graphic_diagnostics <- function(y, distros, prefix, group) {
  distros_summ <- get_distro_sum(distros)

  title <- paste0(group, ' sample: ', length(y), ' trees')
  distrstats <- ggarrange(plot_distros(y, distros_summ),
                          plot_cdf(y, distros_summ),
                          plot_qq(y, distros_summ), nrow = 1, common.legend = TRUE)
  distrstats <- annotate_figure(distrstats, top = text_grob(title, color = 'black',
                                                            face = 'bold', size = 14))
  
  pdf(sprintf('%s_diag.pdf', prefix), width = 10, height = 3, onefile = FALSE)
  print(distrstats)
  dev.off()
  
  for (distro in names(distros)) {
    raw_post <- distros[[distro]][['raw']]
    cl_post <- distros[[distro]][['clean']]
    title <- paste0(group, ' ', distro, ', sample: ', length(y), ' trees')
    
    distsumm <- ggmcmcsummary(raw_post, cl_post, title)
    pdf(sprintf('%s_%s_plots_gg.pdf', prefix, distro), width = 7.8, height = 8.39, onefile = FALSE)
    print(distsumm)
    dev.off()
    
    pac <- ggacplot(gmcmcout, paste0('Raw posterior sample: ', title))
    cl_pac <- ggacplot(cl_post, paste0('Burned and thinned sample: ', title))
    pdf(sprintf('%s_gam_ac.pdf', prefix), width = 7, height = 7)
    print(pac)
    print(cl_pac)
    dev.off()
  }
}
