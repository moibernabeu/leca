# MCMC sampling execution
# Moisès Bernabeu
# Barcelona, November 2023

library(stringr)
library(coda)
library(rjags)
library(dplyr)

# Loading functions ----
# source('/gpfs/projects/bsc40/current/mgil/documents/mammals_dating_newphylome/03_inference/scripts/sampling.R')
# source('/gpfs/projects/bsc40/current/mgil/documents/mammals_dating_newphylome/03_inference/scripts/mcmc_analysis.R')

source('sampling.R')
source('mcmc_analysis.R')

# Loading data
dat <- c()
for (i in c('A', 'B', 'C')) {
  dat <- rbind(dat, read.csv(paste0('../../outputs/modules_vircleaned/TOLDB', i, '_modules_annotation_filtered.tsv'), sep = '\t'))
}

# Setting MCMC parameters
nchains = 3
niter = 10000
thin = 1
unifmax = 100
burnin = 0.1

# Setting the output folder
folder = '../../outputs/inference_vircleaned'
dir.create(folder, recursive = TRUE)

db <- 'TOLDBA'
crit <- 'Three_supergroups'
don <- 'Thermoproteota'

for (db in unique(dat$database)) {
  for (crit in unique(dat$criterion)) {
    for (don in unique(dat$donor)) {
      cat(sprintf('\n\n\n------------------------\n%s %s %s\n------------------------\n', db, crit, don))
      basefold = sprintf('%s/%s/%s/%s/', folder, db, gsub(' ', '_', crit), don)
      prefix_mcmc = sprintf('%s/mcmc/%s_%s', basefold, don, niter)
      prefix_plots = sprintf('%s/plots/%s_%s', basefold, don, niter)
      prefix_tables = sprintf('%s/tables/%s_%s', basefold, don, niter)

      dir.create(paste(basefold, 'mcmc', sep = '/'), folder, recursive = TRUE)
      dir.create(paste(basefold, 'plots', sep = '/'), folder, recursive = TRUE)
      dir.create(paste(basefold, 'tables', sep = '/'), folder, recursive = TRUE)
      
      fdat <- dat %>%
        filter(criterion == crit, database == db,
               donor == don, in_module == 'TRUE')

      y <- fdat$sl
      y <- na.omit(y)
      y[is.infinite(y)] <- NULL
      y <- y[!y == 0]

      # Executing sampling
      system.time(
        gmcmcout <- gmcmcfun(y = y, nchains, niter, 1, unifmax, group = don,
                             prefix = prefix_mcmc, subspl = 1)
      )

      cl_gmcmcout <- window(gmcmcout, start = niter * burnin, end = niter, thin = thin)

      distros <- list('gamma' = list('raw' = gmcmcout, 'clean' = cl_gmcmcout))

      graphic_diagnostics(y, distros, prefix_plots, don)
      
      # Writing numeric diagnostics tables ----
      write.csv(numeric_diagnostics(cl_gmcmcout), file = sprintf('%s_gam.csv', prefix_tables))
    }
  }
}
