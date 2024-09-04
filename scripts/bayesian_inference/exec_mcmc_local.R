# MCMC sampling execution
# Moisès Bernabeu
# Barcelona, November 2023

library(stringr)
library(coda)

# Loading functions ----
# source('/gpfs/projects/bsc40/current/mgil/documents/mammals_dating_newphylome/03_inference/scripts/sampling.R')
# source('/gpfs/projects/bsc40/current/mgil/documents/mammals_dating_newphylome/03_inference/scripts/mcmc_analysis.R')

source('/data/moises/mounted/cluster/documents/mammals_dating_newphylome/03_inference/scripts/sampling.R')
source('/data/moises/mounted/cluster/documents/mammals_dating_newphylome/03_inference/scripts/mcmc_analysis.R')

# source('sampling.R')
# source('mcmc_analysis.R')

# Loading data
args = commandArgs(trailingOnly=TRUE)

datfile = args[1] # '../../06_inference/data/dist_data.RData'
object = args[2] # 'fgnew'
groupcol = args[3] # 'node'
samplecol = args[4] # 'ndist'

# Setting MCMC parameters
nchains = as.numeric(args[5]) # 3
niter = as.numeric(args[6]) # 100000
thin = as.numeric(args[7]) # 1
unifmax = as.numeric(args[8]) # 100
burnin = as.numeric(args[9]) # 0.1
subspl = as.numeric(args[10]) # 1

# Setting the species to
group = args[11] # T

# Setting the output folder
folder = args[12] # ../test/T

# datfile = '../../02_dist_stats/data/dist_data.RData'
# object = 'fgnew'
# groupcol = 'node'
# samplecol = 'ndist'
# 
# # Setting MCMC parameters
# nchains = 3
# niter = 10000
# thin = 1
# unifmax = 100
# burnin = 0.1
# subspl = 0.1
# 
# # Setting the species to
# group = 'T'
# 
# # Setting the output folder
# folder = '../test_ln/T'

prefix_mcmc = sprintf('%s/mcmc/%s_%s', folder, group, str_replace(subspl, '\\.', ''))
prefix_plots = sprintf('%s/plots/%s_%s', folder, group, str_replace(subspl, '\\.', ''))
prefix_tables = sprintf('%s/tables/%s_%s', folder, group, str_replace(subspl, '\\.', ''))

dir.create(folder, recursive = TRUE)
dir.create(paste(folder, 'mcmc', sep = '/'), folder, recursive = TRUE)
dir.create(paste(folder, 'plots', sep = '/'), folder, recursive = TRUE)
dir.create(paste(folder, 'tables', sep = '/'), folder, recursive = TRUE)

# Loading data ---
load(datfile)

# Selecting sample ----
dat <- get(object)
y <- dat[which(dat[, groupcol] == group), samplecol]

# Checking if there is any 0 distance, which crashes the inference and their
# proportion is not high
if (sum(y == 0) > 0) {
  y <- y[-which(y == 0)]
}
y0len <- length(y)

# Getting subsamples ----
if (subspl != 1) {
  y <- sample(y, round(length(y) * subspl))
}
y1len <- length(y)

# y <- y / 2

# Finishing when the subsampling implies less than 50 trees ---
if (length(y) < 50) {
  stop(paste0('Not enough sequences trees to compute the inference ',
              'in this subsampling: ', group, '-', subspl, '\n',
              'Full sample: ', y0len, ' trees. Subsample: ', y1len))
}

# Executing sampling
system.time(
  gmcmcout <- gmcmcfun(y = y, nchains, niter, 1, unifmax, group = group,
                       prefix = prefix_mcmc, subspl = subspl)
)

# system.time(
#   nmcmcout <- nmcmcfun(y = y, nchains, niter, 1, unifmax, group = group,
#                        prefix = prefix_mcmc, subspl = subspl)
# )

system.time(
  lnmcmcout <- lnmcmcfun(y = y, nchains, niter, 1, unifmax, group = group,
                         prefix = prefix_mcmc, subspl = subspl)
)

# Cleaning samples with the burn-in and thinning ----
cl_gmcmcout <- window(gmcmcout, start = niter * burnin, end = niter, thin = thin)
# cl_nmcmcout <- window(nmcmcout, start = niter * burnin, end = niter, thin = thin)
cl_lnmcmcout <- window(lnmcmcout, start = niter * burnin, end = niter, thin = thin)

# Plotting inference results ----
if ('fgnew' %in% ls()) {
  # nodes <- read.csv('/gpfs/projects/bsc40/current/mgil/documents/mammals_dating_newphylome/02_dist_stats/data/node_labs.tsv', sep = '\t', row.names = 3)
  nodes <- read.csv('/data/moises/mounted/cluster/documents/mammals_dating_newphylome/02_dist_stats/data/node_labs.tsv', sep = '\t', row.names = 3)
  title <- nodes[group, 'Group']
} else {
  title <- group
}

distros <- list('gamma' = list('raw' = gmcmcout, 'clean' = cl_gmcmcout),
                'lnorm' = list('raw' = lnmcmcout, 'clean' = cl_lnmcmcout))
                # 'norm' = list('raw' = nmcmcout, 'clean' = cl_nmcmcout))

graphic_diagnostics(y, distros, prefix_plots, title, subspl)

# Writing numeric diagnostics tables ----
write.csv(numeric_diagnostics(cl_gmcmcout), file = sprintf('%s_gam.csv', prefix_tables))
# write.csv(numeric_diagnostics(cl_nmcmcout), file = sprintf('%s_nor.csv', prefix_tables))
write.csv(numeric_diagnostics(cl_lnmcmcout), file = sprintf('%s_lnor.csv', prefix_tables))

