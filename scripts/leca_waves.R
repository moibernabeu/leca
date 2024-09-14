# LECA acquisition waves
# Moisès Bernabeu
# Cocentaina, August, 2023

library(coda)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(tidyverse)

theme_set(theme_bw())

colours <- read.csv('../data/colours.tsv', sep = '\t')[, 2]
names(colours) <- read.csv('../data/colours.tsv', sep = '\t')[, 1]

# niter = 10000
# thin = 10
# burnin = 0.1
# 
# mcmcfiles <- list.files('../outputs/inference/', recursive = TRUE,
#                         pattern = 'gam_mcmc.RData', full.names = TRUE)
# 
# mo <- c()
# for (file in mcmcfiles){
#   load(file)
#   mcmc <- window(get(grep('_1_gam', ls(), value = TRUE)),
#                  start = niter * burnin, end = niter, thin = thin)
# 
#   db <- str_split(file, '/', simplify = TRUE)[, 5]
#   crit <- str_replace(str_split(file, '/', simplify = TRUE)[, 6], '_', ' ')
#   clade <- gsub('_.*', '', grep('_1_gam', ls(), value = TRUE))
#   print(paste(db, crit, clade))
# 
#   for (i in 1:nchain(mcmc)) {
#     mo <- rbind(mo, data.frame(mode = mcmc[[i]][, 'mo'], clade = clade, databse = db, criterion = crit))
#   }
#   remove(list = grep(clade, ls(), value = TRUE))
# }
# 
# write.table(mo, file = '../outputs/inference/mode_posteriors.tsv',
#             sep = '\t', row.names = FALSE, quote = FALSE)


mo <- read.csv('../outputs/inference/mode_posteriors.tsv', sep = '\t')
mod <- mo %>%
  filter(var1 != 0)

mod$criterion <- factor(mod$criterion, levels = c('Three supergroups', 'Five supergroups'))

pdf('../outputs/inference/all_waves.pdf', height = 5.5, width = 7)
ggplot(mod %>% filter(clade != 'Thermoproteota'), aes(var1, colour = clade)) +
  geom_density() +
  facet_grid(databse ~ criterion) +
  xlab('Stem length mode') +
  ylab('Posterior density') +
  labs(colour = 'Donor') +
  scale_colour_manual(values = colours)
dev.off()

clades <- unique(mod$clade)
dbs <- unique(mod$databse)
crits <- unique(mod$criterion)

post_diff_mean <- array(0, dim = c(length(clades), length(clades), length(dbs), length(crits)))

dimnames(post_diff_mean) <- list(clades, clades, dbs, crits)

for (db in dbs) {
  for (crit in crits) {
    for (i in clades) {
      for (j in clades) {
        if (i != j) {
          print(c(db, crit, i, j))
          imod <- mod[(mod$clade == i & mod$databse == db & mod$criterion == crit), 1]
          jmod <- mod[(mod$clade == j & mod$databse == db & mod$criterion == crit), 1]
          ijdf <- imod - jmod
          if (mean(ijdf) > 0) {
            post_diff_mean[i, j, db, crit] <- sum(ijdf > 0) / length(ijdf)
          } else if (mean(ijdf) < 0) {
            post_diff_mean[i, j, db, crit] <- sum(ijdf < 0) / length(ijdf)
          }
        }
      }
    }
  }
}

for (i in clades) {
  for (j in clades) {
    post_diff_mean
  }
}

ord <- mod %>%
  filter(databse == 'TOLDBA', criterion == 'Five supergroups') %>%
  group_by(clade) %>%
  summarise(mod = mean(var1)) %>%
  arrange(-mod)
ord <- ord$clade

mean_post <- apply(post_diff_mean, c(1, 2), median)[ord, ord]
diag(mean_post) <- NA
mean_post[lower.tri(mean_post)] <- NA

pheatmap(mean_post, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE,
         cellwidth = 25, cellheight = 25,
         color = colorRampPalette(c('white', 'steelblue3'))(100),
         na_col = 'white')

pdf('../outputs/inference/prob_heatmap.pdf', width = 8, height = 8)
for (db in dbs) {
  for (crit in crits) {
    ord <- mod %>%
      filter(databse == db, criterion == crit) %>%
      group_by(clade) %>%
      summarise(mod = mean(var1)) %>%
      arrange(-mod)
    ord <- ord$clade
    
    m <- post_diff_mean[ord, ord, db, crit]
    m[lower.tri(m)] <- NA
    diag(m) <- NA
    
    pheatmap(m, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE,
             cellwidth = 25, cellheight = 25,
             color = colorRampPalette(c('white', 'steelblue3'))(100),
             main = paste(db, crit, 'Prob. of being older', sep = ' | '), na_col = 'white')
  }
}
dev.off()
