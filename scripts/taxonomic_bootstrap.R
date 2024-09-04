# Summarise taxonomic bootstrap
# Moisès Bernabeu
# Cocentaina, August 2024

library(tidyverse)
library(patchwork)

theme_set(theme_bw())

db <- 'TOLDBA'

dat <- read.csv(paste0('../outputs/taxonomic_bootstrap/', db, '.tsv'), sep = '\t')

a <- dat %>%
  filter(ML_donor == 'Alphaproteobacteria', gamma > 0) %>%
  ggplot(aes(x = gamma)) +
  geom_histogram(fill = 'steelblue', colour = 'steelblue', alpha = 0.6) +
  xlab('Proportion of BS trees assigned to Gamma') +
  labs(title = 'ML trees assigned to Alpha (111)')

b <- dat %>%
  filter(ML_donor == 'Gammaproteobacteria', gamma > alpha) %>%
  ggplot(aes(x = alpha)) +
  geom_histogram(fill = 'steelblue', colour = 'steelblue', alpha = 0.6) +
  xlab('Proportion of BS trees assigned to Alpha') +
  # geom_vline(xintercept = 0.054) +
  labs(title = 'ML trees assigned to Gamma (80)')

c <- ggplot() +
  geom_histogram(aes(x = na.omit(dat$alpha - dat$gamma)),
                 fill = 'steelblue', colour = 'steelblue', alpha = 0.6) +
  xlab('Alpha BS trees proportion - Gamma BS trees proportion') +
  labs(title = 'ML trees assigned to any donor (359)')

d <- dat %>%
  filter(ML_donor %in% c('Alphaproteobacteria', 'Gammaproteobacteria'), alpha > 0, gamma > 0) %>%
  ggplot(aes(x = alpha - gamma)) +
  geom_histogram(fill = 'steelblue', colour = 'steelblue', alpha = 0.6) +
  xlab('Alpha BS trees proportion - Gamma BS trees proportion') +
  labs(title = 'ML trees assigned to Gamma or Alpha (199)')

a + b + c + d +
  plot_annotation(tag_levels = 'a', theme = theme(text = element_text(face = 'bold')))


orig <- read.csv(paste0('../outputs/mLECA_origins/', db, '_stats.tsv'), sep = '\t')


