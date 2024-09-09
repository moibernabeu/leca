# Behaviour and certainty of the donors
# Moisès Bernabeu
# Cocentaina, August 2024

library(tidyverse)
library(patchwork)
library(ggh4x)

theme_set(theme_bw())

db <- 'TOLDBA'

taxboot <- read.csv(paste0('../outputs/taxonomic_bootstrap/', db), sep = '\t')
origins <- read.csv(paste0('../outputs/mLECA_origins/with_proportions/', db, '_stats.tsv'), sep = '\t')

dat <- taxboot %>%
  left_join(origins, by = c('family' = 'Orthogroup', 'LECA' = 'LECA'))

# LECA bootstrap
a <- ggplot(dat, aes(x = LECA_bs)) +
  geom_histogram(fill = 'steelblue', colour = 'steelblue', alpha = 0.6) +
  ylab('Number of trees') +
  xlab('LECA UltrafastBS')

# LECA-sister UfBS support
b <- ggplot(dat, aes(x = leca_sister_support)) +
  geom_histogram(fill = 'steelblue', colour = 'steelblue', alpha = 0.6) +
  ylab('Number of trees') +
  xlab('LECA-sister UltrafastBS')

# LECA certainty
c <- ggplot(dat, aes(x = LECA_jaccard_mean)) +
  geom_histogram(fill = 'steelblue', colour = 'steelblue', alpha = 0.6) +
  ylab('Number of trees') +
  xlab('LECA congruence with ML LECA')

# Taxonomic bootstrap of the donor
d <- ggplot(dat %>% filter(S1_size > 5), aes(x = btr_ML_donor_prop)) +
  geom_histogram(fill = 'steelblue', colour = 'steelblue', alpha = 0.6) +
  ylab('Number of trees') +
  xlab('Taxonomic bootstrap')

# most abundant taxa (MAT) proportion
e <- ggplot(dat %>% filter(S1_size > 5), aes(x = donor_prop)) +
  geom_histogram(fill = 'steelblue', colour = 'steelblue', alpha = 0.6) +
  ylab('Number of trees') +
  xlab('Proportion of the donor in S1')

# Nestedness of the donor
f <- ggplot(dat, aes(x = S1_in_S2)) +
  geom_bar(fill = 'steelblue', colour = 'steelblue', alpha = 0.6) +
  ylab('Number of trees') +
  xlab('Donor in S2')

# pdf('../outputs/stress_test/histograms.pdf', width = 9, height = 4.5)
a + b + c + d + e + f +
  plot_annotation(tag_levels = 'a') &
  theme(plot.tag = element_text(face = 'bold'))
# dev.off()

1 - ecdf(dat$LECA_bs)(75)
1 - ecdf(dat$leca_sister_support)(75)
1 - ecdf(dat$LECA_jaccard_mean)(0.95)
ecdf(dat$donor_prop)(0.5)
table(dat$S1_in_S2) / length(dat$S1_in_S2)

# Calculating the thresholds and adding to the table
stress_dat <- data.frame(family = dat$family, LECA = dat$LECA, donor = dat$ML_donor, supergroups = dat$supergroups)
stress_dat$leca_sister_support <- cut(dat$leca_sister_support, c(-1, 75, 95, 101), labels = c('Low', 'Medium', 'High'))
stress_dat$btr_ML_donor_prop <- cut(dat$btr_ML_donor_prop, c(-1, 0.75, 0.95, 1.01), labels = c('Low', 'Medium', 'High'))
stress_dat$donor_prop <- cut(dat$donor_prop, quantile((dat %>% filter(S1_size > 5))$donor_prop, c(0, 0.5, 0.75, 1)), include.lowest = TRUE, labels = c('Low', 'Median', 'High'))



stress_dat <- stress_dat %>%
  filter(!str_detect(donor, ';') | donor %in% c('Archaea', 'Bacteria', 'Viruses'))

levs <- list('Low' = c('Low', 'Medium', 'High'),
             'Medium' = c('Medium', 'High'),
             'High' = 'High')

stress_summary <- c()
stress_factors <- c()
for (i in 1:3) {
  for (j in 1:3) {
    for (k in 1:3) {
      for (l in c(3, 5, 7, 9)) {
        lab <- paste0(names(levs)[i], ' BS - ', names(levs)[j], ' TaxBS - ',
                      names(levs)[k], ' DonSup - ', l, ' sg')
        a <- stress_dat %>%
          filter(leca_sister_support %in% levs[[i]],
                 btr_ML_donor_prop %in% levs[[j]],
                 donor_prop %in% levs[[k]],
                 supergroups >= l) %>%
          mutate(supergroups = l) %>%
          group_by(donor, supergroups) %>%
          count() %>%
          ungroup() %>%
          mutate(clade_prop = n / sum(n), label = lab,
                 leca_sister_support = names(levs)[i],
                 btr_ML_donor_prop =  names(levs)[j],
                 donor_prop =  names(levs)[k],
                 supergroups = l)
        stress_summary <- rbind(stress_summary, a)
        stress_factors <- c(stress_factors, lab)
      }
    }
  }
}

stress_summary$label <- factor(stress_summary$label, stress_factors)

stress_summary$leca_sister_support <- factor(stress_summary$leca_sister_support, labels = c('Low', 'Medium', 'High'))
stress_summary$btr_ML_donor_prop <- factor(stress_summary$btr_ML_donor_prop, labels = c('Low', 'Medium', 'High'))
stress_summary$donor_prop <- factor(stress_summary$donor_prop, labels = c('Low', 'Medium', 'High'))

stress_summary %>%
  filter(donor %in% c('Alphaproteobacteria', 'Asgardarchaeota')) %>%
  ggplot(aes(interaction(leca_sister_support, btr_ML_donor_prop, donor_prop, lex.order = TRUE),
             clade_prop, colour = donor, group = donor)) +
  geom_point() +
  geom_line() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_x_discrete(NULL, guide = "axis_nested") +
  coord_flip() +
  facet_grid(supergroups~.)
