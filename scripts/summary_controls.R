# Behaviour and certainty of the donors
# Moisès Bernabeu
# Cocentaina, August 2024

library(tidyverse)
library(patchwork)
library(ggh4x)

theme_set(theme_bw())

db <- 'TOLDBC'

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
      # lab <- paste0(names(levs)[i], ' BS - ', names(levs)[j], ' TaxBS - ',
      #               names(levs)[k], ' DonSup')
      lab <- paste0(names(levs)[i], ' - ', names(levs)[j], ' - ',
                    names(levs)[k])
      for (l in c(3, 5, 7, 9)) {
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
      }
      stress_factors <- c(stress_factors, lab)
    }
  }
}

stress_summary$label <- factor(stress_summary$label, levels = stress_factors)

stress_summary$leca_sister_support <- factor(stress_summary$leca_sister_support, levels = c('Low', 'Medium', 'High'))
stress_summary$btr_ML_donor_prop <- factor(stress_summary$btr_ML_donor_prop, levels = c('Low', 'Medium', 'High'))
stress_summary$donor_prop <- factor(stress_summary$donor_prop, levels = c('Low', 'Medium', 'High'))

stress_summary %>%
  filter(donor %in% c('Alphaproteobacteria', 'Asgardarchaeota')) %>%
  select(donor, supergroups, label, clade_prop) %>% 
  pivot_wider(names_from = donor, values_from = clade_prop) %>%
  mutate(Total = Alphaproteobacteria + Asgardarchaeota) %>%
  pivot_longer(cols = c(Asgardarchaeota, Alphaproteobacteria, Total), names_to = 'donor', values_to = 'clade_prop') %>%
  ggplot(aes(label,
             clade_prop, colour = donor, group = donor)) +
  geom_point() +
  geom_line() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  facet_wrap(supergroups~.) +
  scale_colour_manual(values = c('Alphaproteobacteria' = "#f8766dff",
                                 'Asgardarchaeota' = "#00bfc4ff",
                                 'Total' = "black")) +
  geom_vline(xintercept = 'Medium - High - High', lty = 4, col = 'grey30')

colours <- read.csv('../data/colours.tsv', sep = '\t')[, 2]
names(colours) <- read.csv('../data/colours.tsv', sep = '\t')[, 1]
stress_summary %>%
  group_by(supergroups, label) %>%
  mutate(total_trees = sum(n)) %>%
  ungroup() %>%
  filter(donor %in% c('Alphaproteobacteria', 'Asgardarchaeota')) %>%
  select(donor, supergroups, label, clade_prop, n, total_trees) %>% 
  pivot_wider(names_from = donor, values_from = clade_prop) %>%
  # mutate(Total = Alphaproteobacteria + Asgardarchaeota) %>%
  pivot_longer(cols = c(Asgardarchaeota, Alphaproteobacteria), names_to = 'donor', values_to = 'clade_prop') %>%
  ggplot(aes(label,
             clade_prop, fill = donor, group = donor)) +
  geom_col() +
  geom_point(aes(y = total_trees / 10000), colour = 'grey50') +
  geom_line(aes(y = total_trees / 10000), colour = 'grey50') +
  # geom_line() +
  # theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(axis.text.x = element_blank()) +
  facet_wrap(supergroups~., nrow = 1) +
  scale_fill_manual(values = c(colours, 'Total' = "black")) +
  geom_vline(xintercept = 'Medium - High - High', lty = 4, col = 'grey30') +
  scale_y_continuous(
    
    # Features of the first axis
    name = "Proportion of Alpha/Asgard",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis( trans=~.*10000, name="Number of trees")
  )


# analyses of the donors
donors <- c()
for (db in c('TOLDBA', 'TOLDBB', 'TOLDBC')) {
  taxboot <- read.csv(paste0('../outputs/taxonomic_bootstrap/', db), sep = '\t')
  origins <- read.csv(paste0('../outputs/mLECA_origins/with_proportions/', db, '_stats.tsv'), sep = '\t')
  
  dat <- taxboot %>%
    left_join(origins, by = c('family' = 'Orthogroup', 'LECA' = 'LECA'))
  
  stress_dat <- data.frame(family = dat$family, LECA = dat$LECA, donor = dat$ML_donor, donor_domain = dat$S1_MAB_group_domain, supergroups = dat$supergroups)
  stress_dat$leca_sister_support <- cut(dat$leca_sister_support, c(-1, 75, 95, 101), labels = c('Low', 'Medium', 'High'))
  stress_dat$btr_ML_donor_prop <- cut(dat$btr_ML_donor_prop, c(-1, 0.75, 0.95, 1.01), labels = c('Low', 'Medium', 'High'))
  stress_dat$donor_prop <- cut(dat$donor_prop, quantile((dat %>% filter(S1_size > 5))$donor_prop, c(0, 0.5, 0.75, 1)), include.lowest = TRUE, labels = c('Low', 'Median', 'High'))

  stress_summary <- c()
  stress_factors <- c()
  for (i in 1:3) {
    for (j in 1:3) {
      for (k in 1:3) {
        # lab <- paste0(names(levs)[i], ' BS - ', names(levs)[j], ' TaxBS - ',
        #               names(levs)[k], ' DonSup')
        lab <- paste0(names(levs)[i], ' - ', names(levs)[j], ' - ',
                      names(levs)[k])
        for (l in c(3, 5, 7, 9)) {
          a <- stress_dat %>%
            filter(leca_sister_support %in% levs[[i]],
                   btr_ML_donor_prop %in% levs[[j]],
                   donor_prop %in% levs[[k]],
                   supergroups >= l) %>%
            mutate(supergroups = l) %>%
            group_by(donor, donor_domain, supergroups) %>%
            count() %>%
            ungroup() %>%
            mutate(clade_prop = n / sum(n), label = lab,
                   leca_sister_support = names(levs)[i],
                   btr_ML_donor_prop =  names(levs)[j],
                   donor_prop =  names(levs)[k],
                   supergroups = l)
          stress_summary <- rbind(stress_summary, a)
        }
        stress_factors <- c(stress_factors, lab)
      }
    }
  }
  
  stress_summary$label <- factor(stress_summary$label, levels = stress_factors)
  
  stress_summary$leca_sister_support <- factor(stress_summary$leca_sister_support, levels = c('Low', 'Medium', 'High'))
  stress_summary$btr_ML_donor_prop <- factor(stress_summary$btr_ML_donor_prop, levels = c('Low', 'Medium', 'High'))
  stress_summary$donor_prop <- factor(stress_summary$donor_prop, levels = c('Low', 'Medium', 'High'))
  
  
  a <- stress_summary %>%
    filter(leca_sister_support == 'Medium',
           btr_ML_donor_prop == 'High',
           donor_prop == 'High',
           supergroups == 7,
           !donor %in% c('Viruses', 'Bacteria', 'Archaea')) %>%
    mutate(database = db) %>%
    arrange(-n)
  donors <- rbind(donors, a)
}

stressed_donors <- donors %>%
  select(donor, donor_domain, database, n) %>%
  pivot_wider(names_from = database, values_from = n) %>%
  print(n = 20) %>%
  filter(TOLDBA >= 15 | TOLDBB >= 15 | TOLDBC >= 15)

write.table(stressed_donors[, c(2, 1)], file = '../outputs/LECA_proteomes/selected_groups.tsv', sep = '\t', row.names = FALSE, quote = FALSE)

donors %>%
  select(donor, database, n) %>%
  pivot_wider(names_from = database, values_from = n) %>%
  print(n = 20) %>%
  filter(TOLDBA >= 3 | TOLDBB >= 3 | TOLDBC >= 3) %>%
  na.omit() %>%
  ggplot(aes(TOLDBA, reorder(donor, TOLDBA), colour = TOLDBA >= 15)) +
  geom_point() +
  geom_vline(xintercept = 12)

stress_summary_clades <- c()
stress_factors <- c()
for (i in 1:3) {
  for (j in 1:3) {
    for (k in 1:3) {
      for (l in c(3, 5, 7, 9)) {
        # lab <- paste0(names(levs)[i], ' BS - ', names(levs)[j], ' TaxBS - ',
        #               names(levs)[k], ' DonSup - ', l, ' sg')
        lab <- paste0(names(levs)[i], ' - ', names(levs)[j], ' - ',
                      names(levs)[k])
        a <- stress_dat %>%
          filter(leca_sister_support %in% levs[[i]],
                 btr_ML_donor_prop %in% levs[[j]],
                 donor_prop %in% levs[[k]],
                 supergroups >= l) %>%
          group_by(donor) %>%
          count() %>%
          mutate(leca_sister_support = names(levs)[i],
                 btr_ML_donor_prop =  names(levs)[j],
                 donor_prop =  names(levs)[k],
                 supergroups = l,
                 label = lab)
        stress_summary_clades <- rbind(stress_summary_clades, a)
      }
      stress_factors <- c(stress_factors, lab)
    }
  }
}

stress_summary_clades <- stress_summary_clades %>%
  group_by(leca_sister_support, btr_ML_donor_prop, donor_prop, supergroups, label) %>%
  mutate(clade_prop = n / sum(n)) %>%
  summarise(max_prop = max(clade_prop),
            min_prop = min(clade_prop),
            max_ntree = max(n),
            min_ntree = min(n),
            q25 = quantile(clade_prop, 0.25),
            q50 = quantile(clade_prop, 0.5),
            q75 = quantile(clade_prop, 0.75),
            ndonors = length(unique(donor)))


stress_summary_clades$label <- factor(stress_summary_clades$label, stress_factors)
stress_summary_clades %>%
  pivot_longer(cols = c(max_prop, min_prop, max_ntree, ndonors), names_to = 'variable') %>%
  ggplot(aes(label, value, colour = supergroups)) +
  geom_point() +
  facet_grid(~variable, scales = 'free_x') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  coord_flip() +
  scale_colour_viridis_c(direction = -1)

