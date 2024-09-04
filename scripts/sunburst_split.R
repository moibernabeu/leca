# LECA sunburst
# Moisès Bernabeu
# Cocentaina, August 2024

library(tidyverse)
library(plotly)

dat <- c()
for (db in c('TOLDBA', 'TOLDBB', 'TOLDBC')) {
  for (sg in c(3, 5)) {
    fdat <- read.csv(paste0('../outputs/LECA_proteomes_vircleaned/', db, '_', sg, 'sg_proteome.tsv'), sep = '\t')
    fdat <- fdat %>%
      mutate(db = db, criteria = ifelse(sg == 3, 'Relaxed', 'Strict'))
    dat <- rbind(dat, fdat)
  }
}

aa <- dat %>%
  mutate(donor_domain = ifelse(str_detect(donor_domain, ';'), 'Undetermined', donor_domain)) %>%
  group_by(donor_domain, donor) %>%
  mutate(donor = ifelse(donor_domain == donor, paste0('Mixed ', donor), donor)) %>%
  ungroup() %>%
  group_by(criteria) %>%
  mutate(total_criteria = n() / 3) %>%
  group_by(criteria, origin, donor_domain, donor, total_criteria) %>%
  # summarise(mean_genes = n() / 3, prop = round(n() / 3)) %>%
  summarise(mean_genes = n() / 3, prop = n() / 3 / unique(total_criteria) * 100) %>%
  ungroup()

# RELAXED
aa %>%
  group_by(criteria) %>%
  summarise(sum(prop))

origin <- aa %>%
  filter(criteria == 'Relaxed') %>%
  group_by(label = origin) %>%
  summarise(parent = 'LECA', value = sum(prop))

domain <- aa %>%
  filter(criteria == 'Relaxed', origin == 'Acquisition') %>%
  group_by(label = donor_domain, parent = origin) %>% 
  summarise(value = sum(prop))

donor <- aa %>%
  filter(criteria == 'Relaxed', origin == 'Acquisition') %>%
  mutate(donor = ifelse(str_detect(donor, ';'), 'Undetermined donor', donor)) %>%
  group_by(parent = donor_domain, label = donor) %>%
  summarise(value = sum(prop))

sbdat <- rbind(data.frame(parent = c('', 'Criteria', 'Criteria'),
                          label = c('Criteria', 'Relaxed', 'LECA'),
                          value = c(200, 100, 100)),
               origin, domain, donor) %>%
  filter(value > 2)

plot_relaxed <- plot_ly(sbdat, ids = ~label, labels = ~label, parents = ~parent, values = ~value,
                        type = 'sunburst', branchvalues = 'total',  textinfo= 'label',
                        strokes = 'black', stroke = 'black', size = I(6))
plot_relaxed

save_image(plot_relaxed, file = '../outputs/LECA_proteomes/sunbursnt_consensus_relaxed.pdf')

sbdat <- rbind(data.frame(parent = c(''),
                          label = c('LECA'),
                          value = c(100)),
               origin, domain, donor) %>%
  filter(value >= 2)

plot_relaxed <- plot_ly(sbdat, ids = ~label, labels = ~label, parents = ~parent, values = ~value,
                       type = 'icicle', branchvalues = 'total',  textinfo= 'label',
                       strokes = 'black', stroke = 'black', size = I(16))
plot_relaxed

save_image(plot_relaxed, file = '../outputs/LECA_proteomes/icicle_consensus_relaxed.pdf')

# STRICT
aa %>%
  group_by(criteria) %>%
  summarise(sum(prop))

origin <- aa %>%
  filter(criteria == 'Strict') %>%
  group_by(label = origin) %>%
  summarise(parent = 'LECA', value = sum(prop))

domain <- aa %>%
  filter(criteria == 'Strict', origin == 'Acquisition') %>%
  group_by(label = donor_domain, parent = origin) %>% 
  summarise(value = sum(prop))

donor <- aa %>%
  filter(criteria == 'Strict', origin == 'Acquisition') %>%
  mutate(donor = ifelse(str_detect(donor, ';'), 'Undetermined donor', donor)) %>%
  group_by(parent = donor_domain, label = donor) %>%
  summarise(value = sum(prop))

sbdat <- rbind(data.frame(parent = c('', 'Criteria', 'Criteria'),
                          label = c('Criteria', 'Strict', 'LECA'),
                          value = c(200, 100, 100)),
               origin, domain, donor) %>%
  filter(value >= 2)

plot_strict <- plot_ly(sbdat, ids = ~label, labels = ~label, parents = ~parent, values = ~value,
                       type = 'sunburst', branchvalues = 'total',  textinfo= 'label',
                       strokes = 'black', stroke = 'black', size = I(6))
plot_strict

save_image(plot_strict, file = '../outputs/LECA_proteomes/sunbursnt_consensus_strict.pdf')

sbdat <- rbind(data.frame(parent = c(''),
                          label = c('LECA'),
                          value = c(100)),
               origin, domain, donor) %>%
  filter(value >= 2)

plot_strict <- plot_ly(sbdat, ids = ~label, labels = ~label, parents = ~parent, values = ~value,
                       type = 'icicle', branchvalues = 'total',  textinfo= 'label+percent entry',
                       strokes = 'black', stroke = 'black', size = I(16))
plot_strict

save_image(plot_strict, file = '../outputs/LECA_proteomes/icicle_consensus_strict.pdf')
save_image(plot_strict, file = '../paper/figures_R/raw_plots/figure_2/icicle_consensus_strict.pdf')


aa %>%
  group_by(criteria) %>%
  summarise(unique(total_criteria))
