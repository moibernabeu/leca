# LECA sunburst
# Moisès Bernabeu
# Cocentaina, August 2024

library(tidyverse)
library(plotly)

theme_set(theme_bw())

db = 'TOLDBA'
sg = 3

odir <- '../outputs/mLECA_origins_vircleaned'
dir.create(odir, showWarnings = FALSE, recursive = TRUE)

to_exclude <- read.csv('../data/list_codes_Nucleocytoviricota.second_sister.3SG.stats', sep = '\t', fill = TRUE, header = FALSE)
to_exclude <- to_exclude %>%
  mutate(to_exclude = str_detect(V1, '#'), database = gsub('#', '', V1)) %>%
  filter(to_exclude)

for (db in c('TOLDBA', 'TOLDBB', 'TOLDBC')) {
  for (sg in c(3, 5, 7)) {
    dat <- read.csv(paste0('../outputs/LECA_proteomes_vircleaned/', db, '_', sg, 'sg_proteome.tsv'), sep = '\t')
    
    dat <- dat %>%
      filter(supergroups >= sg)
    
    origin <- dat %>%
      select(label = origin) %>%
      mutate(parent = 'LECA') %>%
      group_by(label, parent) %>%
      count(name = 'value')
    
    domain <- dat %>%
      select(label = donor_domain) %>%
      mutate(label = ifelse(str_detect(label, ';'), 'Undetermined', label)) %>%
      mutate(parent = 'Acquisition') %>%
      group_by(label, parent) %>%
      na.omit() %>%
      count(name = 'value')
    
    donor <- dat %>%
      select(donor_domain, donor) %>%
      filter(!str_detect(donor, ';')) %>%
      rename(parent = donor_domain, label = donor) %>%
      group_by(parent, label) %>%
      mutate(label = ifelse(parent == label, paste0('Mixed ', label), label)) %>%
      count(name = 'value') %>%
      arrange(-value)
    
    sbdat <- rbind(data.frame(parent = '', value = sum(origin$value), label = 'LECA'),
                   origin, domain, donor)
    
    plot <- plot_ly(sbdat, ids = ~label, labels = ~label, parents = ~parent, values = ~value,
                    type = 'sunburst', branchvalues = 'total',  textinfo= 'label+value',
                    strokes = 'black', stroke = 'black', size = I(15)) %>%
      layout(title = paste0(db, ' | ', sg, ' supergroups'))

    save_image(plot, paste0(odir, '/', db, '_', sg, 'sg_sunburst.pdf'))
  }
}

# Stress test
all_summaries <- c()
for (db in c('TOLDBA', 'TOLDBB', 'TOLDBC')) {
  dat <- read.csv(paste0('../outputs/mLECA_origins/', db , '_stats.tsv'), sep = '\t')
  
  db_exclude <- to_exclude %>%
    filter(database == db)
  db_exclude <- db_exclude$V2
  
  dat <- dat %>%
    mutate(support = cut(dat$leca_sister_support, c(0, 70, 95, 101), right = F,
                         labels = c('low', 'moderate', 'high'))) %>%
    filter(!paste0(Orthogroup, '_', LECA) %in% db_exclude)

  for (sg in c(3, 5, 7, 9)) {
    for (sist_supp in c('low', 'moderate', 'high')) {
      # for (leca_size in c('small' = 5, 'medium' = 10, 'big' = 20)) {
      fdat <- dat %>%
        filter(supergroups >= sg, support == sist_supp)
      
      fdat_len <- dim(fdat)[1]

      fdat <- fdat %>%
        select(S1_MAB_group_domain, donor) %>%
        filter(!str_detect(donor, ';')) %>%
        rename(parent = S1_MAB_group_domain, label = donor) %>%
        mutate(label = ifelse(parent == label, paste0('Mixed ', label), label)) %>%
        group_by(parent, label) %>%
        summarise(count = n(), prop = n() / fdat_len) %>%
        mutate(db = db, supergroups = sg, support = sist_supp)

      all_summaries <- rbind(all_summaries, fdat)
      # }
    }
  }
}

all_summaries %>%
  select(!prop) %>%
  filter(supergroups == 7, support == 'high') %>%
  pivot_wider(names_from = db, values_from = count) %>%
  arrange(-TOLDBA) %>%
  mutate(to_incl = (TOLDBA >= 20 | TOLDBB >= 20 | TOLDBC >= 20) & !str_detect(label, 'Mixed')) %>%
  filter(to_incl) %>%
  select(domain = parent, group = label) %>%
  write.table(file = paste0(odir, '/selected_groups.csv'), sep = '\t', quote = FALSE, row.names = FALSE)

a <- c('3sg - low',
       '3sg - moderate',
       '3sg - high',
       '5sg - low',
       '5sg - moderate',
       '5sg - high',
       '7sg - low',
       '7sg - moderate',
       '7sg - high',
       '9sg - low',
       '9sg - moderate',
       '9sg - high')


alph_asg <- all_summaries %>% mutate(x = paste0(supergroups, 'sg - ', support)) %>%
  mutate(x = factor(x, levels = a)) %>%
  filter(label %in% c('Alphaproteobacteria', 'Asgardarchaeota')) %>% 
  ggplot(aes(x, prop, colour = label, group = label)) +
  geom_point() +
  geom_line() +
  facet_grid(~db) +
  labs(title = 'Alpha vs Asgard') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  geom_vline(xintercept = '7sg - high', lty = 4) +
  theme(legend.position = 'bottom') +
  xlab('Stringency level') +
  ylab('Proportion')

# pdf('../paper/figures_R/raw_plots/figure_2/alph_asg.pdf', width = 7, height = 3.5)
alph_asg
# dev.off()

all_summaries %>% mutate(x = paste0(supergroups, 'sg - ', support)) %>%
  mutate(x = factor(x, levels = a)) %>%
  filter(label %in% c('Alphaproteobacteria', 'Gammaproteobacteria')) %>%
  ggplot(aes(x, prop, colour = label, group = label)) +
  geom_point() +
  geom_line() +
  facet_grid(~db) +
  labs(title = 'Alpha vs Gamma') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

all_summaries %>% mutate(x = paste0(supergroups, 'sg - ', support)) %>%
  mutate(x = factor(x, levels = a)) %>%
  filter(label %in% c('Alphaproteobacteria', 'Nucleocytoviricota')) %>%
  ggplot(aes(x, prop, colour = label, group = label)) +
  geom_point() +
  geom_line() +
  facet_grid(~db) +
  labs(title = 'Alpha vs Nucleocytoviricota') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
