# LECA sunburst
# Moisès Bernabeu
# Cocentaina, August 2024

library(tidyverse)
library(plotly)

theme_set(theme_bw())

db = 'TOLDBA'
sg = 3

odir <- '../outputs/mLECA_origins_plots'
dir.create(odir, showWarnings = FALSE, recursive = TRUE)

for (db in c('TOLDBA', 'TOLDBB', 'TOLDBC')) {
  for (sg in c(3, 5)) {
    dat <- read.csv(paste0('../outputs/LECA_proteomes/', db, '_', sg, 'sg_proteome.tsv'), sep = '\t')
    
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

