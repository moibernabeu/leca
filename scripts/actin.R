# Actin trees

defs <- read.csv('../../../dbs/profiles/kos_annotation.tsv', sep = '\t')
actin_KOs <- c('K10354', 'K12313', 'K12314', 'K12315', 'K05692', 'K10355')

proteomes <- c()
for (i in c('A', 'B', 'C')) {
  proteomes <- rbind(proteomes, read.csv(paste0('../outputs/LECA_proteomes_vircleaned//TOLDB', i, '_3sg_proteome.tsv'), sep = '\t'))
}

proteomes %>%
  separate_longer_delim(LECA_KOs, ';') %>%
  separate_longer_delim(sister_KOs, ';') %>%
  select(Orthogroup, db, LECA_KOs, sister_KOs, origin, donor, supergroups) %>%
  filter(LECA_KOs %in% actin_KOs, sister_KOs %in% actin_KOs)


