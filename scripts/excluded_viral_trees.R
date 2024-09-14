# Get a list of families to exclude
# Moisès Bernabeu
# Barcelona, September 2024

library(tidyverse)

dat <- c()
for (db in c('TOLDBA', 'TOLDBB', 'TOLDBC')) {
  a <- read.csv(paste0('../outputs/virus_analysis/', db, '.tsv'), sep = '\t')
  dat <- rbind(dat, data.frame(a, database = db))
}

excluded <- dat %>%
  mutate(Orthogroup = paste0(og, '_', LECA)) %>%
  select(database, Orthogroup, supergroups, ntips, leca_tips, sister_tips, viruses,
         nvirseqs, nv_sister_domain, nv_sister, nv_sister_size) %>%
  filter(nv_sister_domain == 'Eukaryota' | nv_sister_domain == '')

excluded %>%
  select(database, Orthogroup) %>%
  write.table(file = '../outputs/virus_analysis/excluded.tsv', sep = '\t', quote = FALSE, row.names = FALSE)
