# ToLDB manufacturing
# Moisès Bernabeu
# Barcelona, July 2024

library(tidyverse)

stats <- read.csv('../data/repdb_meta.tsv', sep = '\t')

supergroups <- c('Chloroplastida' = 'Archaeplastida_plus',
                 'Glaucophyta' = 'Archaeplastida_plus',
                 'Picozoa' = 'Archaeplastida_plus',
                 'Rhodelphidia' = 'Archaeplastida_plus',
                 'Rhodophyta' = 'Archaeplastida_plus',
                 'Cryptophyceae' = 'Archaeplastida_plus',
                 'Telonemia' = 'TRASH',
                 'Haptophyta' = 'TRASH',
                 'Rhizaria' = 'TRASH',
                 'Stramenopiles' = 'TRASH',
                 'Alveolata' = 'TRASH',
                 'Meteora' = 'PHM',
                 'Hemimastigophora' = 'PHM',
                 'Provora' = 'PHM',
                 'Heterolobosea' = 'Discoba',
                 'Euglenozoa' = 'Discoba',
                 'Jakobida' = 'Discoba',
                 'Amoebozoa' = 'Amorphea',
                 'Breviatea' = 'Amorphea',
                 'Apusomonadida' = 'Amorphea',
                 'Opisthokonta' = 'Amorphea',
                 'Mantamonas' = 'CRuMs',
                 'Rigifila' = 'CRuMs',
                 'Collodictyon' = 'CRuMs',
                 'Malawimonadida' = 'Malawimonadida',
                 'Ancyromonadida' = 'Ancyromonadida',
                 'Preaxostyla' = 'Metamonada',
                 'Fornicata' = 'Metamonada',
                 'Parabasalia' = 'Metamonada',
                 'Anaeramoebidae' = 'Metamonada',
                 'Centroplasthelida' = 'TRASH',
                 'Ancoracysta' = 'PHM',
                 'Kathablepharidacea' = 'Archaeplastida_plus',
                 'Palpitomonas' = 'Archaeplastida_plus')

data <- stats %>%
  mutate(supergroup = ifelse(p %in% names(supergroups), supergroups[p],
                             ifelse(c %in% names(supergroups), supergroups[c],
                                    ifelse(p %in% supergroups, p,
                                           ifelse(c %in% supergroups, c, NA)))))

supergroups <- unique(data$supergroup)

supergroups[1]
a <- data %>%
  group_by(c) %>%
  mutate(mean_num_seqs = mean(num_seqs)) %>%
  slice_max(order_by = completeness, n = 15, with_ties = FALSE)

write.table(a, file = '../data/repdb_filtered_15best.tsv', sep = '\t', row.names = FALSE)
