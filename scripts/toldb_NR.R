# Assessment of the NR first hits
# Moisès Bernabeu
# Barcelona, July 2024

library(tidyverse)

blast <- read.csv('../outputs/gnms_affiliation/gnms_NR.blast',
                  sep = '\t', header = FALSE)

# species <- blast %>%
#   mutate(V16 = V15, mnemo = gsub('_.*', '', V1)) %>%
#   separate_wider_delim(cols = V16, delim = '[', names_sep = '_', too_few = 'align_end') %>%
#   mutate(V16 = gsub(']', '', V16_3)) %>%
#   select(!c(V16_1, V16_2, V16_3))
# 
# cat(unique(species$V16), sep = '\n', file = '../tests/species.txt')
# # Then taxonkit

taxonomy <- read.csv('../tests/ncbi_lineage.txt', sep = '\t', header = FALSE) %>%
  filter(V3 != '') %>%
  mutate(V3 = gsub('[a-z]__', '', V3)) %>%
  separate_wider_delim(V3, ';', names = c('d', 'p', 'c', 'o', 'f', 'g', 's'))

toldb <- read.csv('../tests/preliminar_lineage.txt', sep = '\t', header = FALSE)
names(toldb) <- c('mnemo', paste('R_', names(toldb)[2:dim(toldb)[2]], sep = ''))


blast_tax <- blast %>%
  mutate(V16 = V15, mnemo = gsub('_.*', '', V1)) %>%
  separate_wider_delim(cols = V16, delim = '[', names_sep = '_', too_few = 'align_end') %>%
  mutate(V16 = gsub(']', '', V16_3)) %>%
  select(!c(V16_1, V16_2, V16_3)) %>%
  left_join(taxonomy, by = c('V16' = 'V1')) %>%
  left_join(toldb, by = 'mnemo')

blast_tax_fh <- blast_tax %>%
  group_by(V1) %>%
  slice_max(order_by = V12, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(genus = gsub(' .*', '', R_V4), genus_coincidence = g == genus)

write.table(blast_tax_fh, '../outputs/gnms_affiliation/first_hits.tsv', sep = '\t', row.names = FALSE)

fh_summ <- blast_tax_fh %>%
  group_by(mnemo) %>%
  summarise(coincidence_prop = sum(genus_coincidence) / n()) %>%
  arrange(-coincidence_prop) %>%
  print(n=1000)

