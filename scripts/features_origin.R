# Features origin
# Moisès Bernabeu
# Barcelona, September 2024

library(tidyverse)

theme_set(theme_bw())

cols <- c(Archaea = '#7179E6',
          Bacteria = '#E62E46',
          Innovation = '#E67E0B',
          Viruses = '#2E8B32',
          Other = 'grey70')


dat <- read.csv('../outputs/metabolism_vircleaned/donor_KEGG_donors.tsv', sep = '\t')

features_summary <- dat %>%
  group_by(pathway) %>%
  mutate(total_KOs = n()) %>%
  group_by(type, subtype, pathway, donor) %>%
  summarise(donor_prop = n() / unique(total_KOs)) %>%
  arrange(-donor_prop, .by_group = TRUE)

write.table(features_summary, file = '../outputs/metabolism_vircleaned/donor_KEGG_donors_summary.tsv', sep = '\t', quote = FALSE, row.names = FALSE)

features_summary_first <- dat %>%
  group_by(pathway) %>%
  mutate(total_KOs = n()) %>%
  group_by(type, subtype, pathway, donor) %>%
  summarise(donor_prop = n() / unique(total_KOs)) %>%
  arrange(-donor_prop, .by_group = TRUE) %>%
  slice_max(n = 2, order_by = donor_prop)

write.table(features_summary_first, file = '../outputs/metabolism_vircleaned/donor_KEGG_donors_summary_first2.tsv', sep = '\t', quote = FALSE, row.names = FALSE)

dat <- read.csv('../outputs/metabolism_vircleaned/donor_domain_KEGG_donors.tsv', sep = '\t')

features_summary <- dat %>%
  group_by(pathway) %>%
  mutate(total_KOs = n()) %>%
  group_by(type, subtype, pathway, donor) %>%
  summarise(donor_domain_prop = n() / unique(total_KOs)) %>%
  arrange(-donor_domain_prop, .by_group = TRUE)

write.table(features_summary, file = '../outputs/metabolism_vircleaned/donor_domain_KEGG_donors_summary.tsv', sep = '\t', quote = FALSE, row.names = FALSE)

features_summary_first <- dat %>%
  group_by(pathway) %>%
  mutate(total_KOs = n()) %>%
  group_by(type, subtype, pathway, donor) %>%
  summarise(donor_domain_prop = n() / unique(total_KOs)) %>%
  arrange(-donor_domain_prop, .by_group = TRUE) %>%
  slice_max(n = 2, order_by = donor_domain_prop)

write.table(features_summary_first, file = '../outputs/metabolism_vircleaned/donor_domain_KEGG_donors_summary_first2.tsv', sep = '\t', quote = FALSE, row.names = FALSE)

# dat %>%
#   select(type, subtype, pathway) %>%
#   unique() %>%
#   mutate(nKOs = gsub('\\)', '', str_split(pathway, ' \\(', simplify = TRUE)[, 2])) %>%
#   write.table(file = '../outputs/metabolism_vircleaned/donor_domain_paths.tsv',
#               sep = '\t', quote = FALSE, row.names = FALSE)

paths <- read.csv('../outputs/metabolism_vircleaned/donor_domain_paths.tsv', sep = '\t') %>%
  filter(In.LECA == 'X')
paths <- paths$pathway

library(see)
ggplot(features_summary %>% filter(type == 'Cellular Processes', pathway %in% paths), aes(donor_domain_prop, pathway, fill = donor)) +
  geom_col() +
  scale_fill_okabeito() +
  labs(title = 'Cellular Processes')

dat <- features_summary %>%
  filter(!(donor == 'Other' & donor_domain_prop >= 0.15), !(donor %in% c('Other', 'Viruses'))) %>%
  group_by(type, subtype, pathway) %>%
  mutate(new_total = sum(donor_domain_prop)) %>%
  filter(new_total >= 0.85) %>%
  mutate(new_prop = donor_domain_prop / new_total) %>%
  select(!c(new_total, donor_domain_prop)) %>%
  filter(pathway %in% paths) %>%
  pivot_wider(names_from = donor, values_from = new_prop, values_fill = 0)

library(ggtern)
ggtern(dat, aes(Innovation, Bacteria, Archaea, colour = type)) +
  geom_point() +
  theme_bw() +
  geom_Tline(Tintercept = 0.5, lty = 4, col = 'grey60') +
  geom_Rline(Rintercept = 0.2, lty = 4, col = 'grey60')

library(ggfortify)
boxplot(dat[, c(4:8)])
a <- princomp(dat[, c(4:8)], cor = TRUE)
autoplot(a, loadings = TRUE, loadings.label = TRUE, data = dat, colour = 'type')

boxplot(dat[, c(4:8)])
a <- princomp(dat[, c(4, 6:8)])
plot(a)
autoplot(a, loadings = TRUE, loadings.label = TRUE, data = dat, colour = 'type')


