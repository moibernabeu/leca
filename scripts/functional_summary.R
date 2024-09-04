# Functional summary
# Moisès Bernabeu
# Cocentaina, August 2024

library(tidyverse)

theme_set(theme_bw())

dat <- read.csv('../outputs/LECA_proteomes/TOLDBA_3sg_proteome.tsv', sep = '\t')
cog_defs <- read.csv('../data/cog-20.def.tab', sep = '\t', header = FALSE)
cog_desc <- read.csv('../data/fun-20.tab', sep = '\t', header = FALSE)

cog_defs <- cog_defs %>%
  mutate(category = lapply(str_split(V2, ''), paste0, collapse = ';')) %>%
  select(COG = V1, category)

dat <- dat %>%
  mutate(donor_domain = ifelse(origin == 'Innovation', 'Innovation', donor_domain),
         donor = ifelse(origin == 'Innovation', 'Innovation', donor),
         donor_domain = ifelse(str_detect(donor_domain, ';'), 'Undetermined', donor_domain)) %>%
  select(Orthogroup, LECA_COG, origin, donor, donor_domain) %>%
  separate_longer_delim(LECA_COG, ';') %>%
  left_join(cog_defs, by = c('LECA_COG' = 'COG')) %>%
  separate_longer_delim(category, ';') %>%
  left_join(cog_desc %>% select(category = V1, description = V3), by = 'category') %>%
  na.omit()

table(dat$origin)
table(dat$donor_domain)

background <- dat %>%
  group_by(category, description) %>%
  summarise(background = n()) %>%
  ungroup() %>%
  mutate(background_prop = background / sum(background))

origin <- dat %>%
  group_by(category, description, donor_domain) %>%
  summarise(origin_n = n()) %>%
  ungroup() %>%
  group_by(donor_domain) %>%
  mutate(origin_prop = origin_n / sum(origin_n)) %>%
  left_join(background, by = c('category', 'description')) %>%
  mutate(ratio = origin_prop / background_prop,
         overexpressed = ratio > 1)

origin$donor_domain <- factor(origin$donor_domain, levels = c('Archaea', 'Bacteria', 'Viruses', 'Undetermined', 'Innovation'))

p <- ggplot(origin, aes(ratio, reorder(description, ratio), colour = overexpressed)) +
  geom_point() +
  geom_vline(xintercept = 1) +
  facet_wrap(~donor_domain, nrow = 1) +
  xlab('Ratio') +
  ylab('COG category') +
  theme(legend.position = 'bottom')

pdf('../outputs/functional_annotation/functional_summary.pdf', width = 10, height = 4.5)
p
dev.off()

