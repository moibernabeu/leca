library(tidyverse)
theme_set(theme_bw())

files <- list.files('../tests/', pattern = 'orthogroups_summary.tsv', full.names = TRUE, recursive = TRUE)

odf <- c()
alldf <- c()
for (file in files) {
  ogs_summ <- read.csv(file, sep = '\t')

  info <- c(str_split(str_split(file, '/', simplify = TRUE)[, 4], '_', simplify = TRUE))

  alldf <- rbind(alldf, ogs_summ %>%
                   mutate(is_leca = (Two_stems == 'true' & Five_sp == 'true' & Three_supergroups == 'true'),
                          database = info[1], inflation = info[2]))

  ogs_summ <- ogs_summ %>%
    mutate(is_leca = (Two_stems == 'true' & Five_sp == 'true' & Three_supergroups == 'true'))

  osdf <- ogs_summ %>%
    mutate(groups = cut(supergroups, breaks = c(3, 5, 5.1, 10), labels = c('3 supergroups', '5 supergroups', 'More'), right = FALSE)) %>%
    group_by(is_leca, groups) %>%
    summarise(LECA_OGs = n(), mean_seqs_no = mean(seqs_no), max_seqs_no = max(seqs_no), min_seqs_no = min(seqs_no)) %>%
    filter(is_leca) %>%
    gather(value = value, key = key, -groups, -is_leca) %>%
    mutate(database = info[1], inflation = info[2])
  
  odf <- rbind(odf, osdf)
}

write.table(odf, file = '../tests/merged_ogs_summary.tsv', sep = '\t', quote = FALSE)

alldf %>%
  group_by(is_leca, database, inflation) %>%
  mutate(is_leca = ifelse(is_leca, 'LECA', 'Not LECA')) %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = is_leca, values_from = n) %>%
  write.table(sep = '\t', quote = FALSE, row.names = FALSE)

alldf %>%
  filter(is_leca) %>%
  group_by(database, inflation, Two_stems, Three_supergroups, Five_supergroups, Seven_supergroups) %>%
  summarise(n = n()) %>%
  gather(key = key, value = value, -database, -inflation, -n) %>%
  filter(value == 'true') %>%
  group_by(database, inflation, key) %>%
  summarise(n = sum(n)) %>%
  pivot_wider(names_from = key, values_from = n) %>%
  write.table(sep = '\t', quote = FALSE, row.names = FALSE)

ggplot(odf %>% filter(key != 'min_seqs_no'), aes(groups, value, fill = inflation)) +
  geom_col(position = position_dodge(width = 1)) +
  geom_text(position = position_dodge2(width = 1), size = 3, aes(label = value)) +
  facet_grid(key~database, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank())

alldf %>%
  mutate(n_gth_500 = seqs_no >500) %>%
  group_by(database, inflation, n_gth_500) %>%
  summarise(n = n()) %>%
  
