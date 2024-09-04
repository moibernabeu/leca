# Assessment of the NR first hits
# Moisès Bernabeu
# Barcelona, July 2024

library(tidyverse)

gnms <- read.csv('../data/gnms_table.tsv', sep = '\t', header = FALSE) %>%
  rename(rename = V1, mnemo = V2)
initial_stats <- read.csv('../outputs/renamed_proteomes/proteomes_stats.tsv', sep = '\t') %>%
  select(!c(format, type)) %>%
  mutate(mnemo = gsub('\\..*', '', file)) %>%
  select(!file) %>%
  set_names(c(paste0('raw_', c('num_seqs', 'sum_len', 'min_len', 'avg_len', 'max_len')), 'mnemo'))
complexity <- read.csv('../outputs/cleaned_proteomes/complexity_analysis.tsv', sep = '\t') %>%
  select(mnemo = proteome, lc_index, raw_lc_prots_prop, masked_lc_prots_prop)
final_stats <- read.csv('../outputs/clustered_proteomes/clustered_proteomes_stats.tsv', sep = '\t') %>%
  select(!c(format, type)) %>%
  mutate(mnemo = gsub('_.*', '', file)) %>%
  select(!file) %>%
  set_names(c(paste0('final_', c('num_seqs', 'sum_len', 'min_len', 'avg_len', 'max_len')), 'mnemo'))
initial_busco <- read.csv('../outputs/busco_raw/busco_summary.tsv', sep = '\t', header = FALSE) %>%
  mutate(mnemo = gsub('\\..*', '', gsub('.*/', '', V1))) %>%
  select(!V1) %>%
  set_names(c(paste0('raw_', c('completeness', 'multi_copy', 'fragmented', 'missing')), 'mnemo'))
final_busco <- read.csv('../outputs/busco_final/busco_summary.tsv', sep = '\t', header = FALSE) %>%
  mutate(mnemo = gsub('\\..*', '', gsub('.*/', '', V1))) %>%
  select(!V1) %>%
  set_names(c(paste0('final_', c('completeness', 'multi_copy', 'fragmented', 'missing')), 'mnemo'))

summary_table <- gnms %>%
  left_join(initial_stats) %>%
  left_join(final_stats) %>%
  left_join(initial_busco) %>%
  left_join(final_busco) %>%
  left_join(complexity) %>%
  mutate(len_diff = raw_num_seqs - final_num_seqs,
         len_diff_prop = len_diff / raw_num_seqs * 100,
         compl_diff = raw_completeness - final_completeness)

write.table(summary_table, file = '../outputs/clustered_proteomes/final_proteomes_summary.tsv',
            sep = '\t', row.names = FALSE, quote = FALSE)  
