library(tidyverse)

dat <- read.csv('../data/gnms_list.txt', sep = '\t')

repdb <- read.csv('../data/present_in_repdp.txt', sep = '\t', header = FALSE) %>%
  select(Proteome.ID = V1, RD_sp = V8)

ep <- read.csv('../data/ep_gnms_species.txt', sep = '\t', header = FALSE) %>%
  rename(c(Proteome.ID = V1, EP_sp = V2))
ncbi <- read.csv('../data/ncbi_gnms_species.txt', sep = '\t', header = FALSE) %>%
  mutate(V1 = gsub('\\..*', '', V1)) %>%
  rename(c(Proteome.ID = V1, NCBI_sp = V2))
up <- read.csv('../data/up_gnms_species.txt', sep = '\t', header = FALSE) %>%
  mutate(V1 = gsub('\\..*', '', V1)) %>%
  rename(c(Proteome.ID = V1, UP_sp = V2))
p10K <- read.csv('../data/p10k_meta.tsv', sep = '\t') %>%
  select(p10k_id, species) %>%
  rename(c(Proteome.ID = p10k_id, P10K_sp = species))

merged_dat <- dat %>%
  left_join(repdb) %>%
  left_join(ep) %>%
  left_join(ncbi) %>%
  left_join(up) %>%
  left_join(p10K)

sp_summary <- c()
for (i in 1:dim(merged_dat)[1]) {
  x <- unique(as.character(na.omit(as.character(merged_dat[i, 6:10]))))
  sp_summary <- rbind(sp_summary, data.frame(id = merged_dat[i, 2], sp = paste(x, collapse = ' | ')))
}

write.table(sp_summary, file = '../data/sp_check_summary.tsv', sep = '\t')
