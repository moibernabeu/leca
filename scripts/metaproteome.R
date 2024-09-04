library(tidyverse)

dbs <- c('TOLDBA', 'TOLDBB', 'TOLDBC')

KOs <- c()
COGs <- c()
origins <- c()
for (db in dbs) {
  ddat <- read.csv(paste0('../outputs/LECA_proteomes_vircleaned/', db, '_3sg_proteome.tsv'), sep = '\t')
  dKOs <- ddat %>%
    select(KO = LECA_KOs, db, supergroups) %>%
    separate_longer_delim(KO, ';')
  KOs <- rbind(KOs, dKOs)
  
  dCOGs <- ddat %>%
    select(COG = LECA_COG, db, supergroups) %>%
    separate_longer_delim(COG, ';')
  COGs <- rbind(COGs, dCOGs)
  
  dorigins <- ddat %>%
    select(KO = LECA_KOs, origin, donor, donor_domain, db, supergroups) %>%
    separate_longer_delim(KO, ';')
  origins <- rbind(origins, dorigins)
}

KOs <- KOs %>%
  group_by(db, KO) %>%
  summarise(supergroups = max(supergroups)) %>%
  ungroup() %>%
  group_by(KO) %>%
  summarise(max_supergroups = max(supergroups), pervasiveness = n(),
            to_keep = max_supergroups >= 5 | pervasiveness >= 2) %>%
  na.omit()

to_keep <- KOs %>%
  filter(to_keep)

metaprot <- data.frame(gene_id = paste('LECA_', 1:length(to_keep$KO), sep = ''),
                       enzyme_accession = to_keep$KO, source = 'KOfam')

# write.table(metaprot, file = '../outputs/LECA_proteomes_vircleaned/metaproteome.tsv', sep = '\t', quote = FALSE, row.names = FALSE)

COGs <- COGs %>%
  group_by(db, COG) %>%
  summarise(supergroups = max(supergroups)) %>%
  ungroup() %>%
  group_by(COG) %>%
  summarise(max_supergroups = max(supergroups), pervasiveness = n(),
            to_keep = max_supergroups >= 5 | pervasiveness >= 2) %>%
  na.omit()

to_keep <- COGs %>%
  filter(to_keep)

metaprot_COGs <- data.frame(gene_id = paste('LECA_', 1:length(to_keep$COG), sep = ''),
                            enzyme_accession = to_keep$COG)

# write.table(metaprot_COGs, file = '../outputs/LECA_proteomes_vircleaned/metaproteome_COGs.tsv', sep = '\t', quote = FALSE)

main_donors <- read.csv('../outputs/mLECA_origins_vircleaned/selected_groups.csv', sep = '\t')[, 2]
defs <- read.csv('../data/KO_description.tsv', sep = '\t')

KO_origins <- origins %>%
  filter(KO %in% metaprot$enzyme_accession) %>%
  mutate(donor = ifelse(origin == 'Innovation', origin,
                        ifelse(donor %in% main_donors, donor, 'Other'))) %>%
  select(!donor_domain) %>%
  group_by(KO) %>%
  mutate(total = n()) %>%
  group_by(KO, origin, donor, total) %>%
  summarise(n_donor = n(), prop_donor = n_donor / unique(total)) %>%
  filter(!is.na(KO)) %>%
  left_join(defs, by = c('KO' = 'ko'))

write.table(KO_origins, '../outputs/LECA_proteomes_vircleaned/metaproteome_KOs_origin.tsv', sep = '\t', quote = FALSE, row.names = FALSE)

cat('', file = '../outputs/metabolism_vircleaned/donor_KEGG_mapper.tsv', append = F)
for (i in unique(KO_origins$donor)) {
  cat(paste0('# ', i, '\n'), file = '../outputs/metabolism_vircleaned/donor_KEGG_mapper.tsv', append = T)
  donKOs <- KO_origins %>%
    filter(donor == i) %>%
    na.omit()
  donKOs <- unique(donKOs$KO)
  cat(paste(paste(i, 1:length(donKOs), sep = '_'), donKOs, sep = "\t"), sep = '\n',
      file = '../outputs/metabolism_vircleaned/donor_KEGG_mapper.tsv', append = T)
}

KO_domains <- origins %>%
  filter(KO %in% metaprot$enzyme_accession) %>%
  mutate(donor_domain = ifelse(origin == 'Innovation', 'Innovation', ifelse(str_detect(donor_domain, ';'), 'Other', donor_domain))) %>%
  group_by(KO) %>%
  mutate(total = n()) %>%
  group_by(KO, origin, donor_domain, total) %>%
  summarise(n_donor = n(), prop_donor = n_donor / unique(total)) %>%
  filter(!is.na(KO)) %>%
  left_join(defs, by = c('KO' = 'ko'))

write.table(KO_origins, '../outputs/LECA_proteomes_vircleaned/metaproteome_KOs_origin_domain.tsv', sep = '\t', quote = FALSE, row.names = FALSE)

cat('', file = '../outputs/metabolism_vircleaned/donor_domain_KEGG_mapper.tsv', append = F)
for (i in unique(KO_domains$donor_domain)) {
  cat(paste0('# ', i, '\n'), file = '../outputs/metabolism_vircleaned/donor_domain_KEGG_mapper.tsv', append = T)
  donKOs <- KO_domains %>%
    filter(donor_domain == i) %>%
    na.omit()
  donKOs <- unique(donKOs$KO)
  cat(paste(paste(i, 1:length(donKOs), sep = '_'), donKOs, sep = "\t"), sep = '\n',
      file = '../outputs/metabolism_vircleaned/donor_domain_KEGG_mapper.tsv', append = T)
}
