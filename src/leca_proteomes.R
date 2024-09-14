# LECA proteomes
# Moisès Bernabeu
# Barcelona, March 2024

library(tidyverse)

odir = '../outputs/LECA_proteomes/'
metabolismdir = '../outputs/metabolism/'

dir.create(odir, showWarnings = FALSE, recursive = TRUE)
dir.create(metabolismdir, showWarnings = FALSE, recursive = TRUE)

to_exclude <- read.csv('../outputs/virus_analysis/excluded.tsv', sep = '\t')

for (db in c('TOLDBA', 'TOLDBB', 'TOLDBC')) {
  LECAOGs_sg <- read.csv(paste0('../outputs/innovations/', db, '_supergroups.tsv'), sep = '\t', col.names = c('Orthogroup', 'supergroups'))
  LECAOGs <- read.csv(paste0('../outputs/functional_annotation/', db, '_LECAOGs_annotation.tsv'), sep = '\t')
  LECACOG <- read.csv(paste0('../outputs/functional_annotation/', db, '_LECAOGs_annotation_COG.tsv'), sep = '\t')

  mLECAOGs <- read.csv(paste0('../outputs/functional_annotation/', db, '_mLECAOGs_annotation.tsv'), sep = '\t')
  mLECAsis <- read.csv(paste0('../outputs/functional_annotation/', db, '_mLECAOGs_sister_annotation.tsv'), sep = '\t')
  mLECACOG <- read.csv(paste0('../outputs/functional_annotation/', db, '_mLECAOGs_annotation_COG.tsv'), sep = '\t')
  mLECAori <- read.csv(paste0('../outputs/mLECA_origins/with_proportions/', db, '_stats.tsv'), sep = '\t')[, -c(24, 25)]

  largeOGs <- read.csv(paste0('../outputs/functional_annotation/', db, '_mLECA_large_KOs.tsv'), sep = '\t')
  largeOGssis <- read.csv(paste0('../outputs/functional_annotation/', db, '_mLECA_large_sister_KOs.tsv'), sep = '\t')
  largeOGsCOG <- read.csv(paste0('../outputs/functional_annotation/', db, '_mLECA_large_COGs.tsv'), sep = '\t')
  largeOGs_ori <- read.csv(paste0('../outputs/mLECA_origins/', db, '_large_stats.tsv'), sep = '\t')[, -c(24, 25)]

  total_OGs <- readLines(paste0('../outputs/innovations/01_', db, '_LECA_OGs.txt'))
  expan_OGs <- readLines(paste0('../outputs/innovations/03_', db, '_LECA_OGs_expanded.txt'))
  expan_large <- unique(str_split(largeOGs_ori$Orthogroup, '_', simplify = TRUE)[, 1])

  innovations <- total_OGs[!(total_OGs %in% c(expan_OGs, expan_large))]

  db_exclude <- to_exclude %>%
    filter(database == db)
  db_exclude <- db_exclude[, 2]

  # Gathering the information of the innovations
  inn_proteome <- LECAOGs %>%
    select(Orthogroup, LECA_KOs = KOs) %>%
    filter(Orthogroup %in% innovations) %>%
    left_join(LECACOG %>% select(Orthogroup, LECA_COG = KOs), by = 'Orthogroup') %>%
    left_join(LECAOGs_sg, by = 'Orthogroup') %>%
    mutate(origin = 'Innovation', sister_KOs = NA, coincide = NA, donor = NA, donor_domain = NA, stems = 2, db = db)

  # Getting the information of the large OGs
  large_aqc <- largeOGs_ori %>%
    mutate(Orthogroup = paste0(Orthogroup, '_', LECA), donor = gsub('_[A-Z]', '', donor)) %>%
    left_join(largeOGs %>% select(Orthogroup, LECA_KOs = KOs), by = 'Orthogroup') %>%
    left_join(largeOGssis %>% select(Orthogroup, sister_KOs = KOs), by = 'Orthogroup') %>%
    mutate(coincide = str_detect(LECA_KOs, sister_KOs) | str_detect(sister_KOs, LECA_KOs)) %>%
    left_join(largeOGsCOG %>% select(Orthogroup, LECA_COG = KOs), by = 'Orthogroup') %>%
    mutate(db = db, origin = 'Acquisition') %>%
    rename(donor_domain = S1_MAB_group_domain) %>%
    filter(!Orthogroup %in% db_exclude)

  # Gathering the information of the acquisitions
  acq_proteome <- mLECAori %>%
    mutate(Orthogroup = paste0(Orthogroup, '_', LECA), donor = gsub('_[A-Z]', '', donor)) %>%
    left_join(mLECAOGs %>% select(Orthogroup, LECA_KOs = KOs), by = 'Orthogroup') %>%
    left_join(mLECAsis %>% select(Orthogroup, sister_KOs = KOs), by = 'Orthogroup') %>%
    mutate(coincide = str_detect(LECA_KOs, sister_KOs) | str_detect(sister_KOs, LECA_KOs)) %>%
    left_join(mLECACOG %>% select(Orthogroup, LECA_COG = KOs), by = 'Orthogroup') %>%
    mutate(db = db, origin = 'Acquisition') %>%
    rename(donor_domain = S1_MAB_group_domain) %>%
    filter(!Orthogroup %in% db_exclude)

  acq_proteome <- rbind(large_aqc, acq_proteome)

  write.table(acq_proteome, file = paste0(odir, '/', db, '_acq_proteome.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)

  # Merging into a summary table of the proteome for the three supergroups
  proteome <- rbind(acq_proteome %>% select(Orthogroup, LECA_KOs, sister_KOs, coincide,
                                            LECA_COG, supergroups, origin, stems, donor,
                                            donor_domain, db),
                    inn_proteome) %>% filter(supergroups >= 3)
  
  for (i in c(3, 5)) {
    sub_proteome <- proteome %>%
      filter(supergroups >= i)
    write.table(sub_proteome, file = paste0(odir, '/', db, '_', i, 'sg_proteome.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)

    kos <- proteome %>%
      filter(supergroups >= i) %>%
      mutate(source = 'KOfam') %>%
      select(gene_id = Orthogroup, enzyme_accession = LECA_KOs, source) %>%
      separate_longer_delim(enzyme_accession, ';')
    
    write.table(kos, file = paste0(metabolismdir, '/', db, '_', i, 'sg_KOori_proteome.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
  }
}

