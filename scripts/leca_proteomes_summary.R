library(dplyr)

dat <- c()
for (db in c('TOLDBA', 'TOLDBB', 'TOLDBC')) {
  a <- read.csv(paste0('../outputs/LECA_proteomes/', db, '_3sg_proteome.tsv'), sep = '\t')
  dat <- rbind(dat, a)
}

# 3 SG
dat %>%
  group_by(db) %>%
  count()

dat %>%
  group_by(origin, db) %>%
  count()

# 5 SG
dat %>%
  filter(supergroups >= 5) %>%
  group_by(db) %>%
  count()

dat %>%
  filter(supergroups >= 5) %>%
  group_by(origin, db) %>%
  count()
