# Stem-length analysis

library(ggplot2)
library(stringr)
library(dplyr)
library(ggridges)
library(tidytext)
library(tidyr)
library(forcats)
library(see)

theme_set(theme_bw())

dist_mode <- function(x) {
  if (length(x) > 1) {
    d <- density(x)
    return(d$x[which.max(d$y)])
  } else {
    return(x)
  }
}

find_module <- function(x) {
  if (length(x) > 2) {
    dens <- density(x)
    densx <- dens$x
    densy <- dens$y
    maxmin <- c(FALSE, diff(diff(dens$y) > 0) != 0, FALSE)
    infl <- c(FALSE, FALSE, diff(diff(diff(dens$y)) > 0) != 0, FALSE)
    
    first_min <- which(dens$x > dist_mode(x) & maxmin)[1]
    second_infl <- which(dens$x > dist_mode(x) & infl)[2]
    
    # plot(dens)
    # points(densx[maxmin], densy[maxmin], pch = 19, col = 'steelblue')
    # points(densx[infl], densy[infl], pch = 19, col = 'darkorange3')
    # abline(v = find_module(x))
    
    # When the inflection point is lower than the first minimum, get the lower value  
    # if (is.na(densx[first_min]) | densx[first_min] > densx[second_infl]) {
    #   return(densx[second_infl])
    # } else {
    #   return(densx[first_min])
    # }
    
    if (is.na(densx[first_min])) {
      return(densx[second_infl])
    } else {
      return(densx[first_min])
    }
  } else {
    return(NA) 
  }
}

to_exclude <- read.csv('../outputs/virus_analysis/excluded.tsv', sep = '\t')

odir <- '../outputs/modules/'
db <- 'TOLDBC'
prefix <- paste0(odir, db)
dir.create(odir)

dat <- read.csv(paste0('../outputs/mLECA_origins/', db, '_stats.tsv'), sep = '\t',
                tryLogical = TRUE)
groups <- read.csv('../outputs/LECA_proteomes/selected_groups.tsv', sep = '\t')[, 2]

db_exclude <- to_exclude %>%
  filter(database == db)
db_exclude <- db_exclude[, 2]

# Plotting the stem lengths and calculating the modules ----
# General filtering for the data
fdat <- dat %>%
  filter(donor %in% groups,
         sl <= quantile(sl, 0.99),
         !paste0(Orthogroup, '_', LECA) %in% db_exclude)

pdf(paste0(prefix, '_stem_lengths_modules.pdf'), width = 8.5, height = 6)
criteria <- c(3, 5)
for (criterion in criteria) {
  pdat <- fdat %>%
    filter(supergroups >= criterion)

  pdat <- pdat %>%
    group_by(donor) %>%
    mutate(label = paste0(donor, ' (', n(), ')'))
  
  modes_sort <- pdat %>%
    group_by(label) %>%
    summarise(mode = dist_mode(sl), first_min = find_module(sl)) %>%
    arrange(-mode)
  
  pdat$label <- factor(pdat$label, levels = modes_sort$label)
  modes_sort$label <- factor(modes_sort$label, levels = modes_sort$label)

  p <- ggplot(pdat, aes(sl)) +
    geom_density() +
    geom_vline(aes(xintercept = mode), data = modes_sort, colour = 'steelblue') +
    geom_vline(aes(xintercept = first_min), data = modes_sort, colour = 'darkorange3', lty = 4) +
    facet_wrap(~label) +
    labs(title = paste0(db, ' | ', criterion, ' supergroups'))
  
  print(p)
}
dev.off()

# Plotting the number of trees per dataset ----
pdf(paste0(prefix, '_no_trees.pdf'), width = 5, height = 7)
fdat %>%
  gather(Three_supergroups, Five_supergroups, key = criterion, value = has_crit) %>%
  filter(has_crit == 'True') %>%
  mutate(database = db) %>%
  group_by(donor, database, criterion) %>%
  summarise(trees_no = n()) %>%
  ggplot(aes(trees_no, reorder(donor, trees_no))) +
  geom_point() +
  geom_segment(aes(x = 0, xend = trees_no, y = donor, yend = donor)) +
  facet_grid(criterion ~ database, scales = 'free_y') +
  xlab('Number of trees') +
  ylab('Clade')
dev.off()

# Extracting modules ----
modules <- fdat %>%
  gather(Three_supergroups, Five_supergroups, key = criterion, value = has_crit) %>%
  mutate(database = db) %>%
  filter(has_crit == 'True') %>%
  group_by(database, donor, criterion) %>%
  mutate(mode = dist_mode(sl), module_threshold = find_module(sl), in_module = sl <= module_threshold) %>%
  select(Orthogroup, LECA, LECA_size, database, criterion, donor, sl, module_threshold, in_module, leca_support, leca_sister_support)

modules_summary <- modules %>%
  group_by(database, criterion, donor) %>%
  summarise(ntrees = n(), trees_in_module = sum(in_module),
            prop_in = trees_in_module / ntrees,
            prop_out = 1 - prop_in)

pdf(paste0(prefix, '_rm_trees_prop.pdf'), width = 5, height = 3.25)
ggplot(modules_summary, aes(prop_out)) +
  geom_density() +
  xlab('Proportion of removed trees') +
  ylab('Density')
dev.off()

# Plotting module trees ----
pdf(paste0(prefix, '_trees_in_modules.pdf'), width = 11, height = 7)
criteria <- c('Three_supergroups', 'Five_supergroups')
for (crit in criteria) {
  pdat <- modules %>%
    filter(criterion == crit) %>%
    group_by(donor) %>%
    mutate(ntrees = n(), mode = dist_mode(sl),  in_trees = sum(in_module),
           label = paste0(donor, ' (', ntrees, '/', in_trees, ')'))
  
  modes_sort <- pdat %>%
    group_by(label) %>%
    summarise(mode = dist_mode(sl), first_min = find_module(sl)) %>%
    arrange(-mode)
  
  pdat$label <- factor(pdat$label, levels = modes_sort$label)
  modes_sort$label <- factor(modes_sort$label, levels = modes_sort$label)

  p <- ggplot(pdat, aes(sl)) +
    geom_density() +
    geom_point(aes(y = 0, colour = in_module), alpha = 0.5) +
    geom_vline(aes(xintercept = mode), lty = 1, colour = 'steelblue') +
    geom_vline(aes(xintercept = module_threshold), lty = 4, colour = 'darkorange3') +
    facet_wrap(~label) +
    labs(title = paste(db, crit, sep = ' | ')) +
    xlab('Stem length') +
    ylab('Density')
  
  print(p)
}
dev.off()

# Annotating modules
annotation <- read.csv(paste0('../outputs/functional_annotation/', db, '_mLECAOGs_annotation.tsv'), sep = '\t')[, c(1, 2)]
kos_definition <- read.table('../../../dbs/profiles/kos_annotation.tsv', sep = '\t', quote = "", 
                             header = 1,
                             stringsAsFactors = FALSE)

modules_df <- modules %>%
  mutate(Orthogroup = paste0(Orthogroup, '_', LECA)) %>%
  left_join(annotation %>% select(Orthogroup, KOs), by = c('Orthogroup' = 'Orthogroup')) %>%
  separate_longer_delim(KOs, ';') %>%
  left_join(kos_definition %>% select(knum, definition), by = c('KOs' = 'knum'))

write.table(modules_df, file = paste0(prefix, '_modules_annotation.tsv'), sep = '\t', row.names = FALSE, quote = FALSE)

modules_df_filt <- modules_df %>%
  filter(!(is.na(KOs)))

write.table(modules_df_filt, file = paste0(prefix, '_modules_annotation_filtered.tsv'), sep = '\t', row.names = FALSE, quote = FALSE)

pdf(paste0(prefix, '_modules_annotation_level.pdf'), width = 8, height = 10)
modules_df %>%
  group_by(criterion, database, donor) %>%
  summarise(trees = n(), in_module = sum(in_module), annotated = sum(!(is.na(definition)))) %>%
  gather(key = 'cat', value = 'ntrees', -criterion, -database, -donor) %>%
  ggplot(aes(reorder(donor, ntrees), ntrees, colour = cat)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_linerange(aes(xmin = donor, xmax = donor, ymax = ntrees, ymin = 0), position = position_dodge(width = 0.5)) +
  coord_flip() +
  scale_colour_okabeito() +
  facet_grid(criterion ~ database) +
  xlab('Number of trees') +
  ylab('Class') +
  labs(colour = 'Set')
dev.off()

modules %>%
  filter(in_module, donor %in% groups) %>%
  ggplot(aes(sl, colour = donor)) +
  geom_density() +
  facet_wrap(~criterion)
