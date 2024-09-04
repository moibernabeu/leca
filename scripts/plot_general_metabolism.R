library(ggkegg)
library(igraph)
library(tidyverse)
library(dplyr)
library(rjson)

palette_df <- read.csv('../data/colours.tsv', sep = '\t', header = FALSE)
palette <- palette_df$V2
names(palette) <- palette_df$V1

add_data <- function(g, mods, summary_function) {
  gotten <- g %>%
    activate(edges) %>%
    data.frame()
  gotten <- grep('M[0-9]{5}', colnames(gotten), value = TRUE)
  
  mods <- mods[names(mods) %in% gotten]
  
  da <- g %>%
    activate(edges) %>%
    select(names(mods)) %>%
    data.frame()
  
  for (m in names(da[, -c(1,2)])) {
    da[da[, m], m] <- mods[m]
  }
  
  if (is.null(dim(da[, -c(1:2)]))) {
    leca <- da[, -c(1:2)]
  } else {
    leca <- apply(da[, -c(1:2)], 1, summary_function)
  }

  g <- g %>%
    activate(edges) %>%
    mutate(LECA = leca)
  
  da <- g %>%
    activate(vertices) %>%
    select(names(mods)) %>%
    data.frame()
  
  for (m in names(da[, -c(1,2)])) {
    da[da[, m], m] <- mods[m]
  }
  
  if (is.null(dim(da[, -c(1:2)]))) {
    leca <- da[, -c(1:2)]
  } else {
    leca <- apply(da[, -c(1:2)], 1, summary_function)
  }
  
  g <- g %>%
    activate(vertices) %>%
    mutate(LECA = leca)
  
  return(g)
}

get_status_vec <- function(met, dat) {
  if (met == 'complete') {
    status_vec <- c('complete')
  } else if (met == 'incomplete') {
    status_vec <- unique(dat$status)
  } else if (met == 'missing') {
    status_vec <- unique(dat$status)[-grep('incomplete', unique(dat$status))]
  } else if (met == '1 block missing') {
    status_vec <- c('1 block missing', 'complete')
  } else {
    status_vec <- c('1 block missing', 'complete')
  }
  return(status_vec)
}

load('../outputs/metabolism/full_graph.RData')
odir <- '../outputs/metabolism_vircleaned/'
databases <- c('TOLDBA', 'TOLDBB', 'TOLDBC')
criteria <- c(3, 5)

plotsdir <- paste0(odir, '/plots')
dir.create(plotsdir, showWarnings = FALSE, recursive = TRUE)

dat <- c()
for (db in databases) {
  for (criterion in criteria) {
    prefix <- paste(odir, db, collapse = '/', sep = '')
    ddat <- read.csv(paste0(odir, '/', db, '_', criterion, 'sg_modules.txt'), sep = '\t')
    dat <- rbind(dat, data.frame(ddat, database = db, criterion = criterion))
  }
}

donors <- read.csv('../outputs/mLECA_origins_vircleaned/selected_groups.csv', sep = '\t')[, 1]
orig <- c()
for (db in databases) {
  for (criterion in criteria) {
    fnm <- paste0('../outputs/LECA_proteomes_vircleaned/', db, '_', criterion, 'sg_proteome.tsv')
    orig <- rbind(orig, data.frame(database = db, criterion = criterion, read.csv(fnm, sep = '\t')))
  }
}

modules_origin <- c()
for (i in 1:dim(dat)[1]) {
  prots <- str_split(dat[i, "gene_caller_ids_in_module"], ',', simplify = TRUE)[1, ]
  module <- dat[i, ] %>%
    select(module_category, module, module_name, enzyme_hits_in_module) %>%
    separate_longer_delim(enzyme_hits_in_module, ',')
  module$proteins <- unlist(lapply(module$enzyme_hits_in_module, function (x) paste(prots[[x]], sep = '', collapse = ';')))

  module <- module %>%
    separate_longer_delim(proteins, ';') %>%
    left_join(orig %>% select(family, origin, donor_domain, donor, criterion, database),
              by = c('proteins' = 'family', 'criterion' = 'criterion', 'database' = 'database'))
  modules_origin <- rbind(modules_origin, module)
}

# Plotting presence absence
for (db in databases) {
  for (crit in criteria) {
    fdat <- dat %>%
      filter(database == db,
             criterion == crit)

    mods <- rep(TRUE, length(unique(fdat$module)))
    names(mods) <- unique(fdat$module)
    
    a <- add_data(g, mods, summary_function = function(x) {sum(x) > 0})
    
    p <- ggraph(a, x = x, y = y, layout = 'manual') +
      geom_edge_link0(width = 0.5, color = 'grey95')+
      geom_edge_link(color = "steelblue", aes(filter = LECA > 0)) +
      geom_node_point(size = 1, color = "steelblue", aes(filter = LECA > 0)) +
      theme_void() +
      theme(plot.margin = unit(c(2, 2, 2, 2), units = 'mm')) +
      labs(title = paste(c(combs[, i], crit, stat), collapse = ' '))
    
    pdf(paste0(plotsdir, '/', bnm, '.pdf'), width = 11, height = 7)
    print(p)
    dev.off()
  }
}

combs <- combn(databases, 2)
for (crit in criteria) {
  for (i in 1:dim(combs)[2]) {
    fdat <- dat %>%
      filter(database %in% combs[, i], criterion == crit) %>%
      count(module)
    bnm <- paste(paste(combs[, i], collapse = '_'), crit, stat, sep = '_')
    
    mods <- rep(TRUE, length(unique(fdat$module)))
    names(mods) <- unique(fdat$module)
    
    a <- add_data(g, mods, summary_function = function(x) {sum(x) > 0})

    p <- ggraph(a, x = x, y = y, layout = 'manual') +
      geom_edge_link0(width = 0.5, color = 'grey95')+
      geom_edge_link(color = "steelblue", aes(filter = LECA > 0)) +
      geom_node_point(size = 1, color = "steelblue", aes(filter = LECA > 0)) +
      theme_void() +
      theme(plot.margin = unit(c(2, 2, 2, 2), units = 'mm')) +
      labs(title = paste(c(combs[, i], crit, stat), collapse = ' '))
    
    pdf(paste0(plotsdir, '/', bnm, '.pdf'), width = 11, height = 7)
    print(p)
    dev.off()
  }
}

combs <- combn(databases, 2)
for (crit in criteria) {
  fdat <- dat %>%
    filter(criterion == crit, status %in% get_status_vec(stat, dat)) %>%
    count(module)
  bnm <- paste(paste('ABC', collapse = '_'), crit, stat, sep = '_')

  mods <- rep(TRUE, length(unique(fdat$module)))
  names(mods) <- unique(fdat$module)
  
  a <- add_data(g, mods, summary_function = function(x) {sum(x) > 0})
  
  p <- ggraph(a, x = x, y = y, layout = 'manual') +
    geom_edge_link0(width = 0.5, color = 'grey95')+
    geom_edge_link(color = "steelblue", aes(filter = LECA > 0)) +
    geom_node_point(size = 1, color = "steelblue", aes(filter = LECA > 0)) +
    theme_void() +
    theme(plot.margin = unit(c(2, 2, 2, 2), units = 'mm')) +
    labs(title = paste(c(combs[, i], crit, stat), collapse = ' '))
  
  pdf(paste0(plotsdir, '/', bnm, '.pdf'), width = 11, height = 7)
  print(p)
  dev.off()
}

# Plotting the pervassiveness across datasets
modules_summ <- dat %>%
  group_by(database, criterion, module) %>%
  summarise() %>%
  ungroup() %>%
  count(module) %>%
  mutate(pervasiveness = n / 6)

pervassiveness <- modules_summ[, 3][[1]]
names(pervassiveness) <- modules_summ[, 1][[1]]

mods <- pervassiveness

summary_function <- function(x) {
  y = max(x[x != FALSE])
  if (is.infinite(y)) {
    return(FALSE)
  } else {
    return(y)
  }
}

a <- add_data(g, pervassiveness, summary_function)

p <- ggraph(a, x = x, y = y, layout = 'manual') +
  # geom_edge_link0(width = 0.3, aes(color = I(fgcolor),
  #                                  filter = type == 'line' & fgcolor != 'none'),
  #                 alpha = 0.5)+
  geom_edge_link0(width = 0.5, color = 'grey95')+
  geom_edge_link(linewidth = 1, aes(colour = I(fgcolor), filter = LECA > 0, alpha = LECA)) +
  geom_node_point(size = 1, aes(colour = I(fgcolor), filter = LECA > 0, alpha = LECA)) +
  # geom_edge_link0(width = 0.5, color = 'grey95')+
  # geom_edge_link(color = "steelblue", aes(filter = LECA > 0, alpha = LECA)) +
  # geom_node_point(size = 1, color = "steelblue", aes(filter = LECA > 0, alpha = LECA)) +
  theme_void() +
  theme(plot.margin = unit(c(2, 2, 2, 2), units = 'mm')) +
  labs(title = 'Total metabolism, intensity = pervasiveness across datasets, up to 2 blocks missing')
p

# Plotting the metabolic path according to origin
origin_summary <- modules_origin %>%
  group_by(module, status, database, criterion) %>%
  mutate(n_total = n(), orig_domain = ifelse(is.na(donor_domain), origin, donor_domain)) %>%
  group_by(type, group, module_name, module, status, database, criterion, orig_domain) %>%
  summarise(n = n(), prop = n / unique(n_total), module_origin = ifelse(prop >= 0.6, orig_domain, 'Unknown')) %>%
  slice_max(order_by = prop, n = 1, with_ties = FALSE) %>%
  arrange(-prop)

toni_origins <- modules_origin %>%
  mutate(dataset = paste0(database, '_', criterion)) %>%
  group_by(type, group, module, module_name, status, present_blocks, total_blocks,
           presence_prop, KOs, dataset) %>%
  mutate(donor = ifelse(origin == 'Innovation', 'Innovation', donor)) %>%
  summarise(donor_domain = paste(unique(donor_domain), sep = '', collapse = ';'),
            donor = paste(unique(donor), sep = '', collapse = ';')) %>%
  pivot_wider(names_from = dataset, values_from = c(donor_domain, donor))

toni_origins <- as.data.frame(toni_origins)

# write.table(toni_origins, file = '../outputs/metabolism/metabolism_summary.tsv',
#             sep = '\t', quote = FALSE, row.names = FALSE)

pdf(paste0(plotsdir, '/total_metabolism_by_donor.pdf'), width = 10, height = 6)
for (db in databases) {
  for (crit in criteria) {
    origins <- origin_summary %>%
      filter(status %in% get_status_vec('missing', dat),
             database == db, criterion == crit)
    mods <- origins[, "module_origin"][[1]]
    names(mods) <- origins[, "module"][[1]]
    
    summary_function <- function(x) {
      y <- unique(na.omit(x))
      if (length(y) > 1) {
        y <- y[y != FALSE]
        z <- paste(sort(y), sep = '', collapse = ';')
        if (z == 'Bacteria;Proteobacteria') {
          return('Bacteria')
        } else if (str_detect(z, ';')) {
          return('Unknown')
        } else {
          return(z)
        }
      } else {
        return(NA)
      }
    }
    
    a <- add_data(g, mods, summary_function)
    
    p <- ggraph(a, x = x, y = y, layout = 'manual') +
      geom_edge_link0(width = 0.5, color = 'grey95')+
      geom_edge_link(aes(filter = !is.na(LECA), colour = LECA)) +
      geom_node_point(size = 1, aes(filter = !is.na(LECA), colour = LECA)) +
      scale_colour_manual(values = palette) +
      scale_edge_colour_manual(values = palette) +
      theme_void() +
      theme(plot.margin = unit(c(2, 2, 2, 2), units = 'mm')) +
      labs(title = paste0('Total metabolism, colour = consensus donor, ', db, ', ', crit, ' supergroups up to 2 blocks missing'))
    print(p)
  }
}
dev.off()
