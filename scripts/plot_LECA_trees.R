# Plotting LECA families identified in modules
# Moisès Bernabeu
# Barcelona, January 2024

# Loading packages ----
library(tidyverse)
library(treeio)
library(ggplot2)
# library(ggpubr)
library(ape)
library(ggtree)
library(see)

theme_set(theme_bw())

source('../../../new_software/general/phygeno/ggcollapse/ggcollapse.R')

get_sp <- function(name) {
  if (str_detect(name, 'Prok\\|') & !str_detect(name, 'rvdb')) {
    return(str_split(gsub('Prok\\|', '', name), '_', simplify = TRUE)[, 2])
  } else if (str_detect(name, 'rvdb')) {
    return(str_split(gsub('Prok\\|', '', name), '_', simplify = TRUE)[, 1])
  } else {
    return(str_split(name, '_', simplify = TRUE)[, 1])
  }
}

annotate_tree <- function(tree, data, get_sp, sp_row=1) {
  require(stringr)
  require(tidyr)
  require(treeio)
  
  row.names(data) <- data[, sp_row]
  
  ttree <- as_tibble(tree)
  seqids <- ttree$label[1:Ntip(ttree)]
  taxids <- unlist(lapply(seqids, get_sp))
  
  tree_dat <- tibble(label = seqids, data[taxids, ])
  otree <- full_join(ttree, tree_dat, by = 'label')
  otree <- as.treedata(otree)
  
  return(otree)
}

trees <- read.csv('../outputs/large_ogs/round_1_LECA_trees.nwk', sep = '\t', header = FALSE)
dat <- read.csv('../data/LECA_project.lng', sep = '\t', header = FALSE)

dat <- dat %>%
  mutate(V2 = gsub('[a-z]{1,4}__', '', V2)) %>%
  separate_wider_delim(V2, ';', names_sep = '_') %>%
  distinct() %>%
  as.data.frame()

cols <- okabeito_colors(1:4)
names(cols) <- unique(dat$V2_4)

odir <- '../outputs/large_ogs/plots'
dir.create(odir)
pdf(paste0(odir, '/leca_trees.pdf'),
    width = 8,
    height = 8)
for (i in 1:dim(trees)[1]) {
  tmp <- tempfile(pattern = 'tmp', tmpdir = odir, fileext = '.txt')
  cat(trees[i, 4], file = tmp)
  tr <- read.nhx(tmp)
  file.remove(tmp)
  
  tr <- annotate_tree(tr, dat, get_sp = get_sp, sp_row = 1)
  
  p <- ggtree(tr, ladderize = TRUE) +
    geom_tippoint(aes(colour = V2_4)) +
    geom_nodepoint(aes(subset = name == 'LECA'), colour = 'white', size = 2.5, pch = 21, fill = 'black') +
    scale_colour_manual(values = cols) +
    labs(title = paste(trees[i, 1], trees[i, 2], sep = ' | '))
  
  # pdf(paste0(odir, '/', trees[i, 1], '_', trees[i, 2], '.pdf'),
  #     width = 8,
  #     height = 8)
  print(p)
}
dev.off()

