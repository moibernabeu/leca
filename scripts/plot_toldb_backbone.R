library(treeio)
library(ggtree)
library(see)
library(tidyverse)

pal_df <- read.csv('../data/palette.txt', sep = '\t', header = FALSE)
pal <- pal_df$V2
names(pal) <- pal_df$V1

tree <- read.newick('../data/LECA_backbone.nwk')
tdat <- read.csv('../data/toldb_sp_distr.tsv', sep = '\t')

tdat <- tdat %>%
  mutate(label = paste0(Division, ' (', ToLDB, ')')) %>%
  select(Division, Stem, Supergroup, tip_label = label)

tdat <- rbind(tdat, c('BaSk', 'Stem 2', 'Metamonada', 'BaSk (0)'), c('Tsukubamonadida', 'Stem 1', 'Discoba', 'Tsukubamonadida (0)')) %>%
  as.data.frame()
row.names(tdat) <- tdat$Division

ggtree(tree) %<+% tdat +
  geom_tiplab(aes(colour = Stem, label = tip_label)) +
  xlim(0, 9) +
  geom_nodelab(aes(subset = label != '', x = branch), geom = 'label') +
  geom_nodelab(aes(label = node), geom = 'label') +
  scale_colour_manual(values = c('Stem 2' = 'darkorange3', 'Stem 1' = 'steelblue'))

ggtree(tree, layout = 'fan', open.angle = 180, size = 0.75) %<+% tdat +
  geom_tiplab(aes(colour = Stem, label = node), offset = 0.15, show.legend = FALSE) +
  xlim_tree(c(0, 12)) +
  geom_hilight(node = 38, to.bottom = TRUE, gradient.direction = 'tr', extend = 0.1, type = 'gradient') +
  geom_hilight(node = c(53, 34), to.bottom = TRUE, gradient.direction = 'rt', extend = 0.1, fill = 'darkorange3') +
  scale_colour_manual(values = c('Stem 2' = 'darkorange3', 'Stem 1' = 'steelblue'))

pdf('../paper/eToL_consensus.pdf', width = 10, height = 6)
ggtree(tree, layout = 'fan', open.angle = 180, size = 0.75) %<+% tdat +
  geom_tiplab(offset = 0.15, size = 2.7, show.legend = FALSE) +
  xlim_tree(c(0, 12)) +
  geom_hilight(node = 38, to.bottom = TRUE, gradient.direction = 'tr', extend = 0.1, type = 'gradient', fill = pal['Stem 1']) +
  geom_hilight(node = c(53, 34), to.bottom = TRUE, gradient.direction = 'rt', extend = 0.1, fill = pal['Stem 2']) +
  # geom_nodelab(aes(label = node), geom = 'label') +
  scale_color_manual(values = pal) +
  geom_hilight(node = c(48, 45, 39, 42, 55, 58, 34, 25, 24),
               fill = c(pal['Archaeplastida_plus'],
                        pal['TRASH'],
                        pal['Discoba'],
                        pal['PHM'],
                        pal['Amorphea'],
                        pal['CRuMs'],
                        pal['Metamonada'],
                        pal['Malawimonadida'],
                        pal['Ancyromonadida']),
               to.bottom = TRUE, extend  = 7.5, alpha = 0.7) +
  geom_cladelab(node = c(48, 45, 39, 42, 55, 58, 34, 25, 24),
                label = c('Archaeplastida_plus',
                          'TRASH',
                          'Discoba',
                          'PHM',
                          'Amorphea',
                          'CRuMs',
                          'Metamonada',
                          'M',
                          'A'),
                offset = 7.5, extend = 0.5, barsize = 2, barcolor = pal[c('Archaeplastida_plus',
                                                                         'TRASH',
                                                                         'Discoba',
                                                                         'PHM',
                                                                         'Amorphea',
                                                                         'CRuMs',
                                                                         'Metamonada',
                                                                         'Malawimonadida',
                                                                         'Ancyromonadida')],
                textcolor = pal[c('Archaeplastida_plus',
                                  'TRASH',
                                  'Discoba',
                                  'PHM',
                                  'Amorphea',
                                  'CRuMs',
                                  'Metamonada',
                                  'Malawimonadida',
                                  'Ancyromonadida')],
                offset.text = 1, angle = c(77, 45, 20, 0, -20, -40, -75, -50, -60), hjust = 0.5, textsize = 3) +
  geom_strip('Rhodophyta', 'Provora', label = 'Stem 1', offset = 9.6, angle = 40, offset.text = 1) +
  geom_strip('Opisthokonta', 'Anaeramoebidae', label = 'Stem 2', offset = 9.6, angle = -52, offset.text = 1)
dev.off()
