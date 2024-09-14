library(tidyverse)
library(ggtree)
library(treeio)

pal_df <- read.csv('../data/colours.tsv', sep = '\t', header = FALSE)
pal <- pal_df$V2
names(pal) <- pal_df$V1

taxonomy <- read.csv('../data/euk_taxonomy_toldb.tsv', sep = '\t')

tr <- read.tree('../outputs/sptree_TOLDBA/concatenate/iqtree_LGG4PMSF/TOLDBA.treefile')

ggtree(tr) %<+% taxonomy +
  geom_node_label(aes(label = node)) +
  geom_tiplab(aes(label = Supergroup, colour = Supergroup), size = 2.5)

ggtree(root(tr, node = 178)) %<+% taxonomy +
  # geom_nodelab(aes(label = node)) +
  geom_tiplab(aes(label = Supergroup, colour = Supergroup), size = 2.5) +
  # geom_nodelab(aes(label = label)) +
  labs(title = 'LG+G4+PMSF')

# tr <- read.tree('../outputs/sptree_TOLDBA/concatenate/iqtree/TOLDBA.bionj')
# ggtree(root(tr, node = 138)) %<+% taxonomy +
#   geom_tiplab(aes(label = Supergroup, colour = Supergroup)) +
#   geom_nodelab(aes(label = node))
