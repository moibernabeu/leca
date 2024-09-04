library(tidyverse)

taxonomy <- read.csv('../data/euk_taxonomy_toldb.tsv', sep = '\t')

taxonomy <- taxonomy %>%
  select(MNEMO, stem = Stem, sg = Supergroup, div = Division, d = d,
         p = p, c = c, o = o, f = f, g=g, s = Species)

out_lng <- c()
for (i in 1:dim(taxonomy)[1]) {
  lng <- paste(names(taxonomy)[-1], taxonomy[i, -1], sep = '__')
  lng <- paste0(lng, collapse = ';')
  out_lng <- rbind(out_lng, data.frame(MNEMO = taxonomy[i, 1],
                                       lng = lng))
}
out_lng

write.table(out_lng, file = '../data/toldb.lng', sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
