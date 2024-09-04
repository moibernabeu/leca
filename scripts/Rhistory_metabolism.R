geom_edge_link0(width = 0.5, color = 'grey95')+
geom_edge_link(color = "steelblue", aes(filter = LECA > 0)) +
geom_node_point(size = 1, color = "steelblue", aes(filter = LECA > 0)) +
theme_void() +
theme(plot.margin = unit(c(2, 2, 2, 2), units = 'mm')) +
labs(title = 'Metaproteome', collapse = ' '))
print(p)
pdf(paste0(plotsdir, '/presence.pdf'), width = 11, height = 7)
for (db in databases) {
for (crit in criteria) {
fdat <- dat %>%
filter(database == db,
criterion >= crit)
mods <- rep(TRUE, length(unique(fdat$module)))
names(mods) <- unique(fdat$module)
a <- add_data(g, mods, summary_function = function(x) {sum(x) > 0})
p <- ggraph(a, x = x, y = y, layout = 'manual') +
geom_edge_link0(width = 0.5, color = 'grey95')+
geom_edge_link(color = "steelblue", aes(filter = LECA > 0)) +
geom_node_point(size = 1, color = "steelblue", aes(filter = LECA > 0)) +
theme_void() +
theme(plot.margin = unit(c(2, 2, 2, 2), units = 'mm')) +
labs(title = paste(c(db, crit, 'supergroups'), collapse = ' '))
print(p)
}
}
pdf(paste0(plotsdir, '/presence.pdf'), width = 11, height = 7)
for (db in databases) {
for (crit in criteria) {
fdat <- dat %>%
filter(database == db,
criterion >= crit)
mods <- rep(TRUE, length(unique(fdat$module)))
names(mods) <- unique(fdat$module)
a <- add_data(g, mods, summary_function = function(x) {sum(x) > 0})
p <- ggraph(a, x = x, y = y, layout = 'manual') +
geom_edge_link0(width = 0.5, color = 'grey95')+
geom_edge_link(color = "steelblue", aes(filter = LECA > 0)) +
geom_node_point(size = 1, color = "steelblue", aes(filter = LECA > 0)) +
theme_void() +
theme(plot.margin = unit(c(2, 2, 2, 2), units = 'mm')) +
labs(title = paste(c(db, crit, 'supergroups'), collapse = ' '))
print(p)
}
}
# Plotting metaproteome
fdat <- read.csv(paste0(odir, '/metaproteome_modules.txt'), sep = '\t') %>%
filter(stepwise_module_is_complete == 'True')
mods <- rep(TRUE, length(unique(fdat$module)))
names(mods) <- unique(fdat$module)
a <- add_data(g, mods, summary_function = function(x) {sum(x) > 0})
p <- ggraph(a, x = x, y = y, layout = 'manual') +
geom_edge_link0(width = 0.5, color = 'grey95')+
geom_edge_link(color = "steelblue", aes(filter = LECA > 0)) +
geom_node_point(size = 1, color = "steelblue", aes(filter = LECA > 0)) +
theme_void() +
theme(plot.margin = unit(c(2, 2, 2, 2), units = 'mm')) +
labs(title = 'Metaproteome')
print(p)
dev.off()
databases
combs <- combn(databases, 2)
combs
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
combs[, i]
1:dim(combs)[2]
fdat <- dat %>%
filter(database %in% combs[, i], criterion == crit) %>%
count(module)
fdat
bnm <- paste(paste(combs[, i], collapse = '_'), crit, stat, sep = '_')
combs[, i]
crit
paste(combs[, i], collapse = '_')
paste(combs[, i], collapse = ' and ')
bnm <- paste(paste(combs[, i], collapse = ' and '), crit, sep = ' supergroups')
bnm
bnm <- paste(paste(combs[, i], collapse = ' and '), crit, 'supergroups', sep = ' ')
bnm
mods <- rep(TRUE, length(unique(fdat$module)))
mods
names(mods) <- unique(fdat$module)
mods
a <- add_data(g, mods, summary_function = function(x) {sum(x) > 0})
p <- ggraph(a, x = x, y = y, layout = 'manual') +
geom_edge_link0(width = 0.5, color = 'grey95')+
geom_edge_link(color = "steelblue", aes(filter = LECA > 0)) +
geom_node_point(size = 1, color = "steelblue", aes(filter = LECA > 0)) +
theme_void() +
theme(plot.margin = unit(c(2, 2, 2, 2), units = 'mm')) +
labs(title = paste(c(combs[, i], crit, stat), collapse = ' '))
print(p)
crit
pdf(paste0(plotsdir, '/TOLDB_union.pdf'), width = 11, height = 7)
pdf(paste0(plotsdir, '/TOLDB_union.pdf'), width = 11, height = 7)
for (crit in criteria) {
for (i in 1:dim(combs)[2]) {
fdat <- dat %>%
filter(database %in% combs[, i], criterion == crit) %>%
count(module)
bnm <- paste(paste(combs[, i], collapse = ' and '), crit, 'supergroups', sep = ' ')
mods <- rep(TRUE, length(unique(fdat$module)))
names(mods) <- unique(fdat$module)
a <- add_data(g, mods, summary_function = function(x) {sum(x) > 0})
p <- ggraph(a, x = x, y = y, layout = 'manual') +
geom_edge_link0(width = 0.5, color = 'grey95')+
geom_edge_link(color = "steelblue", aes(filter = LECA > 0)) +
geom_node_point(size = 1, color = "steelblue", aes(filter = LECA > 0)) +
theme_void() +
theme(plot.margin = unit(c(2, 2, 2, 2), units = 'mm')) +
labs(title = paste(c(combs[, i], crit), collapse = ' '))
print(p)
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
print(p)
}
p <- ggraph(a, x = x, y = y, layout = 'manual') +
geom_edge_link0(width = 0.5, color = 'grey95')+
geom_edge_link(color = "steelblue", aes(filter = LECA > 0)) +
geom_node_point(size = 1, color = "steelblue", aes(filter = LECA > 0)) +
theme_void() +
theme(plot.margin = unit(c(2, 2, 2, 2), units = 'mm')) +
labs(title = paste('ABC', crit, 'supergroups', collapse = ' '))
pdf(paste0(plotsdir, '/TOLDB_union.pdf'), width = 11, height = 7)
combs <- combn(databases, 2)
for (crit in criteria) {
for (i in 1:dim(combs)[2]) {
fdat <- dat %>%
filter(database %in% combs[, i], criterion == crit) %>%
count(module)
bnm <- paste(paste(combs[, i], collapse = ' and '), crit, 'supergroups', sep = ' ')
mods <- rep(TRUE, length(unique(fdat$module)))
names(mods) <- unique(fdat$module)
a <- add_data(g, mods, summary_function = function(x) {sum(x) > 0})
p <- ggraph(a, x = x, y = y, layout = 'manual') +
geom_edge_link0(width = 0.5, color = 'grey95')+
geom_edge_link(color = "steelblue", aes(filter = LECA > 0)) +
geom_node_point(size = 1, color = "steelblue", aes(filter = LECA > 0)) +
theme_void() +
theme(plot.margin = unit(c(2, 2, 2, 2), units = 'mm')) +
labs(title = paste(c(combs[, i], crit), collapse = ' '))
print(p)
}
}
combs <- combn(databases, 2)
for (crit in criteria) {
fdat <- dat %>%
filter(criterion == crit, status %in% get_status_vec(stat, dat)) %>%
count(module)
mods <- rep(TRUE, length(unique(fdat$module)))
names(mods) <- unique(fdat$module)
a <- add_data(g, mods, summary_function = function(x) {sum(x) > 0})
p <- ggraph(a, x = x, y = y, layout = 'manual') +
geom_edge_link0(width = 0.5, color = 'grey95')+
geom_edge_link(color = "steelblue", aes(filter = LECA > 0)) +
geom_node_point(size = 1, color = "steelblue", aes(filter = LECA > 0)) +
theme_void() +
theme(plot.margin = unit(c(2, 2, 2, 2), units = 'mm')) +
labs(title = paste('ABC', crit, 'supergroups', collapse = ' '))
print(p)
}
pdf(paste0(plotsdir, '/TOLDB_union.pdf'), width = 11, height = 7)
combs <- combn(databases, 2)
for (crit in criteria) {
for (i in 1:dim(combs)[2]) {
fdat <- dat %>%
filter(database %in% combs[, i], criterion == crit) %>%
count(module)
bnm <- paste(paste(combs[, i], collapse = ' and '), crit, 'supergroups', sep = ' ')
mods <- rep(TRUE, length(unique(fdat$module)))
names(mods) <- unique(fdat$module)
a <- add_data(g, mods, summary_function = function(x) {sum(x) > 0})
p <- ggraph(a, x = x, y = y, layout = 'manual') +
geom_edge_link0(width = 0.5, color = 'grey95')+
geom_edge_link(color = "steelblue", aes(filter = LECA > 0)) +
geom_node_point(size = 1, color = "steelblue", aes(filter = LECA > 0)) +
theme_void() +
theme(plot.margin = unit(c(2, 2, 2, 2), units = 'mm')) +
labs(title = paste(c(combs[, i], crit), collapse = ' '))
print(p)
}
}
combs <- combn(databases, 2)
for (crit in criteria) {
fdat <- dat %>%
filter(criterion == crit) %>%
count(module)
mods <- rep(TRUE, length(unique(fdat$module)))
names(mods) <- unique(fdat$module)
a <- add_data(g, mods, summary_function = function(x) {sum(x) > 0})
p <- ggraph(a, x = x, y = y, layout = 'manual') +
geom_edge_link0(width = 0.5, color = 'grey95')+
geom_edge_link(color = "steelblue", aes(filter = LECA > 0)) +
geom_node_point(size = 1, color = "steelblue", aes(filter = LECA > 0)) +
theme_void() +
theme(plot.margin = unit(c(2, 2, 2, 2), units = 'mm')) +
labs(title = paste('ABC', crit, 'supergroups', collapse = ' '))
print(p)
}
dev.off()
pdf(paste0(plotsdir, '/TOLDBs_union.pdf'), width = 11, height = 7)
combs <- combn(databases, 2)
for (crit in criteria) {
for (i in 1:dim(combs)[2]) {
fdat <- dat %>%
filter(database %in% combs[, i], criterion == crit) %>%
count(module)
bnm <- paste(paste(combs[, i], collapse = ' and '), crit, 'supergroups', sep = ' ')
mods <- rep(TRUE, length(unique(fdat$module)))
names(mods) <- unique(fdat$module)
a <- add_data(g, mods, summary_function = function(x) {sum(x) > 0})
p <- ggraph(a, x = x, y = y, layout = 'manual') +
geom_edge_link0(width = 0.5, color = 'grey95')+
geom_edge_link(color = "steelblue", aes(filter = LECA > 0)) +
geom_node_point(size = 1, color = "steelblue", aes(filter = LECA > 0)) +
theme_void() +
theme(plot.margin = unit(c(2, 2, 2, 2), units = 'mm')) +
labs(title = bnm)
print(p)
}
}
combs <- combn(databases, 2)
for (crit in criteria) {
fdat <- dat %>%
filter(criterion == crit) %>%
count(module)
mods <- rep(TRUE, length(unique(fdat$module)))
names(mods) <- unique(fdat$module)
a <- add_data(g, mods, summary_function = function(x) {sum(x) > 0})
p <- ggraph(a, x = x, y = y, layout = 'manual') +
geom_edge_link0(width = 0.5, color = 'grey95')+
geom_edge_link(color = "steelblue", aes(filter = LECA > 0)) +
geom_node_point(size = 1, color = "steelblue", aes(filter = LECA > 0)) +
theme_void() +
theme(plot.margin = unit(c(2, 2, 2, 2), units = 'mm')) +
labs(title = paste('ABC', crit, 'supergroups', collapse = ' '))
print(p)
}
dev.off()
# Plotting the pervassiveness across datasets
modules_summ <- dat %>%
group_by(database, criterion, module) %>%
summarise() %>%
ungroup() %>%
count(module) %>%
mutate(pervasiveness = n / 6)
modules_summ
pervassiveness <- modules_summ[, 3][[1]]
pervassiveness
names(pervassiveness) <- modules_summ[, 1][[1]]
mods <- pervassiveness
mods
summary_function <- function(x) {
y = max(x[x != FALSE])
if (is.infinite(y)) {
return(FALSE)
} else {
return(y)
}
}
a <- add_data(g, pervassiveness, summary_function)
a
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
labs(title = 'Total metabolism, intensity = pervasiveness across datasets')
p
# Plotting the metabolic path according to origin
origin_summary <- modules_origin %>%
group_by(module, status, database, criterion) %>%
mutate(n_total = n(), orig_domain = ifelse(is.na(donor_domain), origin, donor_domain)) %>%
group_by(type, group, module_name, module, status, database, criterion, orig_domain) %>%
summarise(n = n(), prop = n / unique(n_total), module_origin = ifelse(prop >= 0.6, orig_domain, 'Unknown')) %>%
slice_max(order_by = prop, n = 1, with_ties = FALSE) %>%
arrange(-prop)
# Plotting the metabolic path according to origin
origin_summary <- modules_origin %>%
group_by(module, database, criterion) %>%
mutate(n_total = n(), orig_domain = ifelse(is.na(donor_domain), origin, donor_domain)) %>%
group_by(type, group, module_name, module, status, database, criterion, orig_domain) %>%
summarise(n = n(), prop = n / unique(n_total), module_origin = ifelse(prop >= 0.6, orig_domain, 'Unknown')) %>%
slice_max(order_by = prop, n = 1, with_ties = FALSE) %>%
arrange(-prop)
modules_origin
modules_origin %>%
group_by(module_category, module, database, criterion)
modules_origin %>%
group_by(module_category, module, database, criterion) %>%
mutate(n_total = n(), orig_domain = ifelse(is.na(donor_domain), origin, donor_domain))
modules_origin %>%
group_by(module_category, module, database, criterion) %>%
mutate(n_total = n(), orig_domain = ifelse(is.na(donor_domain), origin, donor_domain)) %>%
group_by(type, group, module_name, module, status, database, criterion, orig_domain)
modules_origin %>%
group_by(module_category, module, module_name, database, criterion) %>%
mutate(n_total = n(), orig_domain = ifelse(is.na(donor_domain), origin, donor_domain)) %>%
group_by(module_category, module, module_name, database, criterion, orig_domain)
modules_origin %>%
group_by(module_category, module, module_name, database, criterion) %>%
mutate(n_total = n(), orig_domain = ifelse(is.na(donor_domain), origin, donor_domain)) %>%
group_by(module_category, module, module_name, database, criterion, orig_domain) %>%
summarise(n = n(), prop = n / unique(n_total), module_origin = ifelse(prop >= 0.6, orig_domain, 'Unknown'))
modules_origin %>%
group_by(module_category, module, module_name, database, criterion) %>%
mutate(n_total = n(), orig_domain = ifelse(is.na(donor_domain), origin, donor_domain)) %>%
group_by(module_category, module, module_name, database, criterion, orig_domain) %>%
summarise(n = n(), prop = n / unique(n_total), module_origin = ifelse(prop >= 0.6, orig_domain, 'Unknown')) %>%
slice_max(order_by = prop, n = 1, with_ties = FALSE)
modules_origin %>%
group_by(module_category, module, module_name, database, criterion) %>%
mutate(n_total = n(), orig_domain = ifelse(is.na(donor_domain), origin, donor_domain)) %>%
group_by(module_category, module, module_name, database, criterion, orig_domain) %>%
summarise(n = n(), prop = n / unique(n_total), module_origin = ifelse(prop >= 0.6, orig_domain, 'Unknown')) %>%
slice_max(order_by = prop, n = 1, with_ties = FALSE) %>%
arrange(-prop)
# Plotting the metabolic path according to origin
origin_summary <- modules_origin %>%
group_by(module_category, module, module_name, database, criterion) %>%
mutate(n_total = n(), orig_domain = ifelse(is.na(donor_domain), origin, donor_domain)) %>%
group_by(module_category, module, module_name, database, criterion, orig_domain) %>%
summarise(n = n(), prop = n / unique(n_total), module_origin = ifelse(prop >= 0.6, orig_domain, 'Unknown')) %>%
slice_max(order_by = prop, n = 1, with_ties = FALSE) %>%
arrange(-prop)
View(origin_summary)
modules_origin
modules_origin
modules_origin %>%
mutate(dataset = paste0(database, '_', criterion)) %>%
group_by(module_category, module, module_name, path_KO)
modules_origin %>%
mutate(dataset = paste0(database, '_', criterion)) %>%
group_by(module_category, module, module_name, path_KO) %>%
mutate(donor = ifelse(origin == 'Innovation', 'Innovation', donor))
modules_origin %>%
mutate(dataset = paste0(database, '_', criterion)) %>%
group_by(module_category, module, module_name, path_KO) %>%
mutate(donor = ifelse(origin == 'Innovation', 'Innovation', donor)) %>%
summarise(donor_domain = paste(unique(donor_domain), sep = '', collapse = ';'),
donor = paste(unique(donor), sep = '', collapse = ';'))
toni_origins <- modules_origin %>%
mutate(dataset = paste0(database, '_', criterion)) %>%
group_by(module_category, module, module_name, path_KO) %>%
mutate(donor = ifelse(origin == 'Innovation', 'Innovation', donor)) %>%
summarise(donor_domain = paste(unique(donor_domain), sep = '', collapse = ';'),
donor = paste(unique(donor), sep = '', collapse = ';')) %>%
pivot_wider(names_from = dataset, values_from = c(donor_domain, donor))
dataset
toni_origins <- modules_origin %>%
mutate(dataset = paste0(database, '_', criterion)) %>%
group_by(module_category, module, module_name, path_KO) %>%
mutate(donor = ifelse(origin == 'Innovation', 'Innovation', donor)) %>%
summarise(donor_domain = paste(unique(donor_domain), sep = '', collapse = ';'),
donor = paste(unique(donor), sep = '', collapse = ';')) %>%
pivot_wider(names_from = dataset, values_from = c(donor_domain, donor))
modules_origin %>%
mutate(dataset = paste0(database, '_', criterion)) %>%
group_by(module_category, module, module_name, path_KO)
toni_origins <- modules_origin %>%
mutate(dataset = paste0(database, '_', criterion)) %>%
group_by(dataset, module_category, module, module_name, path_KO) %>%
mutate(donor = ifelse(origin == 'Innovation', 'Innovation', donor)) %>%
summarise(donor_domain = paste(unique(donor_domain), sep = '', collapse = ';'),
donor = paste(unique(donor), sep = '', collapse = ';')) %>%
pivot_wider(names_from = dataset, values_from = c(donor_domain, donor))
toni_origins <- as.data.frame(toni_origins)
toni_origins
write.table(toni_origins, file = '../outputs/metabolism/metabolism_summary.tsv',
sep = '\t', quote = FALSE, row.names = FALSE)
pdf(paste0(plotsdir, '/total_metabolism_by_donor.pdf'), width = 10, height = 6)
for (db in databases) {
for (crit in criteria) {
origins <- origin_summary %>%
filter(database == db, criterion == crit)
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
labs(title = paste0('Total metabolism, colour = consensus donor, ', db, ', ', crit, ' supergroups'))
print(p)
}
}
dev.off()
write.table(toni_origins, file = '../outputs/metabolism_vircleaned/metabolism_summary.tsv',
sep = '\t', quote = FALSE, row.names = FALSE)
# toni_origins <-
modules_origin %>%
mutate(dataset = paste0(database, '_', criterion)) %>%
group_by(dataset, module_category, module, module_name, path_KO) %>%
mutate(donor = ifelse(origin == 'Innovation', 'Innovation', donor)) %>%
summarise(donor_domain = paste(unique(donor_domain), sep = '', collapse = ';'),
donor = paste(unique(donor), sep = '', collapse = ';')) %>%
pivot_wider(names_from = dataset, values_from = c(donor_domain, donor))
# toni_origins <-
modules_origin %>%
mutate(dataset = paste0(database, '_', criterion)) %>%
group_by(dataset, module_category, module, module_name) %>%
mutate(donor = ifelse(origin == 'Innovation', 'Innovation', donor)) %>%
summarise(donor_domain = paste(unique(donor_domain), sep = '', collapse = ';'),
donor = paste(unique(donor), sep = '', collapse = ';')) %>%
pivot_wider(names_from = dataset, values_from = c(donor_domain, donor))
modules_origin %>%
mutate(dataset = paste0(database, '_', criterion)) %>%
group_by(dataset, module_category, module, module_name) %>%
mutate(donor = ifelse(origin == 'Innovation', 'Innovation', donor)) %>%
summarise(donor_domain = paste(unique(donor_domain), sep = '', collapse = ';'),
donor = paste(unique(donor), sep = '', collapse = ';')) %>%
pivot_wider(names_from = dataset, values_from = c(donor_domain, donor))
toni_origins <- as.data.frame(toni_origins)
write.table(toni_origins, file = '../outputs/metabolism_vircleaned/metabolism_summary_KOs.tsv',
sep = '\t', quote = FALSE, row.names = FALSE)
toni_origins <- modules_origin %>%
mutate(dataset = paste0(database, '_', criterion)) %>%
group_by(dataset, module_category, module, module_name, path_KO) %>%
mutate(donor = ifelse(origin == 'Innovation', 'Innovation', donor)) %>%
summarise(donor_domain = paste(unique(donor_domain), sep = '', collapse = ';'),
donor = paste(unique(donor), sep = '', collapse = ';')) %>%
pivot_wider(names_from = dataset, values_from = c(donor_domain, donor))
toni_origins <- as.data.frame(toni_origins)
write.table(toni_origins, file = '../outputs/metabolism_vircleaned/metabolism_summary_KOs.tsv',
sep = '\t', quote = FALSE, row.names = FALSE)
toni_origins <- modules_origin %>%
mutate(dataset = paste0(database, '_', criterion)) %>%
group_by(dataset, module_category, module, module_name) %>%
mutate(donor = ifelse(origin == 'Innovation', 'Innovation', donor)) %>%
summarise(donor_domain = paste(unique(donor_domain), sep = '', collapse = ';'),
donor = paste(unique(donor), sep = '', collapse = ';')) %>%
pivot_wider(names_from = dataset, values_from = c(donor_domain, donor))
toni_origins <- as.data.frame(toni_origins)
write.table(toni_origins, file = '../outputs/metabolism_vircleaned/metabolism_summary.tsv',
sep = '\t', quote = FALSE, row.names = FALSE)
