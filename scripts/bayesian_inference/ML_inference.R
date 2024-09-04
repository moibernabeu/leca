# Gamma parameters and mode inference
# Moisès Bernabeu
# Barcelona, October 2022

# Loading libraries ----
library(fitdistrplus)
library(ggplot2)
library(ggpubr)
library(latex2exp)

theme_set(theme_bw())

# Loading data ----
load('../../../02_dist_stats/data/dist_data.RData')
nodes <- read.csv('../../../02_dist_stats/data/node_labs.tsv', sep = '\t', row.names = 3)

# Plotting function ----
plot_diag <- function(ginfer, ninfer, title) {
  dat <- data.frame(y = ginfer['data']$data)

  x <- seq(range(dat)[1], range(dat)[2],  0.001)

  infdat <- data.frame(x = c(x, x),
                       y = c(dgamma(x, ginfer$estimate[1], ginfer$estimate[2]),
                             dnorm(x, ninfer$estimate[1], ninfer$estimate[2])),
                       distr = c(rep('Gamma', length(x)), rep('Normal', length(x))))

  dparams <- list('Gamma' = ginfer$estimate, 'Normal' = ninfer$estimate)

  a <- ggplot(dat) +
    geom_histogram(aes(x = y, y = ..density..), alpha = 0.4, colour = 'black') +
    geom_line(data = infdat, aes(x, y, colour = distr, lty = distr), size = 0.75) +
    scale_colour_manual(values = c('Gamma' = 'steelblue', 'Normal' = 'darkorange3')) +
    scale_linetype_manual(values = c('Gamma' = 1, 'Normal' = 4)) +
    ylab('Density') +
    xlab('Distance') +
    labs(colour = 'Inferred distribution', lty = 'Inferred distribution')

  b <- ggplot(dat, aes(sample = y)) +
    geom_abline(slope = 1) +
    geom_qq(distribution = qnorm, dparams = dparams$Normal, colour = 'darkorange3') +
    geom_qq(distribution = qgamma, dparams = dparams$Gamma, colour = 'steelblue') +
    xlab('Theoretical') +
    ylab('Observed')

  legtext <- paste0('BIC\nGamma: ', round(ginfer$bic, 3),
                    '\nNormal: ', round(ninfer$bic, 3))

  c <- ggplot() +
    stat_ecdf(data = dat, aes(y), geom = 'point') +
    geom_line(aes(x = x, y = pnorm(x, dparams$Normal[1], dparams$Normal[2])),
              colour = 'darkorange3', lty = 4, size = 0.75) +
    geom_line(aes(x = x, y = pgamma(x, dparams$Gamma[1], dparams$Gamma[2])),
              colour = 'steelblue', lty = 1, size = 0.75) +
    ylab('CDF') +
    xlab('Distance') +
    geom_label(aes(x = max(dat$y), y = 0, label = legtext), nudge_y = 0.1,
               vjust = 0, hjust = 1)

  p <- ggarrange(a, b, c, common.legend = TRUE, nrow = 1, align = 'hv')
  p <- annotate_figure(p, top = text_grob(title, color = 'black',
                                     face = 'bold', size = 14))
  print(p)
}

summary_ML <- function(MLobject) {
  out <- c(MLobject$estimate,
           Log_likelihood = MLobject$loglik,
           aic = MLobject$aic,
           bic = MLobject$bic)
  return(out)
}

# Inferring gamma distributions ----
MLgam <- by(fgnew$ndist, fgnew$node, fitdist, distr = 'gamma', method = 'mle')
MLnor <- by(fgnew$ndist, fgnew$node, fitdist, distr = 'norm', method = 'mle')

# pdf('../../outputs/ML/ev_diagnostics.pdf', width = 10, height = 3, onefile = TRUE)
for (node in names(MLgam)) {
  plot_diag(MLgam[[node]], MLnor[[node]], nodes[node, 'Group'])
}
# dev.off()

gamsum <- t(data.frame(lapply(MLgam, summary_ML)))
norsum <- t(data.frame(lapply(MLnor, summary_ML)))

gamsum <- data.frame(gamsum, dist = 'Gamma')
norsum <- data.frame(norsum, dist = 'Normal')

write.csv(file = '../../outputs/ML/ev_gamsum.csv', gamsum)
write.csv(file = '../../outputs/ML/ev_norsum.csv', norsum)

sum(gamsum$bic < norsum$bic) / length(norsum$bic)

# Species to species distances ----
load('../../../02_dist_stats/data/sp2sp_dat.RData')

datc <- datc[-which(datc$ndist == 0), ]
MLgam <- by(datc$ndist, datc$sp_to, fitdist, distr = 'gamma', method = 'mle')
MLnor <- by(datc$ndist, datc$sp_to, fitdist, distr = 'norm', method = 'mle')

# pdf('../../outputs/ML/sp2sp_diagnostics.pdf', width = 10, height = 3, onefile = TRUE)
for (node in names(MLgam)) {
  plot_diag(MLgam[[node]], MLnor[[node]], node)
}
# dev.off()

gamsum <- t(data.frame(lapply(MLgam, summary_ML)))
norsum <- t(data.frame(lapply(MLnor, summary_ML)))

gamsum <- data.frame(gamsum, dist = 'Gamma')
norsum <- data.frame(norsum, dist = 'Normal')

write.csv(file = '../../outputs/ML/sp2sp_gamsum.csv', gamsum)
write.csv(file = '../../outputs/ML/sp2sp_norsum.csv', norsum)

sum(gamsum$bic < norsum$bic) / length(norsum$bic)

lm1 <- lm(data = gamsum,  rate ~ I(1 / shape))
summary(lm1)

plot(gamsum$shape, gamsum$rate)
lines(1:500/100, predict(lm1, newdata = data.frame(shape = 1:500/100)))

library(ggrepel)

ggplot(gamsum, aes(shape, rate)) +
  geom_point() +
  geom_smooth(formula = y ~ I(1 / x^6), method = 'lm') +
  geom_label_repel(aes(shape, rate, label = row.names(gamsum)), alpha = 0.3) +
  geom_label_repel(aes(shape, rate, label = row.names(gamsum)), alpha = 1, fill = NA) +
  xlab(TeX('Mean of $\\pi(\\alpha$ | D)')) +
  ylab(TeX('Mean of $\\pi(\\beta$ | D)'))
     
# Lineage distances ----
load('../../02_dist_stats/data/lineage_dat.RData')

MLgam <- by(ldat$ndist, ldat$node, fitdist, distr = 'gamma', method = 'mle')
MLnor <- by(ldat$ndist, ldat$node, fitdist, distr = 'norm', method = 'mle')

pdf('../outputs/ML/lineage_diagnostics.pdf', width = 10, height = 3, onefile = TRUE)
for (node in names(MLgam)) {
  plot_diag(MLgam[[node]], MLnor[[node]], node)
}
dev.off()

gamsum <- t(data.frame(lapply(MLgam, summary_ML), check.names = FALSE))
norsum <- t(data.frame(lapply(MLnor, summary_ML), check.names = FALSE))

gamsum <- data.frame(gamsum, dist = 'Gamma')
norsum <- data.frame(norsum, dist = 'Normal')

write.csv(file = '../outputs/ML/lineage_gamsum.csv', gamsum)
write.csv(file = '../outputs/ML/lineage_norsum.csv', norsum)
