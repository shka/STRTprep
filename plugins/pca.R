#!/usr/bin/env Rscript

source('Rlib/STRTprepHelper.R', chdir=T)
helper <- STRTprepHelper$new(name='pca', required_packages=c('RColorBrewer'))

library(RColorBrewer)

annotations <- helper$samples$annotations
nreads <- helper$expressions$fluctuated$normalized_levels
options <- helper$options
components <- options$COMPONENTS

if(is.null(options$COLOR)) {
  color_indexes <- factor(rep(1, nrow(annotations)))
  color_palette <- c('black')
} else {
  color_indexes <- factor(annotations[, options$COLOR])
  color_palette <- c(brewer.pal(9, 'Set1'), brewer.pal(7, 'Set2'))
}

if(is.null(options$POINT)) {
  point_indexes <- factor(rep(1, nrow(annotations)))
  point_palette <- c(1)
} else {
  point_indexes <- factor(annotations[, options$POINT])
  point_palette <- 1:nlevels(point_indexes)
}

if(nrow(nreads) > 3) {
  pca <- prcomp(t(nreads), scale=T)
  pca$loadings <- t(pca$sdev*t(pca$rotation))
  pca$proportions <- pca$sdev^2/nrow(nreads)
  pca$scores <- pca$x*sqrt(ncol(nreads)/(ncol(nreads)-1))
  save(pca, file=sprintf('%s.RData', helper$output_prefix), compress='gzip')
                                        #
  pdf(sprintf('%s.pdf', helper$output_prefix), pointsize=6,
      width=1.115*components, height=1.115*components)
  scores <- pca$scores
  colnames(scores) <-
    sprintf("%s (%2.1f%%)", colnames(scores), pca$proportions*100)
  par(mar=c(0, 0, 0, 0), mgp=c(2.5, 1, 0))
  pairs(scores[, 1:components],
        col=color_palette[color_indexes],
        pch=point_palette[point_indexes],
        gap=1)
  dev.off()
                                        #
  pdf(sprintf('%s.legend.pdf', helper$output_prefix), pointsize=6,
      width=1.115, height=1.115)
  plot(1, 1, axes=F, pch=NA, xlab='', ylab='')
  if(is.null(options$POINT)) {
    legend('topleft', xpd=T, bty='n', ncol=1,
           legend=levels(color_indexes),
           pch=rep(1, nlevels(color_indexes)),
           col=color_palette[1:nlevels(color_indexes)])
  } else {
    legend('topleft', xpd=T, bty='n', ncol=1,
           legend=c(levels(color_indexes), levels(point_indexes)),
           pch=c(rep(1, nlevels(color_indexes)), point_palette),
           col=c(color_palette[1:nlevels(color_indexes)],
                 rep('black', nlevels(point_indexes))))
  }
  dev.off()
} else {
  warning('Skipped - significant features were less than 3.')
  unlink(sprintf('%s.%s', helper$output_prefix, c('RData', 'pdf', 'legend.pdf')))
}


