#!/usr/bin/env Rscript

source('Rlib/STRTprepHelper.R', chdir=T)
helper <- STRTprepHelper$new(name='pca', required_packages=c('RColorBrewer'))

library(RColorBrewer)

pca_by_spearman <- function(dat) {
    nr <- nrow(dat)
    nc <- ncol(dat)
    vname <- colnames(dat)
    heikin <- colMeans(dat)
    bunsan <- apply(dat, 2, var)
    sd <- sqrt(bunsan)
    r <- cor(dat, method='spearman')
    result <- eigen(r)
    eva <- result$values
    eve <- result$vectors
    contr <- eva/nc*100
    cum.contr <- cumsum(contr)
    fl <- t(sqrt(eva)*t(eve))
    fs <- scale(dat)%*%eve*sqrt(nr/(nr-1))
    names(heikin) <- names(bunsan) <- names(sd) <- rownames(r) <- colnames(r) <- rownames(fl) <- colnames(dat)
    names(eva) <- names(contr) <- names(cum.contr) <- colnames(fl) <- colnames(fs) <- paste("PC", 1:nc, sep="")
    return(structure(list(mean=heikin, variance=bunsan, standard.deviation=sd, r=r, factor.loadings=fl, eva=eva, contribution=contr, cum.contribution=cum.contr, factor.scores=fs), class="pca"))
}

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
  pdf(sprintf('%s.pdf', helper$output_prefix), pointsize=6)
  ppar <- par()
  width_label <- max(strwidth(levels(point_indexes), units='inches'),
                     strwidth(levels(color_indexes), units='inches'))
  dev.off()
  pca <- pca_by_spearman(t(nreads))
  pdf(sprintf('%s.pdf', helper$output_prefix), pointsize=6,
      width=1*components+6*ppar$csi+width_label,
      height=1*components)
  par(mar=c(0, 0, 0, 0), mgp=c(2.5, 1, 0))
  scores <- pca$factor.score
  colnames(scores) <-
    sprintf("%s (%2.1f%%)", colnames(scores), pca$contribution)
  pairs(scores[, 1:components],
        col=color_palette[color_indexes],
        pch=point_palette[point_indexes],
        oma=c(2, 2, 2, 2+6+width_label/ppar$csi), gap=1)
  if(is.null(options$POINT)) {
    legend(1, 1, xpd=T, bty='n', ncol=1, xjust=1, yjust=1, x.intersp=1,
           legend=levels(color_indexes),
           pch=rep(1, nlevels(color_indexes)),
           col=color_palette[1:nlevels(color_indexes)])
  } else {
    legend(1, 1, xpd=T, bty='n', ncol=1, xjust=1, yjust=1, x.intersp=1,
           legend=c(levels(color_indexes), levels(point_indexes)),
           pch=c(rep(1, nlevels(color_indexes)), point_palette),
           col=c(color_palette[1:nlevels(color_indexes)],
                 rep('black', nlevels(point_indexes))))
  }
  message(width_label)
  dev.off()
  save(pca, file=sprintf('%s.RData', helper$output_prefix), compress='gzip')
} else
    warning('Skipped - significant features were less than 3.')


