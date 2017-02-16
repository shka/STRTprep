#!/usr/bin/env Rscript

source('Rlib/STRTprepHelper.R', chdir=T)
helper <- STRTprepHelper$newPlugin(
  name='heatmap_diffexp',
  required_packages=c('renozao/pkgmaker@develop', 'renozao/NMF'))
options <- helper$options
samples <- helper$samples

if(helper$comparison != 'global') {
  annotations <- samples$annotations[, union(options$ANNOTATIONS, 'CLASS')]
  if(is.null(options$LABELS)) {
    annotations[, 'CLASS'] <- samples$annotations[, 'CLASS']
  } else {
    annotations[, 'CLASS'] <-
      factor(samples$annotations[, 'CLASS'], labels=options$LABELS)
  }
} else { annotations <- samples$annotations[, options$ANNOTATIONS] }

if(is.null(options$COLUMN_ORDER_BY)) {
    colv <- TRUE
} else {
    colv <- order(as.vector(samples$annotations[, options$COLUMN_ORDER_BY]))
}

solarized <- c('#dc322f', '#859900', '#268bd2', '#b58900', '#cb4b16', '#6c71c4', '#d33682', '#2aa198', '#002b36')

colors <- list()
for (i in 1:ncol(annotations)) {
  col <- colnames(annotations)[i]
  if (is.factor(annotations[, col])) {
    annotations[, col] <- factor(as.character(annotations[, col]))
    colors[[col]] <- rep(solarized, length=length(levels(annotations[, col])))
  } else
    colors[[col]] <- rep(solarized, length=ncol(annotations))[i]
}

nreads <- helper$expressions$offset$significant$normalized_levels

distfun <- function(x) as.dist((1-cor(t(x), method='spearman'))/2)
clustfun <- function(d) hclust(d, method='ward.D2')

if(nrow(nreads) > 3) {
  library(NMF)
  nmf.options(grid.patch=T)
  heatmap <- aheatmap(
    log10(nreads), scale='row', Colv=colv,
    hclustfun=clustfun, distfun=distfun,
    annCol=annotations, annColors=colors,
    fontsize=6, layout='dlmL|dalm',
    filename=sprintf('%s.pdf', helper$output_prefix))
  save(heatmap,
       file=sprintf('%s.RData', helper$output_prefix), compress='gzip')
} else
  warning('Skipped - significant features were less than 3.')

###

### Local Variables:
### eval: (setq tmppath (file-name-directory (buffer-file-name)))
### eval: (add-to-list 'exec-path (concat tmppath "../.homebrew/bin"))
### eval: (setq default-directory (concat tmppath ".."))
### End:
